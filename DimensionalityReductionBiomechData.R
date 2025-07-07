# motion capture data segmentation technique
# the goal of this code is to take markered mocap data provided by Driveline's open biomechanics project and complete dimensionality reduction on it
# this means to take the lengthy 810x165 time-series dataset (that's only one swing of approx 800!) and extract as much information as possible from it, while also reducing it as much as possible
# I want to do this so the data can be more easily understood (which features are most important?) and for ease of modelling
# most mocap data is heavily redundant (there are 10 markers on the bat alone) so doing this cleanly is of massive importance

# The techniques used in this code are based off of Liu et al, 2006, with more in-depth descriptions of the methods provided by Barbic et al, 2004
# the data can be found here: https://github.com/drivelineresearch/openbiomechanics/tree/main/baseball_hitting

# An overview of the method:
# First, I will perform principal feature analysis (PFA) on the full dataset, by completing a Principal Component Analysis (PCA), then performing a k means clustering wherein the features
# that eventually get selected appropriately cover all clusters used, while having the most coverage and least avg distance to the cluster centers

# testing fully through with a single swing first
hit080_10 <- read.csv("/Users/ryangunther/Downloads/000080_282.csv")
# 080 is the hitter's unique player ID
# 10 is swing number ten in the session
View(hit080_10)
str(hit080_10)

### (added at the end. first pass was only on one swing) tie multiple swings from the same hitter together and determine if the funcs still work well ----
hit080_01 <- read.csv("/Users/ryangunther/Downloads/hit080_01.csv")
hit080_02 <- read.csv("/Users/ryangunther/Downloads/hit080_02.csv")
hit080_04 <- read.csv("/Users/ryangunther/Downloads/hit080_04.csv")
hit080_06 <- read.csv("/Users/ryangunther/Downloads/hit080_06.csv")
concat_df_080 <- rbind(hit080_10,hit080_01, hit080_02, hit080_04, hit080_06)
View(concat_df_080)



# important to note that the global "0" points in the dataset are: ----
# z: the ground
# y: the edge of home plate closest to the RH batter's box
# x: the front of home plate
# time series data points are in metres from those 0 points
# pos x is towards pitcher's mound
# pos y is towards RHB box
# pos z is up towards the sky



# perform PFA on the dataset ----
# originally used only hit080_10 in the next line, changed to concat_df_080 (matrix of multiple rbind'ed swings) later on
Y  <- as.matrix(concat_df_080)
Yc <- scale(Y, center=TRUE, scale=FALSE)
View(Yc)
covY <- cov(Yc)
covY
# below should be the principal components
eig <- eigen(covY)
# below should be the var explained by each principal component
lambda <- eig$values
lambda # they're already ranked in order of most variance explained
V <- eig$vectors # principal directions that the data will be projected/rotated onto

# find the ones that explain the most var
cumvar <- cumsum(lambda) / sum(lambda) #calculate cum var explained
cumvar
q <- which(cumvar >= .99)[1]
q
Aq <- V[, 1:q] # select the first q PCs
Aq            # 165 x q mat

rownames(Aq) <- colnames(Y)

W <- abs(Aq)
View(W)
# (12 clusters on the concat matrix since 10 PCs):
K <- q + 2

# now begins the k means portion where I am determining which components have the most coverage and least distance from the cluster centers
set.seed(42)
km <- kmeans(W, centers = K, nstart = 150)
km
clust <- km$cluster
centers <- km$centers


feat_names <- rownames(Aq)
feat_mark <- sub("_[XYZ]$", "", feat_names)


# distance to centroid of each cluster (used as tiebreaker)
dists <- rowSums((W - centers[clust, ])^2)

# prepare empty DF
markers <- unique(feat_mark)
marker_info <- data.frame(
  marker = markers,
  cover_count = integer(length(markers)),
  total_dist = numeric(length(markers)),
  stringsAsFactors = FALSE
)

# fill in DF with clusters participated in and total distance to their centroids, then work backwards deleting the least important row of this df until I reach a minimum dataset
for (i in seq_along(markers)) {
  m <- markers[i]
  idxs <- which(feat_mark == m)
  cls <- unique(clust[idxs])
  marker_info$cover_count[i] <- length(cls)
  marker_info$total_dist[i] <- sum(dists[idxs])
}
View(marker_info)
# now have a table showing how many clusters a marker appears in (max 3 because 3 dimensions (x,y,z)) and the total dist to centroid
# sorting this provides the most useful markers per the text:

# "a marker that appears in more distinct clusters is considered to be more
# important. To break ties between markers, we prefer those whose
# sum of the square distances from the marker’s features to their
# cluster mean is minimal."



# sort by least important
marker_info <- marker_info[order(marker_info$cover_count, -marker_info$total_dist), ]
marker_info



# Q: why does it drop markers that are "more useful" than the final kept ones?
# A: because each cluster needs at least one representative
all_clusters <- sort(unique(clust))
keep <- marker_info$marker

for (m in marker_info$marker) {
  cand <- setdiff(keep, m)
  # below: check if each cluster still has >= 1 rep
  covered <- sort(unique(unlist(
    lapply(cand, function(x) unique(clust[feat_mark == x]))
  )))
  # if all still covered, candidate list becomes new keep list, continue loop
  if (identical(covered, all_clusters)) keep <- cand
}

principal_markers <- keep
principal_markers
str(principal_markers)
# goal is stated as:
# "We continue removing the least important markers as
# long as every cluster is still covered. This process is repeated until
# no more markers can be eliminated. "

# Hence the list of just 6 markers
# and the selected 6 not aligning with the top of marker_info leaderboard

# "Marker5" "RWRA"    "RMKNE"   "LANK"    "LASI"    "Marker4" are the 6 PFAs that cover all 12 clusters 
# in laymen's terms, these are the Right wrist, Right inner knee, Left ankle, front of Left hip, and two markers on the bat
# this link has a photo of every marker on the body with their labels: https://github.com/drivelineresearch/openbiomechanics/tree/main/baseball_hitting
# also view a different set of code on my github ("quick 2d view") to quickly visualize all the markers in the x and z dimensions with special labels for these principal markers



### Section 3.2 - Piecewise Linear Modeling ----
# attempt to replicate manual ppca method from barbic et al

# retained 6 markers x 3 dims each = 18 features
# I now have the reduced dataset. I am going to complete probabilistic PCA (PPCA) on it
# this means that of the 18 features remaining (x,y,z components of 6 principal markers), perform singular value decomposition (SVD) 
# and then again order by variance explained. Some markers (enough to get over a variance explained threshold "tau" are retained as the relevant dimensions moving forward within ppca)
# the additional dimensions (nonrelevant/"out group") are treated as isotropic ("same in all directions") Gaussian noise
# the goal here is to find the var explained "over and above" the isotropic noise level in the data
# hence why we see the line Λ_r <- pmax(λ[1:r] - σ2, 0) below, we are subtracting the noise from the relevant dimensions

xyz_cols <- unlist(lapply(principal_markers, function(m) {
  paste0(m, c("_X","_Y","_Z"))
}))
reduced_mat <- concat_df_080[, xyz_cols]
View(reduced_mat)


# ppca func to run on an already reduced dataset, in this case reduced_mat
compute_PPCA_C <- function(data_mat, tau = 0.90) {
  n <- nrow(data_mat); p <- ncol(data_mat)
  mu <- colMeans(data_mat)
  D <- data_mat - matrix(mu, n, p, byrow = TRUE)
  
  # singular value decomp portion (for more info see this link: https://www.geeksforgeeks.org/machine-learning/singular-value-decomposition-svd/)
  sv <- svd(D); Σ <- sv$d; V <- sv$v
  λ <- numeric(p); λ[1:length(Σ)] <- Σ^2 / (n - 1)
  
  # compare variance explained to tau, this generates the "in crowd" which will be the signal
  # remaining few markers are the isotropic noise
  # determine signal "over and above" mean noise by subtracting in line 193
  cumv <- cumsum(λ) / sum(λ)
  idx <- which(cumv >= tau)
  r <- if (length(idx) == 0) p else idx[1]
  
  σ2 <- if (r < p) mean(λ[(r + 1):p]) else 0
  V_r <- V[, 1:r, drop = FALSE]
  
  # guard against negative inside sqrt
  # subtract the noise
  Λ_r <- pmax(λ[1:r] - σ2, 0)
  W <- V_r %*% diag(sqrt(Λ_r), r, r)
  C <- W %*% t(W) + σ2 * diag(p)
  list(mu = mu, C = C)
}


### now "slide the window" to begin building local linear models ----
# this next section involves training an algorithm to use the function above, which determines avg covariances across some number of rows, and implement it on only a subset of the rows of the data
# (we will only be using the most reduced dataset in these calculations i.e. taking the product of ppca above)
# the model trains on the first K (50 for now) frames, and determines an avg cov matrix for those frames
# 50 frames while recording at 360 Hz is 5/36 of one second, the full dataset for an average swing by this hitter is 800 frames or roughly 2.25 seconds
# after training, the algo "looks ahead" at the next T (25 for now) frames to determine if the next set is "different enough" beyond a certain threshold to be considered a separate motion type
# if yes, the algorithm splits into segments. If not, it adds delta (5 as of now) frames to its training set, learns the avg cov mat for that, and repeats
# this nets entirely data-driven segments of the swing


N <- nrow(reduced_mat)
N
break_pts <- c()
# recall from barbic et al that Mahalanobis distance is calculated as such:
# 1/T summation t(xi-xbar) %*% solve(C) %*% (xi-xbar)

# Mahalanobis dist (H) is measuring in units of SD how far the distribution of
# K+1 to K+T changes relative to 1 to K

detect_ppca_cuts <- function(reduced_mat,
                             K     = 50,
                             T     = 25,
                             delta = 5,
                             tau   = 0.95,
                             q     = 0.99) {
  
  # have to have everything in matrix form for future matrix mult
  reduced_mat <- as.matrix(reduced_mat)
  storage.mode(reduced_mat) <- "double"
  
  N <- nrow(reduced_mat)
  p <- ncol(reduced_mat)
  
  # set chisq threshold of avg mahalanobis dist
  infl <-1000
  chi_thresh <- infl * qchisq(q, df = T * p) / T
  
  # first frame after a cut is always the start of a new segment, which now as of last updates, starts training at 0 after a cut instead of retaining the existing learned data
  cuts <- 1L   
  s <- 1L
  
  while (s + K + T - 1 <= N) {
    
    # run ppca calc on training block. starts at first 50 frames, tests on 25, if no need for a cut, advances training block by delta, tries again
    seg <- reduced_mat[s:(s+K-1), , drop = FALSE]
    ppca <- compute_PPCA_C(seg, tau)
    mu <- ppca$mu
    Cinv <- solve(ppca$C)
    
    # test on testing block of 25 frames
    probe <- reduced_mat[(s+K):(s+K+T-1), , drop = FALSE]
    diffs <- sweep(probe, 2, mu)
    H_vals <- rowSums((diffs %*% Cinv) * diffs)
    H_avg <- mean(H_vals)
    
    # cut or add delta and try again
    if (H_avg > chi_thresh) {
      cut_pt <- s + K
      cuts <- c(cuts, cut_pt)
      s <- cut_pt
    } else {
      s <- s + delta
    }
  }
  # add a cut at the very end to close the final segment
  c(cuts, N + 1L)
}

forward_breaks <- detect_ppca_cuts(reduced_mat)
forward_breaks
### IMPORTANT:
# going from one swing's dataset to the next occur at 811, 1621, 2393, 3266
# concat df had breaks at 1  561  791 1381 1606 2141 2366 3046 3246 4001 4240
# since the testing window is 25 frames, any cut that is within 25 frames of 811, 1621, 2393, 3266 is properly cutting between datasets

# every break is completely in line with expectations except one:
# 1 = beginning
# 561 = rotational portion of 1st swing
# 791 = changing from first to second swing
# 1381 = rotational portion of second swing
# 1606 = changing from second swing to third
# 2141 = rotational portion of third swing
#########  2366 - oddly enough the test portion of this window ends 2 frames before the break from third to 4th swing*** ###
######### (2366+25=2391 but break is at 2393****) ###
# 3046 = rotational portion of 4th swing
# 3246 = changing from 4th swing to 5th
# 4001 = rotational portion of 5th swing
# 4240 ending

# tricky swing in question has about 55 swings of recoil, perhaps that's the problem, need to compare to other swing recoils
# first swing recoil looks to be 43 frames
# second swing recoil looks to be 35 frames
# third swing is the mystery
# fourth swing recoil is 0-2 frames at most
# fifth and final swing recoil is 32 frames

# this might actually be a correctly segmented break going the opposite direction (recoil) than the normal swing




### sliding the window backwards ----
# sometimes the motion at the beginning of the sequence is very stable and the end is very full of movement
# hence to avoid training on edge cases run the same algorithm but in reverse
rev_breaks <- c()

# reverse version of window slide to make cuts:

detect_ppca_cuts_backward <- function(reduced_mat,
                                      K      = 50,   # training length
                                      T      = 25,   # testing length
                                      delta  = 5,
                                      tau    = 0.95,
                                      q      = 0.99,
                                      infl   = 1000)
{
  # ensure matrix form
  reduced_mat <- as.matrix(concat_df_080)
  storage.mode(reduced_mat) <- "double"
  
  N <- nrow(reduced_mat)
  p <- ncol(reduced_mat)
  
  # chi sq thresh for avg mahalanobis dist
  chi_thresh <- infl * qchisq(q, df = T * p) / T
  
  # beginning point of reverse which is really the end point of the data
  cuts <- N + 1L      
  e <- N
  
  while (e - K - T + 1 >= 1) {
    
    # ppca on the last K frames of the current segment
    seg_idx <- (e - K + 1):e
    seg <- reduced_mat[seg_idx, , drop = FALSE]
    ppca <- compute_PPCA_C(seg, tau)
    mu <- ppca$mu
    Cinv <- solve(ppca$C)
    
    # testing forward the next T frames
    probe_idx <- (e - K - T + 1):(e - K)
    probe <- reduced_mat[probe_idx, , drop = FALSE]
    diffs <- sweep(probe, 2, mu)
    H_vals <- rowSums((diffs %*% Cinv) * diffs)
    H_avg <- mean(H_vals)
    
    # determine segment break or add delta and re-iterate
    if (H_avg > chi_thresh) {
      cut_pt <- e - K
      cuts <- c(cuts, cut_pt)
      e <- cut_pt
    } else {
      e <- e - delta
    }
  }
  # add frame number 1 to the cuts to close the last (first in real terms) segment
  c(cuts, 1L)
}


rev_breaks <- detect_ppca_cuts_backward()
rev_breaks

# the reverse window sliding function cuts at a subset of all the same spots as the forward version
# numbers won't align perfectly since forward cut will happen 1-25 frames BEFORE the true break, and reverse will happen 1-25 frames AFTER

### commenting this part out, since the cuts aren't at identical spots, they're on opposite sides of the true break, but do show the same behaviour
# reverse is essentially functioning as a sanity check
## combine forward and backward cuts, de-duplicate, sort
# all_breaks <- sort(unique(c(1, break_pts, rev_break)))
# adding a few as a test run
# all_breaks <- sort(unique(c(461,711,251, all_breaks)))
# all_breaks


### for tomorrow: ----

### June 13 early morning update
# ppca func is now splitting properly
# need to give more thought to tau
# K, T, delta in a good place right now
# T is pretty small so it can pick up miniscule motions (the whole swing is usually about 40 frames so this is necessary)


# segmenting, feature vectors, PCA, k means clustering ----

# determine largest eucl. dist of points in a cluster to its centroid

# writing the code so that the feature vectors enter the divisive splitting algorithm, if their dist from means is greater than tol, run k means clustering with k=2
# this will separate the fv's into smaller bundles to then be tested again to see if they are similar enough to become a leaf or need to split further
# this concept would make more sense if the motion data repeating similar motions over and over, in practice, in a baseball swing, very few segments will get bundled unless tau for ppca is tuned very low and we have many segments of similar motions
# if <= 2 FV's in the current bundle, the below func just computes euclid. dist between the scaled feature vectors and determines if they are within or beyond the tolerance level
divisive_split_idx <- function(x, tol, ids, depth = 0) {
  x <- as.matrix(x)
  n <- nrow(x)
  cen <- colMeans(x)
  
  # 1-point leaf: just make a leaf
  if (n == 1) {
    spread <- 0
    cat(strrep("  ", depth),
        "depth=", depth, "  n=1  spread=", sprintf("%.3f\n", spread), sep = "")
    cat(strrep("  ", depth+1), "→ leaf\n", sep = "")
    return(list(ids))
  }
  
  # 2-point: compare euclid dist from one fv to mean of both, if dist > tol, split into 2 separate 1 FV leaves, if < tol, 2 FV's in same leaf
  if (n == 2) {
    spread <- sqrt(sum((x[1,] - cen)^2))
    cat(strrep("  ", depth),
        "depth=", depth, "  n=2  spread=", sprintf("%.3f\n", spread), sep = "")
    
    if (spread > tol) {
      # split into two singles
      return(list(ids[1], ids[2]))
    } else {
      cat(strrep("  ", depth+1),
          sprintf("→ 2-row leaf (spread = %.3f)\n", spread), sep = "")
      return(list(ids))
    }
  }
  
  # general case: n > 2
  # compare dist from furthest fv to mean of all within that bundle, if > tol, run kmeans clustering with nstart=100, split into bundles along kmeans clustering lines
  dists <- sqrt(rowSums((x - matrix(cen, n, ncol(x), byrow=TRUE))^2))
  spread <- max(dists)
  cat(strrep("  ", depth),
      "depth=", depth, "  n=", n, "  spread=", sprintf("%.3f\n", spread), "\n", sep = "")
  
  if (spread > tol) {
    km <- kmeans(x, centers = 2, nstart = 100)
    out <- list()
    for (cl in 1:2) {
      sel <- which(km$cluster == cl)
      out <- c(
        out,
        divisive_split_idx(x[sel, , drop=FALSE],
                           tol,
                           ids[sel],
                           depth + 1)
      )
    }
    return(out)
  } else {
    cat(strrep("  ", depth+1), "→ leaf\n", sep = "")
    return(list(ids))
  }
}


# determine avg cov and mu for each segment, assemble feature vecs
# (feature vectors will be # of segments columns x 189 rows) because the covariance matrices are 18x18
# the 189 rows are composed of 171 variances or covariances, pulled from the upper triang of each cov mat, and 18 column means 171+18=189
seg_models <- lapply(seq_len(length(forward_breaks) - 1), function(j) {
  idx <- forward_breaks[j]:(forward_breaks[j + 1] - 1)
  compute_PPCA_C(reduced_mat[idx, ], tau = 0.75)
})

# assemble feature vector as per f = [wvT, (1-w)µT]T
# "difference in distribution between the segments is mainly due to the covariance matrix."
# "Mean vector has very little impact"
w <- 0.8
p <- ncol(reduced_mat)



feature_mat <- t(vapply(seg_models, function(m) {
  v <- m$C[upper.tri(m$C, diag = TRUE)]
  # total length p(p+3)/2 = 189
  c(w * v, (1 - w) * m$mu)
}, numeric(p * (p + 3) / 2)))
# calling feature_mat should emit (len(forward_breaks)-1) rows of 189 col vectors
# first 171 cols are all the variances and covariances with the segment
# final 18 figures are the within segment means of each feature
str(feature_mat)

# if <=2 segments, skip kmeans clustering, since they don't add any value
# if 3+ segments, complete PCA and kmeans clustering:

set.seed(1)
tol <- 10

# center and scale, apply PCA
pca_out <- prcomp(feature_mat, center = TRUE, scale. = TRUE)

var_explained <- cumsum(pca_out$sdev^2) / sum(pca_out$sdev^2)
var_explained
k95 <- which(var_explained >= 0.95)[1]
k95
cat("Retaining", k95, "PCs (", round(var_explained[k95]*100,1), "% of variance)\n")
reduced_mat95 <- pca_out$x[, 1:k95, drop = FALSE]
reduced_mat95
# should produce "segments" number of rows x "retained PCs" number of columns

n_seg <- nrow(reduced_mat95)
n_seg
ids <- seq_len(n_seg)
# if only 1 segment (no breaks in the motion), that piece is the only feature vector to describe the motion
# if 2 segments, calc the euclid dist and compare to tolerance level

behaviour_label <- divisive_split_idx(x = reduced_mat95, tol = tol, ids = ids, depth = 0)
print(behaviour_label)

# incredible result. Worked spectacularly, algorithm is grouping all the rotational portions of the swings together, 
# and all the "transition from end swing X - beginning swing Y" portions together
# except swing 1 which is fine
# really really successful result



### motion reconstruction ----

keep_cols <- unlist(lapply(principal_markers, function(m) {
  grep(paste0("^", m, "_"), colnames(concat_df_080))
}))
keep_cols
# 3. Build X and Z directly on the ORIGINAL (frames × vars) matrix:
X <- concat_df_080[, keep_cols, drop = FALSE]
Z <- concat_df_080[, -keep_cols, drop = FALSE]
X
str(Z)


# set each segment and frame within that segment to a cluster ID
# for segments:
seg_cluster_id <- integer(n_seg)
seg_cluster_id

for (cluster_i in seq_along(behaviour_label)) {
  segs <- behaviour_label[[cluster_i]]
  seg_cluster_id[segs] <- cluster_i
}
segs # not like the tiktok term :p
seg_cluster_id

# for frames:
n_frames <- nrow(concat_df_080)
model_id <- integer(n_frames)
model_id

for (j in seq_len(n_seg)) {
  fr <- forward_breaks[j]:(forward_breaks[j+1] - 1)
  model_id[fr] <- seg_cluster_id[j]
}
fr
model_id

# constructing B
other_cols <- setdiff(seq_len(ncol(concat_df_080)), keep_cols)
cluster_frames <- split(seq_len(nrow(concat_df_080)), model_id)
cluster_frames

Bs <- vector("list", length(cluster_frames))
mu_x_lst <- vector("list", length(cluster_frames))
mu_z_lst <- vector("list", length(cluster_frames))


for (cid in seq_along(cluster_frames)) {
  fr <- cluster_frames[[cid]]
  
  # pull cluster-specific blocks and transpose so rows = coords, cols = frames
  Xc <- t(as.matrix(X[fr, , drop = FALSE]))   # (3k) × n_f
  Zc <- t(as.matrix(Z[fr, , drop = FALSE]))   # 3(m−k) × n_f
  
  mu_x <- rowMeans(Xc);  mu_z <- rowMeans(Zc)
  Xc   <- sweep(Xc, 1, mu_x, "-")
  Zc   <- sweep(Zc, 1, mu_z, "-")
  
  # have to consider what to do if mat is non-invertible
  XtX  <- Xc %*% t(Xc)
  
  Bs[[cid]]       <- Zc %*% t(Xc) %*% solve(XtX)
  mu_x_lst[[cid]] <- mu_x
  mu_z_lst[[cid]] <- mu_z
}

# display the weights that will be used for blind reconstruction later
Bs
# looks great

# reconstruct this swing:

reconstructed <- matrix(NA_real_, nrow = nrow(concat_df_080), ncol = ncol(concat_df_080))
reconstructed


for (f in seq_len(nrow(concat_df_080))) {
  cid <- model_id[f]
  xf <- as.numeric(concat_df_080[f, keep_cols ])
  zhat <- Bs[[cid]] %*% (xf - mu_x_lst[[cid]]) + mu_z_lst[[cid]]
  
  full <- numeric(ncol(concat_df_080))
  full[ keep_cols ] <- xf
  full[ other_cols ] <- zhat
  reconstructed[f, ] <- full
}

colnames(reconstructed) <- colnames(concat_df_080)
View(reconstructed)
View(concat_df_080)
# room for improvement: get detect_ppca_cuts to cut exactly at the frame that violates tau instead of at s+K (minor change)
# investigate further the one recoil segment, although spent lots of time on it and it seems to be a legitimate motion on its own
