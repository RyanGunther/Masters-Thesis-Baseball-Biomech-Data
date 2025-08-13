# 1.Bind swings (R or L) to match the template swing's # of columns and column names
# 2.Compute global principal markers across all swings
# 3.Segment each swing locally (adding a cut at beginning and end of each swing)
# 4.Global divisive clustering will determine the # of unique behaviour types and determine the behaviour pattern of each swing
# 5. LASSO via glmnet to determine if these fv's have predictive value


source("dynamic_PPCA_pipeline.R")

# have used the file shown below as template extensively, am confident with its shape and data, hence using it as the template
params <- list(
  swing_dir = "/Users/ryangunther/Downloads/ALL_OBP_FILES/content/ALL_CSV_FILES",
  template = "000080_000282_71_188_R_010_1070.csv",
  handed = "R",
  var_thresh = 0.95,
  seed = 20250728,
  tol_global = 0.4
)

# some csvs are missing columns (often bat markers)
# some csvs have col names Marker1_X while others have MARKER1_X
# some csvs have random 0s found in the middle of the data. I am not overly familiar with the type of mocap markers used but I think it could be a calibration or charging issue

# match all csv cols to the template cols
read_as_template <- function(path, tmpl_hdr, tmpl_lc) {
  df <- tryCatch(read.csv(path, check.names = FALSE),
                 error = function(e) return(NULL))
  if (is.null(df)) return(NULL)
  hdr_lc <- tolower(names(df))
  if (!setequal(hdr_lc, tmpl_lc)) return(NULL)
  names(df) <- tmpl_hdr[match(hdr_lc, tmpl_lc)]
  df[ , tmpl_hdr, drop = FALSE]
}

# bind all swings together
# this will take our individual 600-1000 row, 165 col csvs and bind/stack them all on top of each other creating a 165 col, massive length row csv
bind_swings <- function(dir, handed, template, verbose = TRUE) {
  tmpl_hdr <- names(read.csv(template, nrows = 1, check.names = FALSE))
  tmpl_lc  <- tolower(tmpl_hdr)
  pat <- sprintf("_[%s%s]_.*\\.csv$", handed, tolower(handed))
  files <- list.files(dir, pattern = pat, full.names = TRUE)
  
  good_mats <- vector("list", length(files))
  keep_idx <- logical(length(files))
  bad_info <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    df <- read_as_template(files[i], tmpl_hdr, tmpl_lc)
    if (is.null(df)) {
      hdr  <- tolower(names(read.csv(files[i], nrows = 1, check.names = FALSE)))
      miss <- setdiff(tmpl_lc, hdr)
      bad_info[[i]] <- data.frame(
        file = basename(files[i]),
        n_missing = length(miss),
        missing_cols = paste(tmpl_hdr[match(miss, tmpl_lc)], collapse = ", "),
        stringsAsFactors = FALSE
      )
    } else {
      good_mats[[i]] <- as.matrix(df)
      keep_idx[i]    <- TRUE
    }
  }
  
  list(
    big_mat = if (any(keep_idx)) do.call(rbind, good_mats[keep_idx]) else NULL,
    good_mats = good_mats[keep_idx],
    good_files = files[keep_idx],
    bad_files = do.call(rbind, bad_info[!keep_idx])
  )
}

# column name remap while ignoring case
apply_template_header <- function(df, tmpl_hdr, tmpl_lc) {
  hdr_lc <- tolower(names(df))
  names(df) <- tmpl_hdr[ match(hdr_lc, tmpl_lc) ]
  df
}


swing_template <- file.path(params$swing_dir, params$template)
tmpl_hdr <- names(read.csv(swing_template, nrows = 1, check.names = FALSE))
tmpl_lc <- tolower(tmpl_hdr)
bind_out <- bind_swings(params$swing_dir, params$handed, swing_template)

if (is.null(bind_out$big_mat))
  stop("No swings passed the template‑matching step – nothing to process.")

cat("Bound swings:", length(bind_out$good_files), "\n")
# all Rhand swings that have full data and all columns bound together properly

# find global PMS
# this may take some additional context to understand, looking at DimensionalityReductionBiomechData file on my github will explain further
# in short, clustering and PFA determines the most valuable clusters to retain that explain enough % of the variance in the data

set.seed(params$seed)
global_pm <- select_principal_markers(bind_out$big_mat,
                                      var_thresh = params$var_thresh,
                                      seed = params$seed)
cat("Global PMs:", paste(global_pm, collapse = ", "), "\n")
# global PMs across all Rhand swings are just 6 (i.e. 18 dimensions):
# Marker10, Marker7, LMKNE, RWRB, RTHI, T10 (X,Y,Z dim for each one)

# LOCAL segmenting
# now want to segment the swings by using the global PMs and determining the covariance matrix of each one with each other within groups of frames of a swing
# more context for this also found in the DRBD.R file on github
# essentially there is an algorithm that checks a small window of frames/rows in a swing, learns the cov structure, then checks against the next 25 rows to see if it's similar within a certain threshold
# if similar, extend the learning window and try again. If different enough, cut the data and start a new segment within that swing

seg_list <- mapply(
  process_swing_fixedpm,
  input = bind_out$good_files,
  MoreArgs = list(
    pm = global_pm,
    tmpl_hdr = tmpl_hdr,
    tmpl_lc  = tmpl_lc
  ),
  SIMPLIFY = FALSE
)

names(seg_list) <- basename(bind_out$good_files)
seg_list
# making progress here, but need to double check the naming conventions
# naming conventions good

pool_segments <- function(seg_list) {
  do.call(rbind, lapply(seq_along(seg_list), function(i) {
    nseg <- nrow(seg_list[[i]]$feature_matrix)
    data.frame(
      swing_id = i,
      seg_id = seq_len(nseg),
      swing_file = names(seg_list)[i],
      seg_list[[i]]$feature_matrix,
      row.names  = NULL
    )
  }))
}

seg_pool <- pool_segments(seg_list)
View(seg_pool)
fv_mat   <- as.matrix(seg_pool[ , -(1:3)])
cat("Total segments pooled:", nrow(fv_mat), "\n")
# 865 segments in existence
# this will be reduced later since swings that follow a X-Y-X pattern would have 3 segs but reduce down to just 2 practical segs
# also, swings can be split into 2 segs, but have similar enough feat vecs to be clustered into the same behaviour type, result in Seg1-Seg2 type swings being pushed back together into 1-1 (aka 1) behaviour type swing

# experiment with different tols to get a feel for how many clusters
# roughly 10 is ideal, first pass gave 200+ different cluster IDs aka very little clustering at all
# need a more feasible number that actually finds like-valued feat vecs 
global_divisive <- function(X, tol) {
  ids <- seq_len(nrow(X))
  out <- divisive_split_idx(X, tol = tol, ids = ids)
  lab <- integer(nrow(X))
  for (cid in seq_along(out$groups))
    lab[ out$groups[[cid]] ] <- cid
  list(cluster_id = lab, split_tree = out$tree_lines)
}

# explore cluster counts quickly
diag_nclusters <- function(X, cand = c(0.1,0.2,.3,.4,.5,1,2)) {
  sapply(cand, function(t) length(unique(global_divisive(X, tol = t)$cluster_id)))
}

cat("Clusters per tol:\n")
print(diag_nclusters(fv_mat))
# happy with 9 clusters hence tol of 0.4


glob_cl <- global_divisive(fv_mat, tol = 0.4)

#fit the columns in my desired order
seg_pool$cluster_id <- glob_cl$cluster_id
seg_pool <- seg_pool[ c(1:3, ncol(seg_pool), 4:(ncol(seg_pool) - 1)) ]
View(seg_pool)



# list of vectors of frame indices grouped by cluster ID
frames_by_cluster <- function(sw) {
  cuts <- sw$breaks_forward
  cid <- sw$cluster_global
  nseg <- length(cid)
  
  split(seq_len(cuts[nseg + 1] - 1),
        rep(cid, diff(cuts)))
}


# Use the same mean+cov construction as what was used in the building of the feature vectors initially
# however, in the building of the design mat, I am gonna leave out the means hence the toggle feature in the function's signature
# same math as before, but allow turning off the means
cluster_fv <- function(dat, frames, pm, tau = 0.75, w = 0.8, include_means = FALSE) {
  xyz_cols <- unlist(lapply(pm, \(m) paste0(m, c("_X","_Y","_Z"))))
  X <- as.matrix(dat[frames, xyz_cols, drop = FALSE])
  m <- compute_PPCA_C(X, tau)

  # vectorize covariance (upper triangle incl. diag)
  v_cov <- m$C[upper.tri(m$C, diag = TRUE)]

  # if wanting to still include the means, use previous weights of 0.8 for covs and 0.2 for means
  if (include_means) {
    c(w * v_cov, (1 - w) * m$mu)
  } else {
    # if not including means, no need for linear weights scaling the data
    v_cov
  }
}



# build one row per swing x global behaviour i.e. if a swing follows a 1-2-1 pattern,
# it would now have one row for all the 1 type behaviours which hypothetically be frames 1-450 and 600-900

pool_cluster_level <- function(seg_list, include_means = FALSE) {
  rows <- lapply(seq_along(seg_list), function(i) {
    sw <- seg_list[[i]]
    # frames for each global cluster visited by this swing
    cuts <- sw$breaks_forward
    cid  <- sw$cluster_global
    fr_by_c <- split(seq_len(cuts[length(cuts)] - 1), rep(cid, diff(cuts)))
    
    do.call(rbind, lapply(names(fr_by_c), function(cid_k) {
      fv <- cluster_fv(sw$reconstructed,
                       frames = fr_by_c[[cid_k]],
                       pm     = sw$principal_markers,
                       include_means = include_means)
      data.frame(
        swing_id   = i,
        cluster_id = as.integer(cid_k),
        swing_file = names(seg_list)[i],
        t(fv), row.names = NULL, check.names = FALSE
      )
    }))
  })
  # stack all swings together in one big df
  do.call(rbind, rows)
}

# now build the feat vecs on a per swing basis
# important: since glmnet would struggle with different length predictors, and since
# comparing a swing with beh type 1-5 is not apples2apples with a swing of 2-9, i'm going to need to entirely flesh out the fv's
# hence, each feat vec is gonna be length 9x171
# 9 is the # of clusters
# 171 is the upper tri of the 18x18 cov mat

# therefore, if a swing follows a 1-5 swing pattern:
# it will have 171 nonzero values, followed by (3x171) 0s, 171 more nonzero vals, then (4x171) 0s
# the 2-9 swing would have 0s, nonzeros, 0s, 0s, 0s, 0s, 0s, 0s, nonzeros

make_block_X_by_swing <- function(cluster_pool_cov,
                                  key_col     = "swing_file",
                                  cluster_col = "cluster_id") {
  # cov-only: these will be just V* columns
  fv_cols <- grep("^(V)?[0-9]+$", names(cluster_pool_cov), value = TRUE)
  levels  <- sort(unique(cluster_pool_cov[[cluster_col]]))  # actual cluster labels (e.g., 1..9)
  K <- length(levels); p <- length(fv_cols)
  
  swings <- unique(cluster_pool_cov[[key_col]])
  X <- matrix(0, nrow = length(swings), ncol = K * p)
  rownames(X) <- basename(swings)
  colnames(X) <- as.vector(outer(paste0("cid", levels, "_"), fv_cols, paste0))
  
  for (i in seq_along(swings)) {
    sub <- cluster_pool_cov[ cluster_pool_cov[[key_col]] == swings[i], , drop = FALSE ]
    if (!nrow(sub)) next
    for (r in seq_len(nrow(sub))) {
      j <- match(sub[[cluster_col]][r], levels)         # which cluster block
      idx <- ((j - 1) * p + 1):(j * p)
      X[i, idx] <- as.numeric(sub[r, fv_cols])
    }
  }
  X
}


for (i in seq_len(nrow(seg_pool))) {
  s <- seg_pool$swing_id[i]
  k <- seg_pool$seg_id[i]
  seg_list[[s]]$cluster_global[k] <- seg_pool$cluster_id[i]
}

cluster_pool <- pool_cluster_level(seg_list)
View(cluster_pool)
# glmnet lasso regression




# first need to join with metadata
# build one-row-per-swing design (cov-only, block by cluster)
cluster_pool_cov <- cluster_pool   # (clarity)
View(cluster_pool_cov)

# remember that the fv's in the below matrix will look different from those in seg_pool by a factor of 0.8!
# in seg pool there are scaled down to make way for the 0.2 * mean vals. Those won't exist here since we're keeping only the covs

X_bySwing <- make_block_X_by_swing(cluster_pool_cov)
View(X_bySwing)

# outcomes & folds keyed by swing (basename)
source("build_df_full.R")
df_full <- build_df_full(
  swing_dir = params$swing_dir,
  handed    = params$handed,
  meta_path = "/Users/ryangunther/Downloads/metadata(hitting).csv"
)

# need to intelligently join feat vecs/ files to correct exit velo / bat speed and user values

df_full$key <- basename(df_full$file)

exit_lookup <- setNames(df_full$exit_velo_mph_x,     df_full$key)
bat_lookup  <- setNames(df_full$bat_speed_mph_max_x, df_full$key)
user_lookup <- setNames(df_full$user,                df_full$key)

swing_keys <- rownames(X_bySwing)
head(swing_keys)
y_ev  <- exit_lookup[swing_keys]
y_bs  <- bat_lookup[swing_keys]
user  <- user_lookup[swing_keys]

ok <- complete.cases(y_ev, y_bs, user)
ok
X_bySwing <- X_bySwing[ok, , drop = FALSE]
y_ev <- y_ev[ok]; y_bs <- y_bs[ok]; user <- user[ok]
# don't want to let users bleed across training/testing sets. hence maintaining all swings from a unique user into one fold
foldid <- as.integer(factor(user))


# scale-only (no centering) so zeros stay zero (could allow centering as it's functionally the same, but it's visual much more difficult to see which behaviours a feat vec contains)
sds <- apply(X_bySwing, 2, sd, na.rm = TRUE)
sds[!is.finite(sds) | sds == 0] <- 1    # guard

X_scale_only <- sweep(X_bySwing, 2, sds, "/")
View(X_scale_only)
set.seed(2025)

cv_ev <- cv.glmnet(X_scale_only, y_ev, alpha = 1,
                   foldid = foldid, standardize = FALSE)
cv_bs <- cv.glmnet(X_scale_only, y_bs, alpha = 1,
                   foldid = foldid, standardize = FALSE)

cat("Exit-velo RMSE:", round(sqrt(min(cv_ev$cvm)), 2), "\n")
cat("Bat-speed RMSE:", round(sqrt(min(cv_bs$cvm)), 2), "\n")

cv_ev
cv_bs

# some nice predictive value found in the feature vectors to predict bat speed
# a good amount better than the baseline
# fv's not useful for predicting exit velo but there are other factors at play like contact quality (as in topping/bottoming the ball) that wouldn't be accounted for here

# biggest takeaway is that there are more efficient and less efficient ways to move, which are captured in the covariance structure of a dimensionality-reduced dataset!






