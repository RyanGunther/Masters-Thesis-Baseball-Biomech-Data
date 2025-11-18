# dimensionality reduction pipeline for all of OBP dataset's hitter C3D files (converted to CSV for analysis here)

`%||%` <- function(a, b) if (!is.null(a)) a else b




defaults_ctrl <- list(
  var_thresh_pfa = 0.99,
  tau_ppca_cut = 0.95,
  tau_ppca_feature = 0.75,
  K_train = 50,
  T_probe = 25,
  delta_shift = 5,
  q_chisq = 0.99,
  infl_chisq = 1000,
  tol_divisive = 5,
  w_feature = 0.8,
  seed = 42
)

update_tuning <- function(a, b) if (is.null(a)) b else a

merge_ctrl <- function(user, dflt = defaults_ctrl)
{ modifyList(dflt, user) }


# probabilistic PCA – returns mean vector and covariance matrix C
compute_PPCA_C <- function(data_mat, tau = 0.90) {
  n  <- nrow(data_mat);  p <- ncol(data_mat)
  mu <- colMeans(data_mat)
  D  <- data_mat - matrix(mu, n, p, byrow = TRUE)
  
  sv <- svd(D); Σ <- sv$d; V <- sv$v
  λ  <- numeric(p);  λ[seq_along(Σ)] <- Σ^2 / (n - 1)
  
  cumv <- cumsum(λ) / sum(λ)
  idx  <- which(cumv >= tau)
  r    <- if (length(idx) == 0) p else idx[1]
  
  σ2   <- if (r < p) mean(λ[(r + 1):p]) else 0
  V_r  <- V[, 1:r, drop = FALSE]
  Λ_r  <- pmax(λ[1:r] - σ2, 0)
  W    <- V_r %*% diag(sqrt(Λ_r), r, r)
  C    <- W %*% t(W) + σ2 * diag(p)
  list(mu = mu, C = C)
}

# principal‑feature analysis --> principal markers
select_principal_markers <- function(Y, var_thresh = 0.99, seed = 42) {
  Yc   <- scale(Y, center = TRUE, scale = FALSE)
  eig  <- eigen(cov(Yc))
  λ    <- eig$values
  cumv <- cumsum(λ) / sum(λ)
  q    <- which(cumv >= var_thresh)[1]
  Aq   <- eig$vectors[, 1:q, drop = FALSE]
  rownames(Aq) <- colnames(Y)
  
  W    <- abs(Aq)
  Kc   <- q + 2
  
  set.seed(seed)
  km   <- kmeans(W, centers = Kc, nstart = 150)
  
  clust  <- km$cluster
  centers<- km$centers
  feat_m <- sub("_[XYZ]$", "", rownames(Aq))
  
  dists <- rowSums((W - centers[clust, ])^2)
  markers <- unique(feat_m)
  info <- data.frame(marker = markers,
                     cover  = integer(length(markers)),
                     total  = numeric(length(markers)))
  for (i in seq_along(markers)) {
    m  <- markers[i]
    idx<- which(feat_m == m)
    info$cover[i] <- length(unique(clust[idx]))
    info$total[i] <- sum(dists[idx])
  }
  info <- info[order(info$cover, -info$total), ]
  
  keep <- info$marker
  allc <- sort(unique(clust))
  for (m in info$marker) {
    cand <- setdiff(keep, m)
    covered <- sort(unique(unlist(
      lapply(cand, function(x) unique(clust[feat_m == x]))
    )))
    if (identical(covered, allc)) keep <- cand
  }
  keep
}



build_global_markers <- function(dir, pattern = "_R_.*\\.csv$",
                                 var_thresh = 0.99, seed = 42) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) stop("No right‑handed swing files found.")
  
  message("Binding ", length(files), " R‑handed swings …")
  big_mat <- do.call(rbind, lapply(files, function(f)
    as.matrix(read.csv(f, check.names = FALSE))))
  
  message("Running PFA on the bound matrix …")
  pm <- select_principal_markers(big_mat,
                                 var_thresh = var_thresh,
                                 seed       = seed)
  pm
}






# forward ppca window slide
detect_ppca_cuts <- function(reduced_mat, ctrl) {
  ctrl <- merge_ctrl(ctrl)
  reduced_mat <- as.matrix(reduced_mat)
  storage.mode(reduced_mat) <- "double"
  K     <- ctrl$K_train
  Tp    <- ctrl$T_probe
  del   <- ctrl$delta_shift
  tau   <- ctrl$tau_ppca_cut
  q     <- ctrl$q_chisq
  infl  <- ctrl$infl_chisq
  
  N     <- nrow(reduced_mat)
  p     <- ncol(reduced_mat)
  chi_th<- infl * stats::qchisq(q, df = Tp * p) / Tp
  
  cuts  <- 1L
  s     <- 1L
  while (s + K + Tp - 1 <= N) {
    seg   <- reduced_mat[s:(s+K-1), , drop = FALSE]
    ppca  <- compute_PPCA_C(seg, tau)
    mu    <- ppca$mu
    Cinv  <- solve(ppca$C)
    
    probe <- reduced_mat[(s+K):(s+K+Tp-1), , drop = FALSE]
    diff  <- sweep(probe, 2, mu)
    Havg  <- mean(rowSums((diff %*% Cinv) * diff))
    
    if (Havg > chi_th) {
      cut_pt <- s + K
      cuts   <- c(cuts, cut_pt)
      s      <- cut_pt
    } else {
      s      <- s + del
    }
  }
  c(cuts, N + 1L)
}

# backward ppca sanity check
ppca_breaks_backward <- function(reduced_mat, ctrl = list()) {
  reduced_mat <- as.matrix(reduced_mat)
  storage.mode(reduced_mat) <- "double"
  ctrl <- merge_ctrl(ctrl)
  
  K     <- ctrl$K_train
  Tp    <- ctrl$T_probe
  del   <- ctrl$delta_shift
  tau   <- ctrl$tau_ppca_cut
  q     <- ctrl$q_chisq
  infl  <- ctrl$infl_chisq
  
  N     <- nrow(reduced_mat)
  p     <- ncol(reduced_mat)
  chi_th<- infl * stats::qchisq(q, df = Tp * p) / Tp
  
  cuts  <- N + 1L
  e     <- N
  while (e - K - Tp + 1 >= 1) {
    seg   <- reduced_mat[(e-K+1):e, , drop = FALSE]
    ppca  <- compute_PPCA_C(seg, tau)
    mu    <- ppca$mu
    Cinv  <- solve(ppca$C)
    
    probe <- reduced_mat[(e-K-Tp+1):(e-K), , drop = FALSE]
    diff  <- sweep(probe, 2, mu)
    Havg  <- mean(rowSums((diff %*% Cinv) * diff))
    
    if (Havg > chi_th) {
      cut_pt <- e - K
      cuts   <- c(cuts, cut_pt)
      e      <- cut_pt
    } else {
      e      <- e - del
    }
  }
  c(cuts, 1L)
}

divisive_split_idx <- function(x, tol, ids, depth = 0, tree_env = NULL) {
  if (is.null(tree_env)) tree_env <- list(lines = character())
  
  append_line <- function(txt) tree_env$lines <- c(tree_env$lines, txt)
  
  x   <- as.matrix(x)
  n   <- nrow(x)
  cen <- colMeans(x)
  
  if (n == 1) {
    append_line(sprintf("%sdepth=%d  n=1  spread=0.000\n%s→ leaf",
                        strrep("  ", depth), depth, strrep("  ", depth + 1)))
    return(list(groups = list(ids),
                tree_lines = tree_env$lines))
  }
  
  if (n == 2) {
    spread <- sqrt(sum((x[1, ] - cen)^2))
    
    append_line(sprintf("%sdepth=%d  n=2  spread=%.3f",
                        strrep("  ", depth), depth, spread))
    if (spread > tol) {
      append_line(sprintf("%s→ split into two singles",
                          strrep("  ", depth + 1)))
      groups <- list(ids[1], ids[2])
    } else {
      append_line(sprintf("%s→ 2‑row leaf (spread = %.3f)",
                          strrep("  ", depth + 1), spread))
      groups <- list(ids)
    }
    return(list(groups = groups,
                tree_lines = tree_env$lines))
  }
  
  
  ## n > 2 case
  dists  <- sqrt(rowSums((x - matrix(cen, n, ncol(x), byrow = TRUE))^2))
  spread <- max(dists)
  append_line(sprintf("%sdepth=%d  n=%d  spread=%.3f",
                      strrep("  ", depth), depth, n, spread))
  
  if (spread > tol) {
    km <- kmeans(x, centers = 2, nstart = 100)
    groups <- list()
    for (cl in 1:2) {
      sel   <- which(km$cluster == cl)
      res   <- divisive_split_idx(x[sel, , drop = FALSE],
                                  tol, ids[sel],
                                  depth + 1, tree_env)
      groups <- c(groups, res$groups)
    }
  } else {
    append_line(sprintf("%s→ leaf", strrep("  ", depth + 1)))
    groups <- list(ids)
  }
  list(groups = groups,
       tree_lines = tree_env$lines)
}


# reconstruction
reconstruct_swing <- function(df, principal_markers, forward_breaks, ctrl = list()) {
  ctrl <- merge_ctrl(ctrl)
  keep_cols <- unlist(lapply(principal_markers,
                             function(m) grep(paste0("^", m, "_"), colnames(df))))
  X   <- as.matrix(df[, keep_cols, drop = FALSE])
  Z   <- as.matrix(df[, -keep_cols, drop = FALSE])
  
  other_cols <- setdiff(seq_len(ncol(df)), keep_cols)
  
  n_seg <- length(forward_breaks) - 1L
  seg_models <- lapply(seq_len(n_seg), function(j) {
    idx <- forward_breaks[j]:(forward_breaks[j+1]-1)
    compute_PPCA_C(X[idx, ], tau = 0.75)   # fixed 0.75 for now
  })
  
  w     <- 0.8
  p     <- ncol(X)
  feature_mat <- t(vapply(seg_models, function(m) {
    v <- m$C[upper.tri(m$C, diag = TRUE)]
    c(w * v, (1 - w) * m$mu)
  }, numeric(p * (p + 3) / 2)))
  
  # pca + divisive clustering
  if (nrow(feature_mat) > 2) {
    pc <- prcomp(feature_mat, center = TRUE, scale. = TRUE)
    k95<- which(cumsum(pc$sdev^2) / sum(pc$sdev^2) >= 0.95)[1]
    red<- pc$x[, 1:k95, drop = FALSE]
  } else {
    red <- feature_mat
  }
  ids <- seq_len(nrow(red))
  tol_use <- ctrl$tol_split %||% ctrl$tol_divisive
  split_out <- divisive_split_idx(red, tol = tol_use, ids = ids)
  behaviour_label <- split_out$groups
  split_tree      <- split_out$tree_lines
  
  
  # frame‑level cluster ID
  seg_cluster_id <- integer(n_seg)
  for (cid in seq_along(behaviour_label))
    seg_cluster_id[behaviour_label[[cid]]] <- cid
  model_id <- integer(nrow(df))
  for (j in seq_len(n_seg))
    model_id[forward_breaks[j]:(forward_breaks[j+1]-1)] <- seg_cluster_id[j]
  
  # implement matrix math for least squares
  cluster_frames <- split(seq_len(nrow(df)), model_id)
  Bs <- mu_x <- mu_z <- vector("list", length(cluster_frames))
  for (cid in seq_along(cluster_frames)) {
    fr  <- cluster_frames[[cid]]
    Xc  <- t(X[fr, ])
    Zc  <- t(Z[fr, ])
    mu_x[[cid]] <- rowMeans(Xc)
    mu_z[[cid]] <- rowMeans(Zc)
    Xc  <- sweep(Xc, 1, mu_x[[cid]], "-")
    Zc  <- sweep(Zc, 1, mu_z[[cid]], "-")
    Bs[[cid]] <- Zc %*% t(Xc) %*% solve(Xc %*% t(Xc))
  }
  # reconstruct all frames
  rec <- matrix(NA_real_, nrow(df), ncol(df))
  for (f in seq_len(nrow(df))) {
    cid <- model_id[f]
    xf  <- X[f, ]
    zhat<- Bs[[cid]] %*% (xf - mu_x[[cid]]) + mu_z[[cid]]
    full<- numeric(ncol(df))
    full[keep_cols]  <- xf
    full[other_cols] <- zhat
    rec[f, ] <- full
  }
  colnames(rec) <- colnames(df)
  list(
    reconstructed   = rec,
    split_tree      = split_tree,
    feature_matrix  = feature_mat,
    behaviour_label = behaviour_label
  )
}


# process either a .csv file or df/matrix
process_swing <- function(input, ctrl = list()) {
  ctrl <- merge_ctrl(ctrl)
  if (is.character(input) && length(input) == 1L) {
    dat <- read.csv(input, check.names = FALSE)
  } else if (is.data.frame(input) || is.matrix(input)) {
    dat <- as.data.frame(input, check.names = FALSE)
  } else {
    stop("`input` must be a file path or a data.frame / matrix.")
  }
  
  # principal markers
  Ymat <- as.matrix(dat)
  pm   <- select_principal_markers(Ymat, var_thresh = ctrl$var_thresh_pfa,
                                   seed = ctrl$seed)
  
  xyz_cols <- unlist(lapply(pm, function(m) paste0(m, c("_X","_Y","_Z"))))
  red_mat  <- dat[, xyz_cols]
  
  # forward cuts
  cuts_fwd <- detect_ppca_cuts(red_mat, ctrl)
  
  # reconstruction & feature mat etc
  rec_out  <- reconstruct_swing(dat, pm, cuts_fwd)
  rec_mat  <- rec_out$reconstructed
  split_tree <- rec_out$split_tree
  
  
  list_out <- list(
    principal_markers = pm,
    breaks_forward    = cuts_fwd,
    behaviour_label   = rec_out$behaviour_label,
    feature_matrix    = rec_out$feature_matrix,    # heavy
    reconstructed     = rec_out$reconstructed      # heavy
  )
  class(list_out) <- "swing_obj"    # <- give it the custom print behaviour
  list_out
  
}

# now process an entire folder of .csv swings instead of just one
process_folder <- function(dir, pattern = "\\.csv$", ctrl = list()) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  res   <- lapply(files, process_swing, ctrl = ctrl)
  names(res) <- basename(files)
  res
}



process_swing_fixedpm <- function(input,
                                  pm,
                                  tmpl_hdr,
                                  tmpl_lc,
                                  ctrl = list()) {
  ctrl <- merge_ctrl(ctrl)
  
  # read file or DF
  dat <- if (is.character(input) && length(input) == 1L)
    read.csv(input, check.names = FALSE)
  else
    as.data.frame(input, check.names = FALSE)
  
  # force header names to template case
  dat <- apply_template_header(dat, tmpl_hdr, tmpl_lc)
  
  # subset to global PM coordinates
  xyz_cols <- unlist(lapply(pm, \(m) paste0(m, c("_X","_Y","_Z"))))
  red_mat  <- dat[ , xyz_cols, drop = FALSE]
  
  # unchanged
  cuts_fwd <- detect_ppca_cuts(red_mat, ctrl)
  rec_out  <- reconstruct_swing(dat, pm, cuts_fwd, ctrl)
  
  structure(list(
    principal_markers = pm,
    breaks_forward    = cuts_fwd,
    behaviour_label   = rec_out$behaviour_label,
    feature_matrix    = rec_out$feature_matrix,
    reconstructed     = rec_out$reconstructed
  ), class = "swing_obj")
}



# show what I want to see when I call something like big_out$hit080_10.csv in the batch run file

print.swing_obj <- function(x, ...) {
  cat("<swing_obj>\n")
  cat("$principal_markers :", paste(x$principal_markers, collapse = ", "), "\n")
  cat("$breaks_forward    :", paste(x$breaks_forward, collapse = ", "), "\n")
  cat("$behaviour_label   :", 
      paste(vapply(x$behaviour_label, toString, ""), collapse = " | "), "\n")
  invisible(x)
}



# get a swing summary
summarise_swings <- function(big_list) {
  data.frame(
    file          = names(big_list),
    n_PM          = vapply(big_list, function(x) length(x$principal_markers), 1L),
    n_breaks      = vapply(big_list, function(x) length(x$breaks_forward) - 1L, 1L),
    n_segments    = vapply(big_list, function(x) length(x$behaviour_label),    1L),
    stringsAsFactors = FALSE
  )
}


