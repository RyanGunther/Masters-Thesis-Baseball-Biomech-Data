# quick 2d plot that I use to sanity check what's happening at that point in the swing:
# (will also upload my 3d plotting script that I use for more involved visualizations)
# R studio will allow me to plot up to 100 frames of the swing at one time
plot_frames_single <- function(data, frame_idxs, label_markers = NULL,
                               dot_col = "grey30", label_col = "red",
                               xlim_fixed = c(-1.25, 1), ylim_fixed = c(0, 2)) {
  ## find X / Z columns (can use Y alternatively if needed)
  x_cols <- grep("_X$", names(data), value = TRUE)
  z_cols <- sub("_X$", "_Z", x_cols)
  
  ## indices of markers to label
  label_idx <- integer(0)
  if (!is.null(label_markers)) {
    idx <- match(toupper(label_markers), toupper(sub("_X$", "", x_cols)))
    keep <- !is.na(idx)
    label_idx      <- idx[keep]
    label_markers  <- label_markers[keep]
    
  }
  
  op <- par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
  on.exit(par(op), add = TRUE)
  
  for (fi in frame_idxs) {
    x_vals <- as.numeric(data[fi, x_cols])
    z_vals <- as.numeric(data[fi, z_cols])
    
    plot(
      x_vals, z_vals,
      pch  = 16, col = dot_col,
      xlab = "X location (m)",
      ylab = "Z (height, m)",
      main = paste("Frame", fi),
      xlim = xlim_fixed,
      ylim = ylim_fixed
    )
    
    ## highlight & label requested markers
    if (length(label_idx)) {
      points(x_vals[label_idx], z_vals[label_idx], pch = 19, col = label_col)
      text(
        x_vals[label_idx], z_vals[label_idx],
        labels = label_markers,
        pos = 4, offset = 0.25, cex = 0.8, col = label_col
      )
    }
  }
}


currentplot <- read.csv("000080_000282_71_188_R_001_1012.csv")
plot_frames_single(currentplot, frame_idxs = 650:720)



