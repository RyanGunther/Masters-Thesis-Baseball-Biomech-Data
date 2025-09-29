library(plotly)

list_markers <- function(data) sub("_X$", "", grep("_X$", names(data), value = TRUE))

# going to plot in interactive 3d so I can rotate and the view the image from all angles
# should also be able to connect points i.e. recreate simple limb and bat structure

plot_frames_3d <- function(data, frame_idxs,
                           label_markers = NULL,
                           connections = NULL,
                           dot_col = "grey30",
                           label_col = "red",
                           line_col = "blue",
                           line_width = 3,
                           marker_size = 3) {
  # find marker base names i.e. without _x/y/z suffix
  x_cols <- grep("_X$", names(data), value = TRUE)
  y_cols <- sub("_X$", "_Y", x_cols)
  z_cols <- sub("_X$", "_Z", x_cols)
  base_names <- sub("_X$", "", x_cols)
  
  # label each index
  label_idx <- integer(0)
  if (!is.null(label_markers)) {
    idx <- match(toupper(label_markers), toupper(base_names))
    keep <- !is.na(idx)
    label_idx     <- idx[keep]
    label_markers <- label_markers[keep]
  }
  
  plots <- vector("list", length(frame_idxs))
  for (k in seq_along(frame_idxs)) {
    fi <- frame_idxs[k]
    x_vals <- as.numeric(data[fi, x_cols, drop = TRUE])
    y_vals <- as.numeric(data[fi, y_cols, drop = TRUE])
    z_vals <- as.numeric(data[fi, z_cols, drop = TRUE])
    
    # might not use the hover names but going to include just in case
    p <- plot_ly(
      x = x_vals, y = y_vals, z = z_vals,
      type = "scatter3d", mode = "markers",
      marker = list(color = dot_col, size = marker_size),
      text = base_names, hoverinfo = "text"
    ) %>%
      layout(
        scene = list(
          xaxis = list(title = "X (m)"),
          yaxis = list(title = "Y (m)"),
          zaxis = list(title = "Z (m)")
        ),
        title = paste("Frame", fi)
      )
    
    if (length(label_idx)) {
      p <- add_trace(
        p,
        x = x_vals[label_idx],
        y = y_vals[label_idx],
        z = z_vals[label_idx],
        type = "scatter3d", mode = "markers+text",
        text = label_markers, textposition = "top right",
        marker = list(color = label_col, size = marker_size + 2),
        showlegend = FALSE
      )
    }
    
    # include connections between points
    if (!is.null(connections)) {
      for (conn in connections) {
        if (length(conn) != 2) next
        idx <- match(toupper(conn), toupper(base_names))
        if (any(is.na(idx))) next
        p <- add_trace(
          p,
          x = x_vals[idx], y = y_vals[idx], z = z_vals[idx],
          type = "scatter3d", mode = "lines",
          line = list(color = line_col, width = line_width),
          text = NULL, hoverinfo = "none", showlegend = FALSE
        )
      }
    }
    
    plots[[k]] <- p
  }
  
  if (length(plots) == 1) plots[[1]] else plots
}


### usage
currentplot <- read.csv("000099_000219_72_187_L_005_947.csv")
list_markers(currentplot) # for connections sake

plots <- plot_frames_3d(
  currentplot,
  frame_idxs = 432:449,
  label_markers = c("L_Wrist", "R_Wrist"),
  connections = list(
    # these are optional but if viewing from github, take a look through the list_markers and pick what you want
    # for a static visualization of what each marker name means, see here: https://github.com/drivelineresearch/openbiomechanics/tree/main/baseball_hitting
    c("Marker1", "Marker2"),
    c("Marker1", "Marker3"),
    c("RSHO", "RELB"),
    c("LSHO", "LELB"),
    c("RELB", "RFIN"),
    c("LELB", "LFIN"),
    c("RASI", "RKNE"),
    c("LASI","LKNE"),
    c("RKNE","RANK"),
    c("LKNE","LANK"),
    c("RANK","RTOE"),
    c("LANK","LTOE")
  )
)

# in this case plot 1 corresponds to the first plot of the range indicated in line 96
# hence plot 1 is really frame 432
plots[[1]]
plots[[5]]
plots[[10]]
plots[[14]]
plots[[18]]

# can be viewed in Rstudio in the viewer section. Also possible to pop out to a browser tab where you have more space and ability to bounce between plots