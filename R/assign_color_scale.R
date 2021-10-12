#' @title assing color scale
#' @description creates a color scale for an input data.frame based on their
#' max/min value (unless force_min / force_max) are specified. This color scale
#' will be appended as a column to the graph.
#' @param d a data.frame with the column 'val'
#' @param colors vector. The three colors used for the color scale. 
#' @param force_min numeric. Force the minimum value of the color scale.
#' @param force_max numeric. Force the maximum value of the color scale.
#' @param force_length_out integer. Force unique colors in scale.
#' @param color_not_in_scale string. Color for items not on the scale.
#' 
#' @family graphs
#' @export

assign_color_scale <- function(d, colors = c("yellow",'white',"red"), force_min = NULL, force_max = NULL, force_length_out = NULL, color_not_in_scale = 'blue'){
  
  require(RColorBrewer)
  # check input
  stopifnot(c('val') %in% colnames(d))
  col_func <- colorRampPalette(colors)
  
  # setup intervals (default is using information in the data)
  min_val = as.numeric(ifelse(is.null(force_min), floor(min(d$val)), force_min))
  max_val = as.numeric(ifelse(is.null(force_min), ceiling(max(d$val)), force_max))
  length_out = as.numeric(ifelse(is.null(force_min), length(unique(d$val)), force_length_out))
  
  # get color intervals and chop up input data  
  color_interval <- seq(min_val, max_val, length.out = length_out)
  color_labels <-  col_func(length(color_interval))[-1]
  d$color <- as.character(cut(d$val, breaks = color_interval, labels = color_labels))
  
  if (any(is.na(d$color))) warnings(paste('some indexes were not assigned to a color. These are colored as ',color_not_in_scale ,'. Increase force_min/force_max to resolve this.' ))
  d$color[is.na(d$color)] <- color_not_in_scale # blue indicates not assigned
  
  
  return(d)
  
}
