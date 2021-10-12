#' @title setup graph legend
#' @description will setup a graph legend according to a decided 
#' threshold (min_color to max_color) values. Also uses a colorfuncn to generate
#' color scheme. This can be generated with \code{assign_color_scale}.
#' @param min_color value. Minimum value for lowest color.
#' @param max_color value. Maximum value for highest color.
#' @param colorfunc function. A function that returns a colorscale in a data.frame
#' @param draw_every numeric. How often should text be drwan?
#' @param text character. Title of the legend.
#' 
#' @note values not in colco scale witll me removed in warning. This can
#' be fixed by changing min_color/max_color.
#' 
#' @export

setup_graph_legend <- function(min_color, max_color, colorfunc, draw_every = 0.5, text = 'Log(WT/MT)'){
  
  # integrate data
  ticks <- data.frame(val = seq(floor(min_color), ceiling(max_color), by = 0.01))
  legend <- data.frame(
    index = 1,
    ticks = ticks$val,
    color = colorfunc(ticks)
  )
  
  # remove blue (unassigned colors)
  errors <- legend[legend$color.color == 'blue',]
  if (nrow(errors) > 0){
    warning('edges were removed since min_color/max_color is not sufficient to cover edge cases.')
    legend <- legend[!legend$color.color %in% 'blue',] 
  }
  
  # generate legend by stacking square points
  plt <- ggplot(legend, aes(x=index, y = ticks, label = round(ticks,3))) +
    geom_point(size = 15, shape=22, color = legend$color.color, fill = legend$color.color) + 
    geom_text(data = legend[legend$ticks %% draw_every == 0,], size = 3) +
    annotate(x = 1, y = max(legend$ticks)*1.15, label = text, geom = 'text', fontface = 'bold') +
    theme_void()
  
  
  return(plt)
  
}