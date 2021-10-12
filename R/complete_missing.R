#' @title Complete missing
#' @description Take a long data.frame and adds in missing values that would otherwise not be shown, when 
#' transforming the data to wide. This can primarily be used for plotting heatmaps.
#' @param df data.frame
#' @param id_var reshape2 id_var
#' @param time_var reshape2 time_var
#' @param value_var String. value column
#' @param complete what should the missing values be replaced with?
#' @param return_as string. Either "long" or "wide"
#' 
#' @family plotting
#' @export


complete_missing <- function(df, id_var = "Prey", time_var = "Sheet", value_var = 'LogFC', complete = 0, return_as = 'long'){
  
  # first convert to wide
  df <- df[,c(id_var, time_var, value_var)]
  df_wide <- reshape(df, idvar = id_var, timevar = time_var, direction = 'wide')
  df_wide[is.na(df_wide)] <- complete
  if (return_as == 'wide') return(df_wide)
  
  # convert back to long
  df_long <- melt(df_wide)
  colnames(df_long) <- c(id_var, time_var, value_var)
  df_long[[time_var]] <- gsub(paste0(value_var,'.'),'',df_long[[time_var]])
  if (return_as == 'long') return(df_long)
  
}


