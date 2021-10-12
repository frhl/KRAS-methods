#' @title Normalize data
#' @description Normalize a table of label or label free data. 
#' The function will calculate the \code{type} e.g. median of
#' all numeric columns and subtract the median from the 
#' intensity value.
#' @author flassen
#' @export

normalize_intensities <- function(table = NULL, type = 'median'){
  
  columns_numeric <- sapply(table, is.numeric)
  table_numeric <- table[, columns_numeric]
  values_type <- as.data.frame(t(sapply(table_numeric, function(x) do.call(type, list(na.omit(x))))))
  values_type <- values_type[rep(seq_len(nrow(values_type)), each = nrow(table)), ]
  table_normalized <- table_numeric - values_type
  table[, columns_numeric] <- table_normalized
  return(table)

}




