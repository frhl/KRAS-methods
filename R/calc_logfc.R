#' @title Transform data
#' @description Calculates log fold change by substracting
#' a case column with a bait column. Assumes the data is
#' already in log2 intensities.
#' @param df a data.frame with bait and controls
#' @author flassen
#' @note assummes log2(bait1), log2(control1), log2(bait2).. 
#' @family processing
#' @export

calc_logfc <- function(df){
  
  dfNum <- sapply(df, is.numeric)
  dfNum <- df[, dfNum]
  pairs <- (ncol(dfNum))/2
  stopifnot(pairs==round(pairs))
  
  for (i in 1:pairs){
    colMock = dfNum[, (i*2)]
    colBait = dfNum[, (i*2)-1]  # log2(bait) - log2(control)
    dfNum[[paste0('rep',i)]] = colBait - colMock
  }
  
  return(cbind(df, dfNum[,grepl('rep', colnames(dfNum))]))
}