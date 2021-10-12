#' @title collapse workbook
#' @description collapses a list of data.frames.
#' @param workbook a list of data.frames
#' @family excel
#' @export

collapse_workbook <- function(workbook){

  df <- as.data.frame(do.call(rbind, workbook))
  rownames(df) <- NULL
  return(df)

}