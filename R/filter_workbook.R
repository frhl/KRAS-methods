#' @title filter workbook
#' @description filters a workbook given certain criteria using dplyr::filter.
#' @param workbook a list of data.frames
#' @param ... see ?dplyr
#' @export

filter_workbook <- function(workbook, ...){
  result <- lapply(workbook, function(wb){
    dplyr::filter(wb, ...)
  })
  return(result)
}




