#' @title read excel sheets
#' @description reads in excel sheets alongside their sheet-names as a list of tibbles.
#' @param path the path to the excel file.
#' @return a list of tibbles
#' @export
#' @family excel

read_excel_sheets <- function(path){
  if (!file.exists(path)) stop(paste(file,'does not exist!'))
  sheets <- readxl::excel_sheets(path)
  data <- lapply(sheets, function(sheet) readxl::read_excel(path, sheet))
  names(data) <- sheets
  if (length(unique(sheets)) != length(sheets)) warning('Non-unique excel sheets names! Future indexing may be erroneous.')
  return(data)
}