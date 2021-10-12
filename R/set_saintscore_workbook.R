#' @title set a notbook as a saintscore notbook
#' @description dafafa
#' @param workbook a list of data.frames
#' @export


set_saintscore_workbook <- function(workbook){
  
  result <- lapply(workbook, function(wb){

    # read in the data
    stopifnot('FoldChange' %in% colnames(wb))
    stopifnot('Bait' %in% colnames(wb))
    stopifnot('Prey' %in% colnames(wb))
    stopifnot('BFDR' %in% colnames(wb))
    stopifnot('SaintScore' %in% colnames(wb))
    
    # get LogFC column
    wb$LogFC <- log2(wb$FoldChange)
    
    # extract releveant columns
    wb <- wb[,c('Bait','Prey', 'LogFC', 'SaintScore', 'BFDR', 'Sheet')]
    return(wb)
    
  })
  
  return(result)
  
}




