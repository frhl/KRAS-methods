#' @title get the right ms data columns
#' @description a function used to handle Oxford Mass Spectromoetry intensity files.
#' @param d a data.frame (typically an output from Mascott). This should ideally
#' contain colnames like this "Intensity 37", Intensity 38", "Intensity 39". Where
#' the first string indicates the name of column and the number indicates the well.
#' @param what name of column that is used for regex
#' @param master a master sheet that is used for subsetting. Should contain
#' at least one column 'well' that indicates the tube/sample of interest.
#' @param id what is the id column? Wil be keept alongside regexed columns
#' 
#' @family processing
#' @export


grep_ms_data <- function(d, what = 'Intensity', master, id = "Gene names"){
  
  # read master sheets which contains experiment well IDs
  #master <- read.csv('inst/extdata/kras_wells_master.csv')
  #master <- master[!is.na(master$experiment),]
  #master$experiment <- gsub(' ', '_', master$experiment)
  
  # read MS data
  #d <- fread('../data/proteinGroups Kras.txt')
  #str(d)
  #ncol(d)
  stopifnot(id %in% colnames(d))
  
  # preliminary remove cols
  stopifnot(id %in% colnames(d))
  d_subset <- d[, grepl(what, colnames(d)), with = F]
  d_subset <- d_subset[,grepl('[0-9]+', colnames(d_subset)), with = F]
  
  # only keep columns also in master sheet
  names <- colnames(d_subset)
  matches <- as.integer(unlist(regmatches(names, gregexpr("[[:digit:]]+", names))))
  newnames <- merge(data.frame(names, matches), master, by.x = 'matches', by.y = 'well')
  d_new <- d_subset[,c(newnames$names), with = F]
  colnames(d_new) <- newnames$experiment
  d_new[[id]] = d[[id]]
  
  return(d_new)
}
