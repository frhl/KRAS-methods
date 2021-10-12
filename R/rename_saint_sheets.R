# a helper function for renameing lists
rename_saint_sheets <- function(x){
  x <- gsub('SAINTexpress\\ +', '',x)
  x <- gsub('^\\ ', '', x)
  x <- gsub('\\ $', '', x)
  x <- gsub(' ','_',x)
  return(x)
}
