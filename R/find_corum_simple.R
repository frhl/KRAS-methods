#' @title find corum complexes in data
#' @description finds all corum complexes in data
#' @param genes a vector of genes
#' @param species what species? 
#' @export

find_corum_simple<- function(genes, species = 'Human'){
  
  # find complexes
  corum <- all_corum_complexes
  corum <- corum[corum$organism %in% species,]
  corum <- corum[corum$genes %in% genes,]
  complexes <- unique(corum$complex)
  lst <- lapply(complexes, function(x) intersect(corum$genes[corum$complex %in% x],genes))
  names(lst) <- complexes
  lst <- null_omit(lst)
  return(lst)
  
}