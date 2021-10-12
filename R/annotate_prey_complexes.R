annotate_prey_complexes <- function(preys, db, n = 2){
  
  stopifnot(is.vector(preys))
  stopifnot(length(preys) > 0)
  stopifnot('genes' %in% colnames(db))
  stopifnot('complex' %in% colnames(db))
  
  prey_df <- data.frame(genes = preys)
  merge_df <- merge(prey_df, db)
  count_complexes <- table(merge_df$complex)
  complexes <- names(count_complexes)[count_complexes >= n]
  
  data_complexes <- lapply(complexes, function(complex){
    prey_group <- merge_df[merge_df$complex %in% complex,]
    prey_prey_df <- as_interactor_complex(prey_group$genes)
    prey_prey_df <- prey_prey_df[!duplicated_interaction(prey_prey_df),]
    prey_prey_df <- prey_prey_df[prey_prey_df$int1 !=  prey_prey_df$int2,]
    return(prey_prey_df)
  })
  
  names(data_complexes) <- complexes
  return(data_complexes)
}




