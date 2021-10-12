get_prey_complexes <- function(preys, db){
  
  stopifnot(is.vector(preys))
  stopifnot(length(preys) > 0)
  stopifnot('genes' %in% colnames(db))
  stopifnot('complex' %in% colnames(db))
  
  prey_df <- data.frame(genes = preys)
  merge_df <- merge(prey_df, db)
  count_complexes <- sort(table(merge_df$complex))
  
  table_complexes <- lapply(rev(names(count_complexes)), function(current_complex){
    
    complexes_total_genes <- db$genes[db$complex %in% current_complex]
    complexes_sample_genes <- merge_df$genes[merge_df$complex %in% current_complex]
    complexes_total_genes_n <- length(complexes_total_genes)
    complexes_sample_genes_n <- length(complexes_sample_genes)
    complexes_pct_found <- round(complexes_sample_genes_n / complexes_total_genes_n,3)*100
    d <- data.frame(complex = current_complex, n_found =  complexes_sample_genes_n, n_total = complexes_total_genes_n, pct = complexes_pct_found,
                    genes_found = paste(complexes_sample_genes, collapse = ';'), genes_total = paste(complexes_total_genes, collapse = ';'))
    return(d)
  })
  
  table_complexes <- do.call(rbind, table_complexes)
  return(table_complexes)
}

