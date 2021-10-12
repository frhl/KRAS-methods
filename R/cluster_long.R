cluster_long <- function(long, formula = as.formula(gene_lysine ~ comparison)){
  
  long_copy <- long
  long$gene <- NULL
  wide <- dcast(long, formula)
  wide_na_replaced <- wide
  wide_na_replaced[is.na(wide_na_replaced)] <- 0
  wide_mat <- wide_na_replaced
  wide_mat$gene_lysine <- NULL
  clust <- hclust(dist(wide_mat))
  levels <- as.character(wide_na_replaced$gene_lysine)[clust$order]
  long_copy$gene_lysine <- factor(long_copy$gene_lysine, levels = levels)
  return(long_copy)
  
}

