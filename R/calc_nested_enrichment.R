#' @title calculate nested enrichment
#' @description use hypergeometric overlap to calculate the one-tailed 
#' enrichment of pathways.
#' @param data a list of data.frames. Expects this data.frame to contain
#' the columns 'gene' and 'significant' 
#' @param pathways a data.frame with at least columns gene and geneset
#' @param bait the bait
#' @param intersectN should thw two total populations be merged?
#' @param saintscore the saintscore threshold (greater-than or equal).
#' @param logfc the logfc threshold (greater-than or equal)
#' @export

calc_nested_enrichment <- function(data, pathways, bait = 'KRAS', intersectN = F){
  
  if (!all(c('gene','geneset') %in% colnames(pathways))) stop('require gene/geneset column!')
  if (!all(c('gene', 'significant') %in% colnames(data)))
  
  # Find out what complexes are enriched in the data
  enrichment <- lapply(names(data), function(cur_sheet){
    cur_data <- data[[cur_sheet]]
    #cur_data <- data.frame(gene = cur_data$Prey, significant = cur_data$SaintScore >= saintscore & cur_data$LogFC >= logfc)
    #print(cur_sheet)
    pathways <- unique(signature_pathways$geneset)
    my_enrichment <- lapply(pathways , function(pathway_name){
      cur_pathway <- signature_pathways
      cur_pathway$significant <- cur_pathway$geneset %in% pathway_name
      cur_pathway <- data.frame(gene = cur_pathway$gene, significant = cur_pathway$significant)
      hyper <- calc_hyper(cur_data, cur_pathway, intersectDf = data.frame(intersectN = intersectN), bait = bait)
      hyper$statistics$list_name <- pathway_name
      return(hyper$statistics)
    })
    enrichment_sorted <- do.call(rbind, my_enrichment)
    enrichment_sorted <- enrichment_sorted[(order(enrichment_sorted$pvalue)),]
    enrichment_sorted$FDR <- stats::p.adjust(enrichment_sorted$pvalue, method = 'fdr')
    return(enrichment_sorted)
  })
  
  return(enrichment)
  
}