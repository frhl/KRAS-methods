# what genesets are enriched in the WT versus MT.

library(ggplot2)
library(genoppi)

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
in_sheets <- in_sheets[grepl('starv', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data)
wt <- data$WT_Starvation
mt <- data[2:4]

# set background to prey of WT (starvation)
background <- unique(wt$Prey)

# MT minus WT
mt_minus_wt <- lapply(mt, function(x){

  x <- x[! x$Prey %in% background,]
  sig = data.frame(gene = unique(x$Prey[x$LogFC >= 4 & x$SaintScore > 0.8]), significant = T)
  insig = data.frame(gene = unique(x$Prey[! x$prey %in% sig$gene]), significant = F)
  d <- as.data.frame(rbind(sig, insig))
  return(d)
  
})

lapply(mt_minus_wt, function(x) table(x$significant))


# external database
go_table <- goa_cc_table
colnames(go_table) <- c('genes','id','geneset')
terms <- unique(go_table$geneset)
head(go_table)

result_conditional <- lapply(names(mt_minus_wt), function(sheet){
  
  print(sheet)

  go_enrichment <- do.call(rbind, lapply(terms, function(term){
    
    #print(term)
    # get apex2 dataset
    df <- mt_minus_wt[[sheet]]
    
    # get go subset
    go_df <- go_table[go_table$geneset %in% term,]
    go_df$significant <- TRUE
    go_df <- go_df[,c('genes','significant')]
    
    # calculate overlap
    hypergeom <- calc_hyper(df, go_df, bait = 'KRAS', intersectDf = data.frame(intersectN = F))
    outdf <- hypergeom$statistics
    outdf$list_name <- sheet
    outdf$dataset <- term
    outdf$comparison <- 'SaintScore >= 0.8, LogFC > 4 (Conditional Enrichment)'
    outdf$overlap_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
    return(outdf)
  }))
  
  go_enrichment$FDR <- stats::p.adjust(go_enrichment$pvalue, method = 'fdr')
  go_enrichment <- go_enrichment[,c(1,8,9,2:7,11,10)]
  go_enrichment <- go_enrichment[order(go_enrichment$pvalue), ]
  go_enrichment <- go_enrichment[nchar(go_enrichment$overlap_genes) > 1,]

  return(go_enrichment)
  
})

names(result_conditional) <- names(mt_minus_wt)
write_xlsx(result_conditional,'derived/tables/210406_goa_cc_table_conditional_enrichment_mt_starv_minus_wt_logfc4.xlsx')

# heatmap of enrichment pathways
hm_data <- do.call(rbind, lapply(result_conditional, function(x) data.frame(go = x$dataset, pvalue = x$pvalue, experiment = x$list_name)))
exclude <- unlist(lapply(unique(hm_data$go), function(x) all(hm_data$pvalue[hm_data$go %in% x] > 0.5)))
hm_data <- hm_data[hm_data$go %in% unique(hm_data$go)[!exclude],]

ggplot(hm_data, aes(y = reorder(go, -log10(pvalue)), x = experiment, fill = -log10(pvalue))) +
  geom_tile() + theme_minimal() + 
  scale_fill_gradient(low = 'white', high = 'red') +
  xlab('Experiment - WT') + 
  ylab('Gene Ontology (CC)') +
  ggtitle('GO Enrichment Analysis', 
          'Experiment (Significant if LogFC >= 4 & SaintScore > 0.8)')

ggsave('derived/plots/210406_goa_cc_table_conditional_enrichment_mt_starv_minus_wt_logfc4.pdf', width = 8, height = 15)
