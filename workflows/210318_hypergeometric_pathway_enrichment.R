# libraries (all are strictly required to run the following code)
devtools::load_all()
library(igraph)
library(qgraph)
library(RColorBrewer)
library(genoppi) # can be downloaded from github. Ask frederik.

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- lapply(data, function(x){
  x$significant <- x$SaintScore >= 0.8 & x$LogFC > 0
  return(x)
})

# genesets
lawrence2014 <- fread('inst/extdata/lawrence2014cancergenes.csv')
cancer_genes <- data.frame(genes = lawrence2014$geneName, significant = TRUE)
neuro <- fread('inst/extdata/Ripke2014_Grove2018_Satterstrom2018.csv', header = T)
neuro_genes <- data.frame(genes = neuro$hgnc_symbol, significant = TRUE)



#############################
# conditional on background #
#############################

# calculate hypergeometric overlap with lawrence cancer genes
result_conditional <- do.call(rbind, lapply(names(data), function(sheet){
  
  df <- as.data.frame(data[[sheet]][,c('Prey','significant')])
  colnames(df)[1] <- 'genes'
  hypergeom <- calc_hyper(df, cancer_genes, bait = 'KRAS', intersectDf = data.frame(intersectN = F))
  outdf <- hypergeom$statistics
  outdf$cancer_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
  outdf$list_name <- sheet
  outdf$comparison <- 'SaintScore >= 0.8, LogFC > 0 (Conditional Enrichment)'
  return(outdf)
  
}))

# move columns
result_conditional$FDR <- stats::p.adjust(result_conditional$pvalue, method = 'fdr')
result_conditional <- result_conditional[,c(1,9,2:7,10,8)]
#write.table(result_conditional, 'derived/tables/210319_conditional_hypergeoemtric_enrichment_analysis_logfc4.txt', row.names = F, quote = F, sep = '\t')

##############################
# global enrichment analysis #
##############################

# get all genes (let's just use InWeb)
all_genes <- unique(c(inweb_table$Gene1, inweb_table$Gene2))
df_all <- data.frame(genes = all_genes, significant = F)

# calculate hypergeometric overlap with lawrence cancer genes
result_global <- do.call(rbind, lapply(names(data), function(sheet){
  
  # use all known proteins as background instead
  df <- as.data.frame(data[[sheet]][,c('Prey','significant')])
  colnames(df)[1] <- 'genes'
  df <- df[df$significant,]
  df <- rbind(df, df_all)
  
  # calculate hypergeometric overlap test
  hypergeom <- calc_hyper(df, cancer_genes, bait = 'KRAS', intersectDf = data.frame(intersectN = F))
  outdf <- hypergeom$statistics
  outdf$cancer_genes <- paste(hypergeom$genes$mylist$successInSample_genes, collapse = ';')
  outdf$list_name <- sheet
  outdf$comparison <- 'SaintScore >= 0.8, LogFC > 4 (Global Enrichment)'
  return(outdf)
  
}))

# move columns
result_global$FDR <- stats::p.adjust(result_global$pvalue, method = 'fdr')
result_global <- result_global[,c(1,9,2:7,10,8)]
#write.table(result_global, 'derived/tables/210319_global_hypergeoemtric_enrichment_analysis_logfc4.txt', row.names = F, quote = F, sep = '\t')

# combine data into one table
result_conditional_minimal <- result_conditional[,c('list_name', 'pvalue')]
result_conditional_minimal$analysis <- 'Conditional Enrichment'
result_global_minimal <- result_global[,c('list_name', 'pvalue')]
result_global_minimal$analysis <- 'Global Enrichment'
result_combined <- rbind(result_conditional_minimal, result_global_minimal)
result_combined$logpvalue <- -log10(result_combined$pvalue)

# set significance thresholds
bonf <- 0.05/8
ggbarplot(result_combined, bonf) + 
  ggtitle('Enrichment of Cancer Genes in KRAS Experiments',
          'Thresholded by SaintScore > 0.8 and LogFC > 0 (Hypergeometric overlap)') +
  facet_wrap(~analysis)


ggsave('derived/plots/210323_hypergeoemtric_enrichment_cancer_analysis_logfc0.pdf', width = 8, height = 6)

