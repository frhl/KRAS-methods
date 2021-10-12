# libraries (all are strictly required to run the following code)
devtools::load_all()
library(igraph)
library(qgraph)
library(RColorBrewer)
library(genoppi) # can be downloaded from github

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- lapply(data, function(x){
  x$significant <- x$SaintScore >= 0.8 & x$LogFC > 4
  return(x)
})


# external database
go_table <- goa_bp_table
colnames(go_table) <- c('genes','id','geneset')
terms <- unique(go_table$geneset)
head(go_table)

result_conditional <- lapply(names(data), function(sheet){
  
  print(sheet)
  count = 0
  
  go_enrichment <- do.call(rbind, lapply(terms, function(term){
    
    # print to screen
    count <<- count + 1
    #print(term)
    #if (count %% 1000) print(paste0(term,': ', count, '/', length(terms)))
    
    # get apex2 dataset
    df <- as.data.frame(data[[sheet]][,c('Prey','significant')])
    colnames(df)[1] <- 'genes'
    
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
  go_enrichment <- go_enrichment[order(go_enrichment$FDR), ]
  return(go_enrichment)
  
})
names(result_conditional) <- names(data)
#write_xlsx(result_conditional,'derived/tables/210319_goa_bp_table_conditional_enrichment_logfc4.xlsx')
