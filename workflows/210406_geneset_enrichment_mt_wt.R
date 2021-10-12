# what genesets are enriched in the WT versus MT.
setwd('~/Projects/06_kessler_rotation/MassSpectrometry/')
library(readxl)
library(ggplot2)
library(genoppi)
devtools::load_all()

# data path
path <- 'inst/extdata/KRAS/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
in_sheets <- in_sheets[grepl('(starv)|(fcs)', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data)
data <- lapply(data, function(x) {
  x$significant <- x$SaintScore >= 0.8 & x$LogFC >= 4
  return(x)
  })



# 
datasets <- list(
  msigdb_h = msigdb_h_table
  #go_bp = goa_bp_table[,c(1,3)],
  #go_cc = goa_cc_table[,c(1,3)],
  #go_mf = goa_mf_table[,c(1,3)],
  #msigdb_c5 = msigdb_c5_table,
  #corum = all_corum_complexes[all_corum_complexes$organism == 'Human',][,c(1,2)],
  #reactome = all_reactome_complexes[,c(1,4)]
)

# run all datasets
for (dataset_name in names(datasets)){
  
  print(dataset_name)
  dataset = datasets[[dataset_name]]
  
  # external database
  go_table <- dataset
  #go_table$GO.ID <- NULL
  colnames(go_table) <- c('genes','geneset')
  #colnames(go_table) <- c('genes','id','geneset')
  terms <- unique(go_table$geneset)
  head(go_table)
  
  result_conditional <- lapply(sort(names(data)), function(sheet){
    
    print(sheet)
    
    go_enrichment <- do.call(rbind, lapply(terms, function(term){
      
      # get apex2 dataset
      cur <- data[[sheet]]
      df <- data.frame(gene = cur$Prey, significant = cur$significant)
      
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
    go_enrichment$how <- ifelse(grepl("Starv", go_enrichment$list_name), 'Starvation', 'FCS')
    go_enrichment <- go_enrichment[nchar(go_enrichment$overlap_genes) > 1,]
    
    # return 
    if (nrow(go_enrichment) == 0){
      return(NULL)
    } else {
      return(go_enrichment)
    }

    
  })
  
  # rename data
  names(result_conditional) <- sort(names(data))
  order <- c('WT_Starvation','G12_Starvation','G13_Starvation','Q61H_Starvation',
             'WT_FCS', 'G12_FCS',  'G13_FCS_Averange', 'Q61H_FCS')
  result_conditional <- result_conditional[order]
  
  # get an overview
  overview <- as.data.frame(t(data.frame(
    dataset = dataset_name,
    dataset_size = length(unique(dataset[,2])),
    analysis = 'Proteome conditional enrichment analysis (One-tailed Hypergeometric test)',
    criteria = 'SaintScore >= 0.8 and LogFC >= 4'
   )))
  overview$V2 <- rownames(overview)
  rownames(overview) <- NULL
  overview <- overview[,c(2,1)]
  overview_sheets <- as.data.frame(cbind('Enrichment Analysis',(matrix(names(result_conditional)))))
  overview_final <- rbind(overview, overview_sheets)
  colnames(overview_final) <- c('Info','Value')
  
  result_conditional[['Description']] <- overview_final
  result_conditional <- result_conditional[c('Description',order)]
  
  # omit NULL
  result_conditional <- null_omit(result_conditional)
  
  # write analysis
  outfile = paste0('derived/KRAS/tables/210901_',dataset_name,'_enrichment_logfc4_s8.xlsx')
  write_xlsx(x = result_conditional, path = outfile)
  
}

#################
### plot data ###
#################

paths <- list.files('derived/tables/', pattern = '210413', full.names = T)
paths <- paths[c(2,3,4,6)]

for (path in paths){
#path <- 'derived/tables/210413_msigdb_h_enrichment_logfc4_saintscore8.xlsx'
#path <- 'derived/'
infile_sheets <- readxl::excel_sheets(path)
result_conditional <- lapply(infile_sheets, function(sheet) read_excel(path, sheet))
names(result_conditional) <- infile_sheets

# heatmap of enrichment pathways
hm_data <- do.call(rbind, lapply(result_conditional, function(x) data.frame(go = x$dataset, pvalue = x$pvalue, FDR = x$FDR, experiment = x$list_name, how = x$how)))
exclude <- unlist(lapply(unique(hm_data$go), function(x) all(hm_data$FDR[hm_data$go %in% x] > 0.2)))
hm_data <- hm_data[hm_data$go %in% unique(hm_data$go)[!exclude],]
hm_data$experiment <- gsub('_Averange', '',hm_data$experiment)

# setup labels
hm_data$label <- ifelse(hm_data$pvalue < 0.05, ifelse(hm_data$FDR < 0.1, '**','*'),'')

# setup rows
hm_data$row <- as.factor(ifelse(grepl('Starv',hm_data$experiment),'Starvation','FCS'))
hm_data$experiment <- unlist(lapply(strsplit(hm_data$experiment, split = '_'), function(x) x[1]))
hm_data$experiment <- factor(hm_data$experiment, levels = c("WT","G12", "G13","Q61H"))
hm_data$how <- factor(hm_data$how, levels = c("Starvation","FCS"))


p <- ggplot(hm_data, aes(y = reorder(go, -log10(FDR)), x = experiment, label = label, fill = -log10(FDR))) +
    geom_tile() + theme_bw() + geom_text() +
    scale_fill_gradient(low = 'white', high = 'red') +
    #scale_fill_gradient(low = 'white', high = unlist(ifelse(x == 'FCS', 'red','blue'))) +
    xlab('Experiment - WT') + 
    ylab('MSigDB Hallmark Genesets') +
    ggtitle(basename(path)) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    facet_wrap(~how)
    #theme(legend.position= unlist(ifelse(x == 'FCS', 'left','right'))) +
    #if (x == "FCS") scale_y_discrete(position = "left") else scale_y_discrete(position = "right")

outname <- gsub('tables','plots',gsub('\\.xlsx','.pdf',path))
ggsave(outname, p, width = 9, height = 10)

}
#ggave(p, file'derived/plots/210412_kras_go_bp_conditional_heatmap.pdf', width = 14, height = 14)

#pdf(, width = 14, height = 14)

#graphics.off()
  



#ggsave('derived/plots/210406_goa_mf_table_conditional_enrichment_mt_wt_logfc4.pdf', width = 8, height = 15)
