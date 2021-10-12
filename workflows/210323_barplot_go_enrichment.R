devtools::load_all()

# plot the results
files <- list.files(path = 'derived/tables/', pattern = '210319', full.names = T)
files <- files[grepl('\\.xlsx',files)]

graphics.off()
for (f in files){
  
  table = paste0(unlist(strsplit(f, split = '_'))[2:3],collapse = '_')
  outfile = paste('derived/plots/210319_',table,'_enrichment_barplots.pdf')
  pdf(outfile, width = 9, height = 6)
  sheets <- excel_sheets(f)
  
  for (sheet in sheets){
    
    # process file
    infile <- as.data.frame(read_xlsx(f, sheet))
    bonf = 0.05/ length(infile$dataset)
    infile <- head(infile, 25)
    infile$experiment <- infile$list_name
    infile$list_name <- infile$dataset
    
    title = paste0(sheet,' (',table,')')
    subtitle = paste0('Conditional Enrichment (LogFC > 4, SS > 0.8)')
    plt <- ggbarplot(infile, bonf) + ggtitle(title, subtitle)
    print(plt)
  }
  graphics.off()
}

#d read_xlsx('derived/tables/210319_goa_bp_table_conditional_enrichment_logfc4.xlsx')


