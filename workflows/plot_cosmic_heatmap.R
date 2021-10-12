devtools::load_all()
library(data.table)

# prepare files
files <- list.files('derived/tables/211012_pancreas_test/', pattern = '.csv', full.names = TRUE)
#files <- files[grepl("(hallmark)|(celltypes)", files)]
newnames <- tools::file_path_sans_ext(basename(files))
groups <- unique(unlist(lapply(strsplit(newnames, split = '_'), function(x) paste0(x[4:length(x)], collapse = '_'))))
condtions <- c('FCS','STV')
analyses <- c('global','conditional')

# load files
dts <- lapply(files, fread)
names(dts) <- newnames

# setup out dirs
out_dir <- 'derived/plots/211012_cosmic_heatmaps'
dir.create(out_dir)

# for each group iterate ocer the repsective files
for (group in groups){
  
  print(group)
  dts_groups <- dts[grepl(group, newnames)]
  
  lapply(analyses, function(analysis){
    print(analysis)
    
    # combine groups
    dts_combined <- do.call(rbind, lapply(names(dts_groups), function(dt_name){
      dt <- dts_groups[[dt_name]]
      dt$name_minimal <- unlist(lapply(strsplit(dt_name, split = '_'), function(x) paste0(x[1:2], collapse = '_')))
      dt$condition <- unlist(lapply(strsplit(dt_name, split = '_'), function(x) x[3]))
      dt$side <-  unlist(lapply(strsplit(dt_name, split = '_'), function(x) x[4]))
      dt$dataset <-  unlist(lapply(strsplit(dt_name, split = '_'), function(x) paste0(x[5:length(x)], collapse = '_')))
      dt$experiment <- dt_name
      return(dt)
    }))
    
    # subset to either conditional or global enrichment analysis
    bool_analysis <- dts_combined$analysis == analysis
    dts_combined <- dts_combined[bool_analysis]
    
    # combine and re-order 
    dts_ordered <- dts_combined[order(dts_combined$pvalue),]
    dts_ordered$list_name <- factor(dts_ordered$list_name, levels = rev(unique(dts_ordered$list_name)))
    
    # P-value for boneferroni correction
    n_pathways <- length(unique(dts_ordered$list_name))
    p_critical <- 0.05 / n_pathways
    dts_ordered$significant <- dts_ordered$pvalue < p_critical
    #dts_ordered$significant <- dts_ordered$FDR < 0.1
    dts_ordered$label <- ifelse(dts_ordered$significant, '*','')
    
    # count up untill we have selected 20 pathways
    n_found = 0
    counter = 0
    while (n_found <= 19){
      pathways <- head(dts_ordered$list_name, n = counter)
      bool_pathways <- dts_ordered$list_name %in% pathways
      dts_subsetted <- dts_ordered[bool_pathways]
      n_found <- length(unique(dts_subsetted$list_name))
      counter = counter + 1
    }
    
    # make title
    side <- unique(dts_subsetted$side)
    dataset <- unique(dts_subsetted$dataset)
    test <- analysis
    condition <- 'STV+FCS'
    analysis_title <- paste(dataset, side, paste0('(',test,' analysis)'))
    analysis_subtitle = '(*) Indicates bonferroni significant'
    
    # setup plotting
    plt <- ggplot(dts_subsetted, aes(x = name_minimal, y = list_name, fill = -log10(pvalue), label = label)) +
      geom_tile() +
      geom_text() +
      scale_fill_gradient(low = 'white', high = 'red') + 
      facet_wrap(~condition) +
      theme_bw() + 
      xlab('experiment') + 
      ylab('GO term') +
      ggtitle(analysis_title, analysis_subtitle)
    
    # setup outpaths
    prefix <- paste(dataset, side, test, sep = '_')
    prefix_ext <- paste(prefix, '.pdf', sep = '')
    out <- file.path(out_dir, prefix_ext)
    write(paste0('writing ', out), stdout())
    ggsave(plot = plt, filename = out, width = 10, height = 6)
    
  })
  
}




#for (d in dts){
#  for (analysis in c('conditional','global')){
#    d_analysis <- d[d$analysis == 'global']
#    selected_pathways <- head(d_analysis$list_name, n_top)
#    d_pathways <- d_analysis[d_analysis$list_name %in% selected_pathways]
#    
#  }
#}





