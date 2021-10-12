devtools::load_all()
library(data.table)

# prepare files
files <- list.files('derived/KRAS/tables/211010/', pattern = '.csv', full.names = TRUE)
newnames <- tools::file_path_sans_ext(basename(files))
groups <- unique(unlist(lapply(strsplit(newnames, split = '_'), function(x) paste0(x[4:length(x)], collapse = '_'))))
condtions <- c('FCS','STV')
analyses <- c('global','conditional')

# load files
dts <- lapply(files, fread)
names(dts) <- newnames

# load genoppi names
gtex_map <- fread('~/Projects/08_genesets/genesets/data/gtex/GTEX.tstat.categories.genoppi.csv')

# setup out dirs
out_dir <- 'derived/KRAS/plots/211011_GTEx_heatmaps'
dir.create(out_dir)

# for each group iterate ocer the repsective files
for (group in groups){
  
  dts_groups <- dts[grepl(group, newnames)]
  
  lapply(analyses, function(analysis){
    
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
    
    # map to prettier labels
    dts_combined$Tissue.genoppi <- dts_combined$list_name
    dts_combined <- merge(dts_combined, gtex_map, all.x = TRUE)
    
    # combine and re-order 
    dts_ordered <- dts_combined[order(dts_combined$pvalue),]
    #dts_ordered$list_name <- factor(dts_ordered$list_name, levels = rev(unique(dts_ordered$list_name)))
    
    # P-value for boneferroni correction
    n_pathways <- length(unique(dts_ordered$list_name))
    p_critical <- 0.05 / n_pathways
    dts_ordered$significant <- dts_ordered$pvalue < p_critical
    dts_ordered$label <- ifelse(dts_ordered$significant, '*','')
    
    # make title
    side <- unique(dts_ordered$side)
    dataset <- unique(dts_ordered$dataset)
    test <- analysis
    condition <- 'STV+FCS'
    analysis_title <- paste(dataset, side, paste0('(',test,' analysis)'))
    analysis_subtitle = '(*) indicates significant at alpha = 0.05 with bonferroni correction'
    
    # setup plotting
    plt <- ggplot(dts_ordered, aes(x = name_minimal, y = Tissue, fill = -log10(pvalue), label = label)) +
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
    ggsave(plot = plt, filename = out, width = 9, height = 10)
    
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





