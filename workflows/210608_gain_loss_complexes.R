# Explore what corum complexes are enriched in the data through 
# networks that are colored by heatmap expression. Compares WT versus MT
# starvartion condition only.


# load these libraries
library(readxl)
library(ggplot2)
library(genoppi)
library(ggrepel)
library(RColorBrewer)
devtools::load_all()

# data path
path <- 'inst/extdata/KRAS/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
in_sheets <- in_sheets[grepl('(starv)|(fcs)', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
#data <- set_saintscore_workbook(data)
#data <- filter_workbook(data, SaintScore > 0.8)
wt <- data[["WT_Starvation"]]
data[["WT_Starvation"]] <- NULL
data[["WT_FCS"]] <- NULL

# rename strange sheet name
names(data)[4] <- "G13_FCS"

# load up/down regulaton data
#reg <- as.data.frame(fread('derived/KRAS/tables/210422_kras_log_wt_mt_ss00_logfc0.csv'))
reg <- as.data.frame(fread('derived/KRAS/tables/210518_kras_log_wt_mt_ss00_logfc0.csv'))
colnames(reg)[2:7] <- c("G12_Starvation", "G13_Starvation", "Q61H_Starvation","G12_FCS",  "G13_FCS", "Q61H_FCS" )
reg$cluster_complete_order <- NULL
reg$Prey <- reg$wt_prey
reg$wt_prey <- NULL

# check input
stopifnot(colnames(reg)[!colnames(reg) %in% 'Prey'] %in% names(data))
stopifnot(names(data) %in% colnames(reg))

combine_data <- function(x1, x2, saintscore = 0.8, combine_by = c('Bait',"Prey")){

    # setup suffixes
    suffix1 <- paste0('.',unique(x1$Sheet))
    suffix2 <- paste0('.',unique(x2$Sheet))
    stopifnot(suffix1 != suffix2)
    
    # combine data
    mrg <- merge(x1, x2, all.x = T, all.y = T, by.x = combine_by, by.y = combine_by, suffixes = c(suffix1, suffix2))
    
    # set NA saintscore values to 0 and subset by saintscoring threshold
    mrg[[paste0('SaintScore', suffix1)]][is.na(mrg[[paste0('SaintScore', suffix1)]])] <- 0
    mrg[[paste0('SaintScore', suffix2)]][is.na(mrg[[paste0('SaintScore', suffix2)]])] <- 0
    #mrg[[paste0('FoldChange', suffix1)]][is.na(mrg[[paste0('FoldChange', suffix1)]])] <- 0
    #mrg[[paste0('FoldChange', suffix2)]][is.na(mrg[[paste0('FoldChange', suffix2)]])] <- 0
    mrg$significant <- mrg[[paste0('SaintScore', suffix1)]] >= saintscore | mrg[[paste0('SaintScore', suffix2)]] >= saintscore
    
    # clean up data
    mrg$gene <- mrg$Prey
    mrg <- mrg[,c('gene','significant', 
                  paste0('SaintScore', suffix1), paste0('SaintScore', suffix2),
                  paste0('FoldChange', suffix1), paste0('FoldChange', suffix2))]
    mrg <- mrg[!duplicated(mrg),]
    mrg <- mrg[! mrg$gene[!mrg$significant] %in% mrg$gene[mrg$significant], ]
    mrg <- mrg[rev(order(mrg$significant)),]
    return(mrg)
}

# Find out what complexes are enriched in the data (basically, 
# what complexes do want to display? This may take some time to run.)
complex_enrichment <- lapply(names(data), function(cur_sheet){
  
  # load data
  print(cur_sheet)
  cur_data <- combine_data(wt, data[[cur_sheet]], 0.8)
  cur_data <- cur_data[,c('gene','significant')]
  
  #cur_data <- data[[cur_sheet]]
  #cur_data <- data.frame(gene = cur_data$Prey, significant = cur_data$SaintScore >= 0.8)
  
  # setup for global analysis
  complexes <- unique(all_corum_complexes$complex)
  background <- data.frame(gene = unique(all_corum_complexes$genes), significant = F)
  cur_data_global <- rbind(cur_data, background[!background$gene %in% cur_data$gene,])
  
  my_enrichment <- lapply(complexes , function(complex_name){
    
    # prepare complexes
    cur_complex <- all_corum_complexes
    cur_complex$significant <- cur_complex$complex %in% complex_name
    cur_complex <- data.frame(gene = cur_complex$genes, significant = cur_complex$significant)
    
    # conditional enrichment
    hyper_cond <- calc_hyper(cur_data, cur_complex, intersectDf = data.frame(intersectN = F), bait = 'KRAS')
    out_cond <- hyper_cond$statistics
    out_cond$list_name <- complex_name
    out_cond$enrichment <- 'conditional'
    
    # global enrichment
    hyper_global <- calc_hyper(cur_data_global, cur_complex, intersectDf = data.frame(intersectN = F), bait = 'KRAS')
    out_global <- hyper_global$statistics
    out_global$list_name <- complex_name
    out_global$enrichment <- 'global'
    out <- rbind(out_cond, out_global)
    return(out)
    
  })
  
  # sort data 
  outdf <- do.call(rbind, my_enrichment)
  outdf <- outdf[(order(outdf$pvalue)),]
  outdf$FDR <- NA
  outdf$FDR[outdf$enrichment == 'global'] <- stats::p.adjust(outdf$pvalue[outdf$enrichment == 'global'], method = 'fdr')
  outdf$FDR[outdf$enrichment != 'global'] <- stats::p.adjust(outdf$pvalue[outdf$enrichment != 'global'], method = 'fdr')
  return(outdf)
})


names(complex_enrichment) <- names(data)
path <- 'derived/KRAS/tables/210603_corum_conditional_and_global_s8.xlsx'
#path <- 'derived/KRAS/tables/'
#path <- 'derived/KRAS/tables/210518_corum_conditional_and_global_s8.xlsx'
#write_xlsx(complex_enrichment, path)
#complex_enrichment <- read_excel_sheets(path)
#extracted_complexes <- lapply(complex_enrichment, function(x){ head(x$list_name, 40)})


outlist <- list()
pdf('derived/KRAS/plots/210603_corum_conditional_and_global_ss8_wtmt_combi.pdf', width = 12, height = 14)
for (name in names(complex_enrichment)){
  
  for (enrichment in c('conditional','global')){
   
    # subset data we want to plot
    saintscore_threshold <- 0.8
    cur_data <- combine_data(wt, data[[name]], saintscore_threshold)
    colnames(cur_data) <- gsub('_Averange','',colnames(cur_data))
    cur_data$Prey <- cur_data$gene
    cur_data$Bait <- 'KRAS'
    save_data <- cur_data
    cur_data <- cur_data[cur_data$significant,]
    cur_data <- cur_data[,c('Bait','Prey','gene')]
    
    # check overlap with reg
    #ints <- cur_data$gene[cur_data$significant]
    #print(paste(sum(ints %in% reg$Prey), '/', length(ints)))
    #print(paste(sum(ints %in% reg$wt_prey), '/', length(ints)))
    
    # get complexes data
    cur_complex <-  complex_enrichment[[name]]
    cur_complex <- cur_complex[cur_complex$enrichment %in% enrichment,]
    cur_complex <- head(cur_complex, 40)
    cur_data_complexes <- as.data.frame(all_corum_complexes[all_corum_complexes$complex %in% cur_complex$list_name & all_corum_complexes$organism == 'Human',])
    cur_complex_list <- lapply(unique(cur_data_complexes$complex), function(x) cur_data_complexes$genes[cur_data_complexes$complex %in% x])
    names(cur_complex_list) <- unique(cur_data_complexes$complex)
    
    # get regulation (colors) for nodes
    cur_reg <- reg[,colnames(reg) %in% c(name,"Prey")]
    colnames(cur_reg) <- c('val', 'Prey')
    values <- cur_reg[cur_reg$Prey %in% cur_data$Prey & cur_reg$Prey %in% unlist(cur_complex_list),]$val
    min_color <- -3 #min(values)*1.2
    max_color <- 3 #max(values)*1.2
    colorize <- function(x) assign_color_scale(x, colors = c('yellow','white','red'), force_min =  min_color, force_max = max_color, force_length_out = 400)
    node_colors <- colorize(cur_reg)
    node_colors <-  node_colors[!is.na(node_colors$color),]
    
    # get regulation (colors) for hulls
    convex_groups <- lapply(names(cur_complex_list), function(xname) {
      x <- cur_complex_list[[xname]]
      prey_in_data <- x[x %in% cur_data$Prey]
      group <- colorize(data.frame(Prey = xname, val = mean(cur_reg$val[cur_reg$Prey %in% prey_in_data])))
      return(group$color)
    })
    names(convex_groups) <- names(cur_complex_list)
    
    # convert to list 
    prey_unique <- unique(node_colors$Prey)
    cur_col_list <- lapply( prey_unique , function(x) node_colors$color[node_colors$Prey %in% x] )
    names(cur_col_list) <- prey_unique
    
    # titles
    title = paste(name, '+ WT_Starvation \n CORUM Enrichment Threshold: SaintScore > 0.8 +',toupper(enrichment))
    area_force = 1.9 #2.3 #2
    repulse_rad =2.9 #2.9
    
    # setup legend
    legend <- suppressWarnings(setup_graph_legend(min_color, max_color, colorize, text = paste(name,'Log(WT/MT)')))
    print(legend)
    
    # setup normal plot 
    make_group_graph(cur_data, cur_complex_list,
                      simple_plot = T, node_colors = cur_col_list,
                      size_group = 9, size_not_group =  4, 
                      area_force = area_force, repulse_rad = repulse_rad , 
                      plot_labels = "NONE", main = title)
    
    # setup plot with node labels
    make_group_graph(cur_data, cur_complex_list, 
                      simple_plot = T, node_colors = cur_col_list,
                      size_group = 9, size_not_group =  4, 
                      area_force = area_force, repulse_rad = repulse_rad , 
                      plot_labels = 'NODES', main = title)

    # setup plots grouped by complexes
    make_group_graph(cur_data, cur_complex_list, convex_groups,
                     convex_hull_borders = adjustcolor('black', alpha.f = 0.03),
                     simple_plot = T, node_colors = cur_col_list,
                     size_group = 9, size_not_group =  4, 
                     area_force = area_force, repulse_rad = repulse_rad , 
                     plot_labels = 'NONE', main = title)
    
    # setup plots grouped by complexes
    make_group_graph(cur_data, cur_complex_list, convex_groups,
                      convex_hull_borders = adjustcolor('black', alpha.f = 0.03),
                      simple_plot = T, node_colors = cur_col_list,
                      size_group = 9, size_not_group =  4, 
                      area_force = area_force, repulse_rad = repulse_rad , 
                      plot_labels = 'NODES', main = title)
    
    # count complexes
    count <- lapply(cur_complex_list,function(x) paste0(sum(x %in% cur_data$Prey),'/',length(x)))
    cur_complex_count_list <- cur_complex_list
    names(cur_complex_count_list) <- paste(names(cur_complex_count_list), count)
    
    # finally plot with group names + count
    make_group_graph(cur_data, cur_complex_count_list, simple_plot = T,
                      node_colors = cur_col_list,
                      size_group = 9, size_not_group =  4, 
                      area_force = area_force, repulse_rad = repulse_rad , 
                      plot_labels = T, main = title)
    
    
    ##################
    # writing files #
    #################
    
    # combine with reg to save
    save_reg <- reg[,c('Prey',name)]
    colnames(save_reg)[2] <- paste0(colnames(save_reg)[2],'_LogFC-WT_Starvation_LogFC')
    out_data <- merge(save_data, save_reg, by = 'Prey')
    out_data$gene <- NULL
    
    # merge with complexes
    complex_mrg <- stack(cur_complex_list)
    colnames(complex_mrg) <- c('Prey','Complex')
    complex_count <- stack(count)
    colnames(complex_count) <- c('Count','Complex')
    out_data_complex <- merge(out_data, merge(complex_count, complex_mrg))
    out_data_complex <- out_data_complex[order(out_data_complex$Complex),]
   
    # merge with complex up/down regulation
    convex_groups_val <- lapply(names(cur_complex_list), function(xname) {
      x <- cur_complex_list[[xname]]
      prey_in_data <- x[x %in% cur_data$Prey]
      group <- colorize(data.frame(Prey = xname, val = mean(cur_reg$val[cur_reg$Prey %in% prey_in_data])))
      return(group$val)
    })
    names(convex_groups_val) <- names(cur_complex_list)
    convex_groups_val_df <- stack(convex_groups_val)
    colnames(convex_groups_val_df) <- c('mean LogFC(MT/WT)','Complex')
    out_data_complex_val <- merge(out_data_complex, convex_groups_val_df) 
    out_data_complex_val$SaintScoreThreshold <- saintscore_threshold
    colnames(out_data_complex_val)
    
    
    # re-organize columns
    out_data_complex_val <- out_data_complex_val[,c('Bait','Prey','significant', 'SaintScoreThreshold',
                            paste0('SaintScore.',name), 'SaintScore.WT_Starvation',
                            paste0('FoldChange.',name), 'FoldChange.WT_Starvation',
                            paste0(name,'_LogFC-WT_Starvation_LogFC'),
                            'mean LogFC(MT/WT)',
                            'Count', 
                            'Complex'
                            )]
                            
    # save file to list
    outfile <- paste0(name,'_',enrichment)
    outlist[[outfile]] <- out_data_complex_val
    
                      
  }

}
graphics.off()

path <- 'derived/KRAS/tables/210603_wt_vs_mt_corum.xlsx'
write_xlsx(outlist, path)








