# data path
path <- 'inst/extdata/KRAS/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
#in_sheets <- in_sheets[grepl('starv', tolower(in_sheets))] # 

# get high confidence KRAS interactors
#hc_int <- fread('inst/extdata/210317_kras_hitproteins.txt')

# effector proteins
effectors <- fread('inst/extdata/KRAS/effector_genes.txt')$gene
effectors <- effectors[! effectors %in% 'LZTR1']
effectors <- effectors[! effectors %in% '']
interactors <- get_inweb_list('KRAS')
interactors <- interactors$gene[interactors$significant]
hc_interactors <- fread('inst/extdata/KRAS/hc_ehc_kras_interactors.csv')
hc_interactors <- hc_interactors$EHC[hc_interactors$EHC != ""]

# Get throughput
int <- as.data.frame(read_xlsx('inst/extdata/KRAS/BIOGRID-GENE-110043-4.3.196.tab3.xlsx'))
int_ht <- int[int$`Experimental System` %in% 'Proximity Label-MS' & int$Throughput == 'High Throughput',] # not labelled
int_ht <- unique(c(int_ht$`Official Symbol Interactor A`, int_ht$`Official Symbol Interactor B`))

# get proteomics data
conds <- c("G12_Starvation", "G13_Starvation", "Q61H_Starvation","WT_FCS", "G12_FCS", "G13_FCS_Averange","Q61H_FCS")
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
wt <- data$WT_Starvation
mts <- data[conds]


df_logfc_difference <- lapply(names(mts), function(mt_name){
  
  mt = mts[[mt_name]]
  
  # merge data by prey
  wt_df <- data.frame(wt_prey = wt$Prey, wt_FC = wt$FoldChange, wt_fdr = wt$BFDR, wt_saintscore = wt$SaintScore)
  mt_df <- data.frame(mt_prey = mt$Prey, mt_FC = mt$FoldChange, mt_fdr = mt$BFDR, mt_saintscore = mt$SaintScore)
  df <- merge(wt_df, mt_df, by.x = 'wt_prey', by.y = 'mt_prey', all.x = T, all.y = T)
  df$cell1 <- 'WT'
  df$cond1 <- 'Starvation'
  df$cell2 <- index(mt_name, 1)
  df$cond2 <- index(mt_name, 2)
                    
  
  # NA.rows are set to zero
  df$mt_FC[is.na(df$mt_FC)] <- 0
  df$wt_FC[is.na(df$wt_FC)] <- 0
  df$mt_saintscore[is.na(df$mt_saintscore)] <- 0
  df$wt_saintscore[is.na(df$wt_saintscore)] <- 0
  df$logfc_mt_wt <- log2( df$mt_FC / df$wt_FC )
  df$experiment <- mt_name
  return(df)
  
})


# only keep proteins with specific saint score
out_res <- df_logfc_difference
res <- do.call(rbind, df_logfc_difference)

# what proteins should we keep in heatmap?
res_ss <- res[res$mt_saintscore >= 0.8 | res$wt_saintscore >= 0.8,] 
unique_proteins <- unique(res_ss$wt_prey)
keep_bool <- unlist(lapply(unique_proteins, function(x) any(abs(res_ss$logfc_mt_wt[res_ss$wt_prey %in% x]) >= 2))) # <----
keep_protein <- unique_proteins[keep_bool]

# subset heatmap
df_res <- res[res$wt_prey %in% keep_protein,]
df_res <- na.omit(df_res)

# convert to wide for clustering
df_wide_pre <- df_res[,c("logfc_mt_wt", "wt_prey", "experiment")]
df_wide <- reshape(df_wide_pre, 
                   direction = "wide",
                   v.names = "logfc_mt_wt",
                   idvar = "wt_prey",
                   timevar = "experiment")
df_wide[is.na(df_wide)] <- 0
#tree <- hclust(dist(df_wide[,2:4]), method = 'complete')
tree <- hclust(dist(df_wide[,2:7]), method = 'complete')
df_wide$cluster_complete_order <- tree$order
nrow(df_wide)

# current saved data is agnostic to logfc. Remember to change above if overwriting!
#write.csv(df_wide, 'derived/KRAS/tables/210518_kras_log_wt_mt_ss00_logfc0.csv', row.names = F, quote = F) 
# save data
#writexl::write_xlsx(df_wide[tree$order,], path = 'derived/tables/210406_wt_versus_mt_heatmap_values.xlsx')

# set order to clustered order
plot_df <- df_res
plot_df$wt_prey <- factor(plot_df$wt_prey, levels = df_wide$wt_prey[tree$order])
plot_df$grid <- 'Wildtype Log2FC (Starvation) / Mutant Log2FC (Starvation / FCS)'
plot_df$experiment <- factor(plot_df$experiment, levels = conds)

# overlay with known interactors
# Setup colors
label_colors <- ifelse(df_wide$wt_prey %in% effectors, 'red', ifelse(df_wide$wt_prey %in% hc_interactors, 'blue','black'))
label_colors <- label_colors[tree$order]
# colors

# play around with zscores
#plot_df$zscore <- plot_df$logfc_mt_wt/sd(plot_df$logfc_mt_wt[is.finite(plot_df$logfc_mt_wt)])
#plot_df <- plot_df[abs(plot_df$zscore) > 1.96,]

#is_interactor <- levels(plot_df$wt_prey) %in% hc_int$Name
#label_colors <- ifelse(is_interactor, 'red', 'black')

# plot data
p1 <- ggplot(plot_df, aes(x = experiment, y = wt_prey, fill = logfc_mt_wt)) +
  geom_tile() +
  theme_bw() + 
  xlab('Experiment') + 
  ylab('Prey protein') + 
  geom_vline(xintercept=3.5) +
  ggtitle('KRAS Heatmap ["complete" clustering]', 'SaintScore > 0.8 & Log difference > 2 for either WT or MT condition') + 
  scale_fill_gradient2(low = 'navyblue', mid = 'white', high = 'firebrick1', midpoint = 0) +
  labs(fill = expression(Log[2]*FC[wt/mt])) +
  theme(axis.text.x = element_text(angle = 0)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = label_colors)
    #panel.border = element_blank()
  ) +
  facet_grid(~grid)
p1
#ggsave('derived/210518_kras_heatmap_wt_mt_logfc.pdf', width = 9, height = 32)
ggsave('derived/210803_kras_heatmap_wt_mt_logfc.pdf', width = 9, height = 32)


## save corresponding heatmap dendogram
pdf('derived/210412_kras_heatmap_dendogram_wt_mt_logfc.pdf', width = 20, height = 6)
labels <-  df_wide$wt_prey[tree$order]
labels[labels %in% hc_int$Name] <- paste0('known ---->',labels[labels %in% hc_int$Name])
plot(tree, 
     labels = labels, 
     cex = 0.5,
     main = 'hclustering of LogFC(WT/MT) across experiments')
graphics.off()



## generate out res
out_res_list <- lapply(df_logfc_difference, function(out_res){
  out_res$prey_type <- ifelse(out_res$wt_prey %in% effectors, 'Effector', ifelse(out_res$wt_prey %in% hc_interactors, 'High-confidence interactor',''))
  out_res <- out_res[,c(1,8:11,2:7,12,14)]
  out_res <- out_res[order(out_res$cell2, out_res$cond2, out_res$logfc_mt_wt),]
  return(out_res)
})
names(out_res_list) <- names(mts)

#write_xlsx(out_res_list, 'derived/tables/210419_heatmap_logfc_difference.xlsx')
