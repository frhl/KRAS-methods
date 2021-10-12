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

# Get throughput groups
int <- as.data.frame(read_xlsx('inst/extdata/KRAS/BIOGRID-GENE-110043-4.3.196.New.xlsx'))
int_ht <- int[int$`Experimental System` %in% 'Proximity Label-MS' & int$Throughput == 'High Throughput',] # not labelled
int_ht <- unique(c(int_ht$`Official Symbol Interactor A`, int_ht$`Official Symbol Interactor B`))
int_lt <- int[int$Throughput == 'Low Throughput',] 
int_lt <- unique(c(int_lt$`Official Symbol Interactor A`, int_lt$`Official Symbol Interactor B`))

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

plot_df[plot_df$wt_prey %in% 'MLLT4',]

# only keep proteins with specific saint score
out_res <- df_logfc_difference
res <- do.call(rbind, df_logfc_difference)

# what proteins should we keep in heatmap?
# keep_bool is for just removing saint score proteins,
# whereas keep_bool_2 is for removing groups of proteins below logfc diff thresholds.
res_ss <- res[res$mt_saintscore >= 0.8 | res$wt_saintscore >= 0.8,] 
unique_proteins <- unique(res_ss$wt_prey)
keep_bool <- unlist(lapply(unique_proteins, function(x) any(abs(res_ss$logfc_mt_wt[res_ss$wt_prey %in% x]) >= 0))) # <----
keep_bool_2 <- unlist(lapply(unique_proteins, function(x) any(abs(res_ss$logfc_mt_wt[res_ss$wt_prey %in% x]) >= 2)))
keep_bool_3 <- unlist(lapply(unique_proteins, function(x) any(abs(res_ss$logfc_mt_wt[res_ss$wt_prey %in% x]) >= 3)))
keep_protein <- unique_proteins[keep_bool]
keep_protein_2 <- unique_proteins[keep_bool_2]
keep_protein_3 <- unique_proteins[keep_bool_3]

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

# cluster the three groups seperately
lt_genes <- int_lt[!int_lt %in% effectors]
ht_genes <- int_ht[! (int_ht %in% int_lt | int_ht %in% effectors) & int_ht %in% keep_protein_2]
rest <- keep_protein[! (keep_protein %in% int_lt | keep_protein %in% effectors | keep_protein %in% int_ht) & keep_protein %in% keep_protein_3]
groups <- list(effectors = effectors, lt = lt_genes, ht = ht_genes, rest = rest)

count <- 0
hlines <- c()
df_wide_grouped <- do.call(rbind, lapply(groups, function(genes){
  
  df_new <- df_wide[df_wide$wt_prey %in% genes, 2:7]
  genes <- df_wide[df_wide$wt_prey %in% genes, 1]
  tree <- hclust(dist(df_new), method = 'complete')
  df_new$cluster_complete_order <- tree$order + count
  df_new$wt_prey <- genes
  count <<- count + length(genes)
  hlines <<- c(hlines, count)
  return(df_new)
  
}))

duplicated(df_wide_grouped$cluster_complete_order)
df_wide_grouped$group <- unlist(lapply(strsplit(rownames(df_wide_grouped), split = '\\.'), function(x) x[1]))

# set order to clustered order
plot_df <- df_res
plot_df <- plot_df[plot_df$wt_prey %in% df_wide_grouped$wt_prey,]
plot_df$wt_prey <- factor(plot_df$wt_prey, levels = df_wide_grouped$wt_prey[df_wide_grouped$cluster_complete_order])
plot_df$grid <- 'Mutant Log2FC / Wildtype Log2FC'

# group order in plot
group_order <- df_wide_grouped$group[df_wide_grouped$cluster_complete_order]
print(group_order)
table(group_order)/sum(table(group_order))


sub_labels <- function(x) {
  x <- gsub('_Starvation','_STV',x)
  x <- gsub('_Averange','',x)
  return(x)
}


plot_df$experiment <- factor(sub_labels(plot_df$experiment), levels = sub_labels(conds))

# overlay with known interactors
# Setup colors
label_colors <- ifelse(df_wide$wt_prey %in% effectors, 'red', ifelse(df_wide$wt_prey %in% hc_interactors, 'blue','black'))
label_colors <- label_colors[tree$order]



# plot data
stopifnot(c('LZTR1','ARAF','RAF1') %in% df_res$wt_prey)
p1 <- ggplot(plot_df, aes(x = experiment, y = wt_prey, fill = logfc_mt_wt)) +
  geom_tile() +
  theme_bw() + 
  xlab('Experiment') + 
  ylab('Prey protein') + 
  geom_hline(yintercept=hlines + 0.5, linetype = 'dashed') +
  geom_vline(xintercept=3.5) +
  ggtitle('KRAS Heatmap ["complete" clustering]', 'SaintScore > 0.8 & Log difference > 3 for either WT or MT condition') + 
  scale_fill_gradient2(low = 'navyblue', mid = 'white', high = 'firebrick1', midpoint = 0) +
  labs(fill = expression(Log[2]*FC[mt/wt])) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color = label_colors)
    #panel.border = element_blank()
  ) +
  #guides(fill=guide_legend(title=expression(Log[2]*FC[wt/mt]))) +
  facet_grid(~grid)


for (w in 3:6){
  ggsave(p1, filename = paste0('derived/KRAS/sizes/210803_kras_heatmap_grouped_wt_mt_logfc_w',w,'.pdf'), width = w, height = 18)
}



## also write file with all hits
