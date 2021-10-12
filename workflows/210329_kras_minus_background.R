# analysis for Andreas Phizer meeting
# - only starvation conditions are considered.
# - using WT as background proteins

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
data <- filter_workbook(data, SaintScore >= 0.8)
wt <- data$WT_Starvation
mt <- data[2:4]

# set background to prey of WT (starvation)
background <- unique(wt$Prey)

# MT minus WT
mt_minus_wt <- lapply(mt, function(x){
  unique(x$Prey[! x$Prey %in% background])
})

# MT overlap
mt_minus_wt_intersect <- Reduce(intersect, mt_minus_wt)
mt_minus_wt_union <- Reduce(union, mt_minus_wt)

# 1) generate heatmap of logFC
mt_logfcs <- do.call(rbind, mt)[,c(2,3,6)]
mt_logfcs_intersect <- mt_logfcs[mt_logfcs$Prey %in% mt_minus_wt_intersect,]
mt_logfcs_union <- mt_logfcs[mt_logfcs$Prey %in% mt_minus_wt_union,]
#clust <- hclust(dist(mt_logfcs$LogFC), method = 'complete')
#plot(clust, labels = mt_logfcs$Prey)
#mt_logfcs$order <- clust$order
p1 <- ggplot(mt_logfcs_intersect, aes(x=Sheet, y = reorder(Prey, LogFC), fill = LogFC, label = round(LogFC, 2))) +
  geom_tile(color = 'black') + 
  geom_text(size = 2) + 
  xlab('Experiment') +
  ylab(expression(Log[2]*" FC")) + 
  ggtitle('MT (intersect) minus WT (All SaintScore >= 0.8)', 'Ordered by LogFC') +
  theme_bw() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12,face="bold"))

p2 <- ggplot(mt_logfcs_union, aes(x=Sheet, y = reorder(Prey, LogFC), fill = LogFC, label = round(LogFC, 2))) +
  geom_tile(color = 'black') + 
  geom_text(size = 2) + 
  xlab('Experiment') +
  ylab(expression(Log[2]*" FC")) + 
  ggtitle('MT (union) minus WT (All SaintScore >= 0.8)', 'Ordered by LogFC') +
  theme_bw() + 
  theme(axis.text.y = element_text(size=6), axis.title=element_text(size=12,face="bold"))

#ggsave('derived/plots/210329_mt_intersect_minus_wt_ss08_heatmap.png', p1, width = 5, height = 14)
#ggsave('derived/plots/210329_mt_union_minus_wt_ss08_heatmap.png', p2, width = 5, height = 20)

# 2) generate network plot for these data points.

# combine origin data into shared interactors of mutans versus WT background
origin_shared <- data.frame(gene = unique(mt_logfcs_union$Prey), significant = TRUE)
origin_background <- data.frame(gene = background, significant = FALSE)
origin <- rbind(origin_shared, origin_background)

# get overlap with GO data bases
go_db <- list(goa_bp = goa_bp_table, goa_cc = goa_cc_table, goa_mf = goa_mf_table)
mt_logfc_geneset <- lapply(go_db, function(db){
  
  # geneset reference 
  reference <- data.frame(gene = db$Gene.symbol, pathway = db$GO.name, significant = TRUE)
  go_overlap <- lapply_calc_hyper(origin, reference, col.by = 'pathway', bait = "KRAS", intersectN = F, verbose = T)
  return(go_overlap)
  
})

#pdf('derived/plots/210329_mt_union_minus_wt_ss08_network_go.pdf', width = 10, height = 10)
# aggregate data
for (p in names(mt_logfc_geneset)){
  
  # go pathways
  go <- mt_logfc_geneset[[p]]
  go <- go[rev(order(go$successInSample_count)),]
  go <- go[1:10,]
  geneset <- strsplit(gsub(' ','', go$successInSampleGenes), split = ';')
  names(geneset) <- go$list_name
  
  # plot data
  inweb <- get_inweb_list('KRAS')
  inweb <- inweb$gene[inweb$significant]
  df <- data.frame(origin = "KRAS", target = origin$gene[origin$significant])
  
  # titles
  title = paste0("MT (union) - WT",'\n(',p,' ordered by top 10 most recurring genes in pathways, SaintScore >= 0.8)')
  
  # plot data
  setup_group_graph(df, geneset, known_interactors = inweb, main = title, size.group = 5)
  
}
#graphics.off()


