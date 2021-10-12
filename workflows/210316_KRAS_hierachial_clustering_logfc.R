library(devtools)
devtools::load_all()

## in files and helper functions

# use tidyverse package to read excel files
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# read in the data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data, SaintScore >= 0.65 & LogFC > 4)
data <- collapse_workbook(data)

# effector genes 
effector_genes <- read.table('inst/extdata/effector_genes.txt', header = T)$gene
data <- data[data$Prey %in% effector_genes,]

# We are only interested in LogFC
logfc_data_pre <- data[,c('Prey','Sheet','LogFC')]
logfc_data_wide <- complete_missing(logfc_data_pre, return_as = 'wide')
logfc_data <- complete_missing(logfc_data_pre, return_as = 'long')

# format for heatmap
row_names <- logfc_data_wide$Prey
logfc_data_wide$Prey <- NULL
logfc_data_wide <- as.data.frame(sapply(logfc_data_wide, as.numeric))
rownames(logfc_data_wide) <- row_names
colnames(logfc_data_wide) <- gsub('LogFC\\.','',colnames(logfc_data_wide))

# setup heatmap
labels <- round(as.matrix(logfc_data_wide), 2)
palette = colorRampPalette(c('white', 'orange'))(n = 200)
pdf('derived/210317_heatmap_kras_effector_genes_saintscore065.pdf', width = 12, height = 10)
gplots::heatmap.2(as.matrix(logfc_data_wide),
                  Colv=FALSE,
                  dendrogram = 'row',
                  cellnote = labels,    # same data set for cell labels
                  #main = bait,          # heat map title
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  margins =c(12,9),     # widens margins around plot
                  col=palette)       # use on color palette defined earlier
graphics.off()

# cluster data
tree <- hclust( dist(logfc_data$LogFC))
logfc_data$Prey <- factor(logfc_data$Prey, levels = unique(logfc_data$Prey[tree$order]))
logfc_data$Sheet <- factor(logfc_data$Sheet, levels = unique(logfc_data$Sheet))
logfc_data$label <- round(logfc_data$LogFC,1)

# heatmap ggplot
plt <- ggplot(logfc_data, aes(x = Sheet, y = Prey, fill = LogFC, label = label)) +
  geom_tile(color = 'black') + 
  xlab('Experiment') + 
  ylab('Prey') + 
  geom_text() + 
  theme_minimal() + 
  ggtitle('LogFC for selected effector protein preys of K-RAS ') + 
  scale_fill_gradient(low = 'white', high = 'orange') +
  theme(axis.text.x = element_text(angle = 90))
plt

ggsave('derived/210316_heatmap_kras_effector_genes_saintscore065.png', width = 8, height = 6)













