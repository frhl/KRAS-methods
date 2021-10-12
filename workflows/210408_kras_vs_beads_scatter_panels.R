# what genesets are enriched in the WT versus MT.

library(readxl)
library(ggplot2)
library(genoppi)
library(ggrepel)
library(RColorBrewer)
devtools::load_all()

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
in_sheets <- in_sheets[grepl('(starv)|(fcs)', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
d <- do.call(rbind, filter_workbook(data))

# effector proteins
effectors <- fread('inst/extdata/effector_genes.txt')$gene
effectors <- effectors[! effectors %in% 'LZTR1']
effectors <- effectors[! effectors %in% '']
interactors <- get_inweb_list('KRAS')
interactors <- interactors$gene[interactors$significant]
hc_interactors <- fread('inst/extdata/hc_ehc_kras_interactors.csv')
hc_interactors <- hc_interactors$EHC[hc_interactors$EHC != ""]

# setup types
d$type <- 'N/A'
d$type[d$Prey %in% interactors] <- 'InWeb_InBioMap Interactor'
d$type[d$Prey %in% hc_interactors] <- 'High Confidence Interactors'
d$type[d$Prey %in% effectors] <- 'Effector Protein' 
d$type <- as.factor(d$type)

# colors
my_colors <- c('red','blue','yellow','black') #brewer.pal(3,"Set1")
names(my_colors) <- levels(d$type)
col_scale <- scale_colour_manual(name = "type",values = my_colors)

# plot data
#d_background <- d[d$type %in% 'N/A',]
#p1 <- ggplot(d_background, aes(x = SaintScore, y = LogFC, color = type, label = Prey)) +
#  geom_point(alpha = 0.3) +
#  col_scale +
#  geom_hline(yintercept=4, linetype = 'dashed', alpha = 0.5) +
#  geom_vline(xintercept=0.8,  linetype = 'dashed', alpha = 0.5) +
#  theme_bw() +
#  ggtitle('KRAS SaintScore Versus LogFC') + 
#  facet_wrap(~Sheet)

# subset data for labels
#d_interactor <- d[d$type %in% 'InWeb_InBioMap Interactor',]
#d_hc_interactor <- d[d$type %in% 'High Confidence Interactors',]
#d_effector <- d[d$type %in% 'Effector Protein',]
#p1 <- p1 + 
#  geom_point(data = d_interactor, size = 2) +
#  geom_point(data = d_hc_interactor, size = 2) +
#  geom_point(data = d_effector,  size = 2)


#pdf('derived/plots/210607_logfc_vs_saintscore_scatterplot.pdf', width = 16, height = 18)
#p1 + geom_text_repel(data = d_effector, color = 'black', size = 2, max.overlaps = 40)
#p1 + geom_text_repel(data = d_interactor, color = 'black', size = 2, max.overlaps = 40)
#p1 + geom_text_repel(data = d_hc_interactor, color = 'black', size = 2, max.overlaps = 40)
#graphics.off()

# setup names
d$row <- as.factor(ifelse(grepl('Starv',d$Sheet),'Starvation','FCS'))
d$experiment <- unlist(lapply(strsplit(d$Sheet, split = '_'), function(x) x[1]))
d$experiment <- factor(d$experiment, levels = unique(d$experiment))

## make jitter plots
d2 <- d[d$SaintScore >= 0.8,]
d2$unif <- runif(nrow(d2))
d2_background <- d2[d2$type %in% 'N/A',]

# calculate p-values for overlap
stat_df <- data.frame(gene = d2$Prey, significant = ifelse(d2$LogFC > 4, T, F), Sheet = d2$Sheet)
lst <- do.call(rbind, lapply(unique(stat_df$Sheet), function(cur_sheet){
  
  # calculate stats
  hyper_effectors <- calc_hyper(stat_df[stat_df$Sheet == cur_sheet,], data.frame(gene = effectors, significant = T), intersectDf = data.frame(intersectN = F), bait= 'KRAS')
  hyper_hc <- calc_hyper(stat_df[stat_df$Sheet == cur_sheet,], data.frame(gene = hc_interactors, significant = T), intersectDf = data.frame(intersectN = F), bait= 'KRAS')
  
  row <- ifelse(grepl('Starv',cur_sheet),'Starvation','FCS')
  experiment <- unlist(lapply(strsplit(cur_sheet, split = '_'), function(x) x[1]))
  
  # rename data
  hyper_effectors <- hyper_effectors$statistics
  hyper_effectors$list_name <- 'Effector Protein'
  hyper_effectors$Sheet <- cur_sheet
  hyper_effectors$row <- row
  hyper_effectors$experiment <- experiment
  hyper_hc <- hyper_hc$statistics
  hyper_hc$list_name <- 'High Confidence Interactors'
  hyper_hc$Sheet <- cur_sheet
  hyper_hc$row <- row
  hyper_hc$experiment <- experiment
  
  return(rbind(hyper_hc, hyper_effectors))
}))

# subset statdf
lst <- lst[order(lst$list_name),]
lst$lab <- paste0("P-value (",lst$list_name,"): ",round(lst$pvalue, 4))
lst$type <- "X"
lst$x <- 0
lst$y <- Inf
lst$row <- factor(lst$row, levels = c('Starvation','FCS'))
lst$experiment <- factor(lst$experiment, levels = c("WT","G12", "G13","Q61H"))

# plot data
p2 <- ggplot(d2_background, aes(x = unif,y = LogFC, color = type, label = Prey)) +
  geom_point(alpha = 0.9, color = 'grey') +
  geom_hline(yintercept=4, linetype = 'dashed', alpha = 0.5) +
  col_scale +
  theme_bw() +
  ggtitle('KRAS SaintScore Versus LogFC jitterplot', 'SaintScore >= 0.8') + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_grid(row~experiment) 

# for each analysis calculate p-value
d2_hc_interactor <- d2[d2$type %in% 'High Confidence Interactors',]
d2_effector <- d2[d2$type %in% 'Effector Protein',]

# plot the data
p2 <- p2 + 
  geom_point(data = d2_hc_interactor, size = 2) +
  geom_point(data = d2_effector,  size = 2) 

pdf('derived/plots/210412_logfc_vs_saintscore_jitter.pdf', width = 10, height = 8)
p2 
p2 + geom_text_repel(data = d2_effector, color = 'black', size = 2, max.overlaps = 40) +
  geom_text(aes(x, y, label=lab), data = lst[lst$list_name %in% 'Effector Protein',], color = 'black', size = 2, vjust = 1.5, hjust = 0, alpha = 0.5)
p2 + geom_text_repel(data = d2_hc_interactor, color = 'black', size = 2, max.overlaps = 40) +
  geom_text(aes(x, y, label=lab), data = lst[lst$list_name %in% 'High Confidence Interactors',], color = 'black', size = 2, vjust = 1.5, hjust = 0, alpha = 0.5)
graphics.off()

# gradually display data
#by = 0.25
#sequence <- seq(min(d2$LogFC), max(d2$LogFC), by = by)
#for (i in sequence){
#  start = i
#  stop = i+by
#  d_labels <- d2[d2$LogFC >= start  & d2$LogFC < stop,]
#  print(nrow(d_labels))
#  p3 <- p2 + geom_text_repel(data = d_labels, color = 'black', size = 2, max.overlaps = 40)
#  print(p3)
#}
#graphics.off()


