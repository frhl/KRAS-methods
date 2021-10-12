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
data <- set_saintscore_workbook(data)
d <- do.call(rbind, filter_workbook(data))

# what classes to keep
keep <- c(
  'activator',
  'Classical Ras Effector',
  'Proximitome',
  'receptor',
  'repressor'
)

# load classes to be plotted
mapping <- read_xlsx('inst/extdata/KRAS/20210826 KRAS interaction network-final.xlsx')
colnames(mapping) <- c('uniprot','Prey', 'class1', 'class2', 'class3')
mapping <- mapping[mapping$class1 %in% keep,]

# merge with original data
d_mrg <- merge(d, mapping, by = 'Prey', all.x = T)
#d_mrg$class1[is.na(d_mrg$class1)] <- 'New'
stopifnot(nrow(d_mrg) == nrow(d))
d <- d_mrg
d$type <- as.factor(d$class1)

# combine with deep proteome intensity
deep_proteome <- fread('derived/KRAS/tables/210408_deep_proteome_min2reps.csv')
colnames(deep_proteome) <- c("Prey",'LFQ')
d <- merge(d, deep_proteome, all.x = T)

# setup names
d$row <- as.factor(ifelse(grepl('Starv',d$Sheet),'Starvation','FCS'))
d$experiment <- unlist(lapply(strsplit(d$Sheet, split = '_'), function(x) x[1]))
d$experiment <- factor(d$experiment, levels = c("WT","G12", "G13","Q61H"))

# colors
#my_colors <- brewer.pal(6,"Set1")
my_colors <- c('blue','red','purple','green','cyan')
#names(my_colors) <- levels(d$type)
col_scale <- scale_colour_manual(name = "type",values = my_colors)
fill_scale <- scale_fill_manual(name = "type",values = my_colors)

## make jitter plots
d2 <- d[d$SaintScore >= 0.8,]
d2$LFQ[is.na(d2$LFQ)] <- 4 + runif(sum(is.na(d2$LFQ)))
#d2_background <- d2[d2$type %in% 'New',]
d2_background <- d2[is.na(d2$type),]

# some annotations
subtitle = 'Deep proteome median LFQ (At least 2 replicates with non-zero LFQ values)'


# plot data
p2 <- ggplot(d2_background, aes(x = LFQ,y = LogFC, color = type, label = Prey)) +
  geom_point(alpha = 1, color = 'grey', size = 2) +
  geom_hline(yintercept=4, linetype = 'dashed', alpha = 0.5) +
  scale_x_continuous(breaks = 6:max(ceiling(d2$LFQ))) +
  theme_bw() +
  xlab('Deep Proteome LFQ') +
  ggtitle('KRAS Scatter Plot of SaintScore versus deep proteome LFQ', subtitle) + 
  facet_grid(row ~ experiment)

# for each analysis calculate p-value
types <- na.omit(unique(d2$type))
type_list <- lapply(types, function(x) d2[d2$type %in% x,])
names(type_list) <- types

p2 <- p2 + geom_point(data = type_list$Proximitome, size = 2.5, fill = 'steelblue1', color = 'black', shape = 21)
p2 <- p2 + geom_point(data = type_list$activator, size = 2.5, fill = 'red', color = 'black', shape = 24)
p2 <- p2 + geom_point(data = type_list$repressor, size = 2.5, fill = 'red', color = 'black', shape = 25)
p2 <- p2 + geom_point(data = type_list$receptor, size = 2.5, fill = 'red', color = 'black', shape = 22)
p2 <- p2 + geom_point(data = type_list$`Classical Ras Effector`, size = 2.5, fill = 'red', color = 'black', shape = 21)

# plot
pdf('derived/KRAS/210921_deep_proteome_apex/210921_saintscore_vs_deep_proteome.pdf', width = 11, height = 8)
p2
p2 + geom_text_repel(data = type_list$Proximitome, color = 'black', size = 2, max.overlaps = 40)
p2 + geom_text_repel(data = type_list$activator, color = 'black', size = 2, max.overlaps = 40)
p2 + geom_text_repel(data = type_list$repressor, color = 'black', size = 2, max.overlaps = 40)
p2 + geom_text_repel(data = type_list$receptor, color = 'black', size = 2, max.overlaps = 40)
p2 + geom_text_repel(data = type_list$`Classical Ras Effector`, color = 'black', size = 2, max.overlaps = 40)
graphics.off()

#
theme_none <- theme(axis.line = element_line(colour = "black"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.border = element_blank(),
                    panel.background = element_blank()) 

# generate the histograms across LFQ data
experiments <- unique(d2_background$experiment)
row <- levels(d2$row)

outdir <- 'derived/KRAS/210927_deep_proteome_apex/'

for (e in experiments){
  for (r in row){
    test <- d2[d2$experiment == e & d2$row == r,]
    test$lfq_hist <- test$LFQ > 5.5
    p1 <- ggplot(test[test$lfq_hist,], aes(x = LFQ)) +
      geom_histogram(binwidth=1.5, color = 'black') +
      theme_minimal() + ylim(0, 600) + theme_none
    outname <- file.path(outdir, paste0(e,'_',r,'_lfq_hist_right.pdf'))
    ggsave(outname, p1, width = 1.25, height = 2.5)
    p2 <- ggplot(test[!test$lfq_hist,], aes(x = LFQ)) +
      geom_histogram(binwidth=1.5, color = 'black') +
      theme_minimal() + ylim(0, 600) + theme_none
    outname <- file.path(outdir, paste0(e,'_',r,'_lfq_hist_left.pdf'))
    ggsave(outname, p2, width = 1.25, height = 2.5)
    
  }
}

# 


# plot data
ggplot(test, aes(x = LFQ,y = LogFC, color = type, label = Prey)) +
  geom_point(alpha = 0.9, color = 'grey') +
  geom_hline(yintercept=4, linetype = 'dashed', alpha = 0.5) +
  scale_x_continuous(breaks = 6:max(ceiling(d2$LFQ))) +
  col_scale +
  theme_bw() +
  xlab('Deep Proteome LFQ') +
  ggtitle('KRAS Scatter Plot of SaintScore versus deep proteome LFQ', subtitle) +
  geom_point(data = d2_ht, size = 2) +
  geom_point(data = d2_lt,  size = 2) +
  geom_point(data = d2_effector,  size = 2) 







