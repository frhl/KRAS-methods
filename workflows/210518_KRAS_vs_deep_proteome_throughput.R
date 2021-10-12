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

# effectors and interactors
effector_df <- read_excel('inst/extdata/KRAS/K-Ras known interactor and putative effectors  new.xlsx', 'Strong interactors')
effector_lst <- lapply(colnames(effector_df), function(x) unique(na.omit(effector_df[[x]])))
names(effector_lst) <- colnames(effector_df)
gene_effectors <- fread('inst/extdata/KRAS/effector_genes.txt')$gene
gene_effectors <- gene_effectors[! gene_effectors %in% 'LZTR1']
gene_effectors <- gene_effectors[! gene_effectors %in% '']
effector_lst$effectors <- gene_effectors
effectors <- effector_lst$effectors
unique_interactors <- effector_lst$`UNIQUE INTERACTORS  Kovalski JR (2019)`
hc_interactors <- fread('inst/extdata/KRAS/hc_ehc_kras_interactors.csv')
hc_interactors <- hc_interactors$EHC[hc_interactors$EHC != "" & !hc_interactors %in% unique_interactors]

# throughput
int <- as.data.frame(read_xlsx('inst/extdata/KRAS/BIOGRID-GENE-110043-4.3.196.New.xlsx'))
int_ht <- int[int$`Experimental System` %in% 'Proximity Label-MS' & int$Throughput == 'High Throughput',] # not labelled
int_ht <- unique(c(int_ht$`Official Symbol Interactor A`, int_ht$`Official Symbol Interactor B`))
int_lt <- int[int$Throughput == 'Low Throughput',] 
int_lt <- unique(c(int_lt$`Official Symbol Interactor A`, int_lt$`Official Symbol Interactor B`))


# setup types
d$type <- 'N/A'
d$type[d$Prey %in% int_ht] <- 'High Throughput'
d$type[d$Prey %in% int_lt] <- 'Low Throughput'
d$type[d$Prey %in% effectors] <- 'Effector Protein' 
d$type <- as.factor(d$type)

# combine with deep proteome intensity
deep_proteome <- fread('derived/KRAS/tables/210408_deep_proteome_min2reps.csv')
colnames(deep_proteome) <- c("Prey",'LFQ')
d <- merge(d, deep_proteome, all.x = T)

# setup names
d$row <- as.factor(ifelse(grepl('Starv',d$Sheet),'Starvation','FCS'))
d$experiment <- unlist(lapply(strsplit(d$Sheet, split = '_'), function(x) x[1]))
d$experiment <- factor(d$experiment, levels = c("WT","G12", "G13","Q61H"))

# colors
my_colors <- c('blue','orange','red','black','orange') #brewer.pal(3,"Set1")
names(my_colors) <- levels(d$type)
col_scale <- scale_colour_manual(name = "type",values = my_colors)

## make jitter plots
d2 <- d[d$SaintScore >= 0.8,]
d2$LFQ[is.na(d2$LFQ)] <- 4 + runif(sum(is.na(d2$LFQ)))
d2_background <- d2[d2$type %in% 'N/A',]

# some annotations
subtitle = 'Deep proteome median LFQ (At least 2 replicates with non-zero LFQ values)'

# plot data
p2 <- ggplot(d2_background, aes(x = LFQ,y = LogFC, color = type, label = Prey)) +
  geom_point(alpha = 0.9, color = 'grey') +
  geom_hline(yintercept=4, linetype = 'dashed', alpha = 0.5) +
  scale_x_continuous(breaks = 6:max(ceiling(d2$LFQ))) +
  col_scale +
  theme_bw() +
  xlab('Deep Proteome LFQ') +
  ggtitle('KRAS Scatter Plot of SaintScore versus deep proteome LFQ', subtitle) + 
  facet_grid(row ~ experiment)
#facet_wrap(~Sheet) 

# for each analysis calculate p-value
d2_ht <- d2[d2$type %in% 'High Throughput',]
d2_lt <- d2[d2$type %in% 'Low Throughput',]
d2_effector <- d2[d2$type %in% 'Effector Protein',]


# overlay points
p2 <- p2 + 
  geom_point(data = d2_ht, size = 2) +
  geom_point(data = d2_lt,  size = 2) +
  geom_point(data = d2_effector,  size = 2) 

# plot
pdf('derived/KRAS/plots/210523_saintscore_vs_deep_proteome_lfq.pdf', width = 12, height = 7)
p2
p2 + geom_text_repel(data = d2_lt, color = 'black', size = 2, max.overlaps = 40)
p2 + geom_text_repel(data = d2_effector, color = 'black', size = 2, max.overlaps = 40)
graphics.off()


