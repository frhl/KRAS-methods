library(genoppi) # my own package for genetic analysis
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
data <- collapse_workbook(data)

# databases
set.seed(1)
inweb <- get_inweb_list('KRAS')
irefindex <- get_irefindex_list('KRAS', n = 1)
ok_genes <- unique(c(inweb$gene[!inweb$significant], irefindex$gene[!irefindex$significant]))
random <- unique(sample(ok_genes, 1000))

hitproteins <- fread('inst/extdata/210317_kras_hitproteins.txt')
hitproteins <- hitproteins$Name



# also manually load biogrid
#x <- fread('inst/extdata/BIOGRID-MV-Physical-4.3.195.mitab.txt')
#

# iterate over data and check overlap
saintscore_seq <- seq(0,1, by = 0.025)
result <- do.call(rbind, lapply(saintscore_seq, function(s){

  data_subset <- data[data$SaintScore >= s,]

  # hit database
  data_hitproteins <- data_subset[data_subset$Prey %in% hitproteins,]
  data_hitproteins <- aggregate(Prey ~ Sheet, data = data_hitproteins, length)
  data_hitproteins$PreyPct <- round( data_hitproteins$Prey / length(hitproteins), 3)*100
  data_hitproteins$database <- 'Hit Proteins (Andreas)'
    
  # inweb database
  data_inweb <- data_subset[data_subset$Prey %in% inweb$gene[inweb$significant],]
  data_inweb <- aggregate(Prey ~ Sheet, data = data_inweb, length)
  data_inweb$PreyPct <- round( data_inweb$Prey / sum(inweb$significant), 3)*100
  data_inweb$database <- 'InWeb'
  
  # irefindex database
  data_irefindex <- data_subset[data_subset$Prey %in% irefindex$gene[irefindex$significant],]
  data_irefindex <- aggregate(Prey ~ Sheet, data = data_irefindex, length)
  data_irefindex$PreyPct <- round( data_irefindex$Prey / sum(irefindex$significant), 3)*100
  data_irefindex$database <- 'IRefIndex'
  
  # random database
  data_random <- data_subset[data_subset$Prey %in% random,]
  data_random <- aggregate(Prey ~ Sheet, data = data_random, length)
  data_random$PreyPct <- round( data_random$Prey / length(random), 3)*100
  data_random$database <- 'Random'
  
  # biogrid database
  
  outdata <- rbind(data_hitproteins, data_inweb, data_irefindex, data_random)
  colnames(outdata) <- c('Experiment', 'Count','Pct','Database')
  outdata$Threshold <- s
  return(outdata)
}))

# plot percentage found versus saint score threshold
plt <- ggplot(result, aes(Pct, Threshold, color = Database)) +
  geom_point() + 
  ylab('SaintScore Threshold') + 
  xlab('Overlap with external database (%)') +
  ggtitle('SaintScores Versus Overlap with external PPI databases', 'Random: 1000 randomly sampled proteins not found in InWeb/IRefIndex') +
  facet_wrap(~Experiment)

plt
ggsave('derived/210317_SaintScore_vs_overlap.pdf', width = 14, height = 10)



## simulation

all_interactors <- inweb$gene[inweb$significant]
all_genes <- inweb$gene
n = 500

genes_sampled <- 500

# simulate random proteins
saintscore_seq <- seq(0,1, by = 0.05)
result_sim <- do.call(rbind, lapply(saintscore_seq, function(s){
  data_subset <- data[data$SaintScore >= s,]
  sim_random <- apply(replicate(n, sample(all_genes, genes_sampled) %in% data_subset$Prey ), 2, sum) / genes_sampled
  df_random <- data.frame(pct = sim_random, threshold = s, dataset = 'random')
  sim_inweb <- apply(replicate(n, sample(all_interactors, 100) %in% data_subset$Prey ), 2, sum) / 100
  df_inweb <- data.frame(pct = sim_inweb, threshold = s, dataset = 'inweb')
  outdata <- rbind(df_random, df_inweb)
  return(outdata)
}))

ggplot(result_sim, aes(pct, threshold, color = dataset)) +
  geom_point()







