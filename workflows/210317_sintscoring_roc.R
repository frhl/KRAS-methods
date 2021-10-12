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

# databases
set.seed(1)
inweb <- get_inweb_list('KRAS')
irefindex <- get_irefindex_list('KRAS', n = 1)
ok_genes <- unique(c(inweb$gene[!inweb$significant], irefindex$gene[!irefindex$significant]))
random <- unique(sample(ok_genes, 1000))

hitproteins <- fread('inst/extdata/210317_kras_hitproteins.txt')
hitproteins <- hitproteins$Name


# plot saintscore as a function of
# proteins found / proteins found among hit proteins 

result <- do.call(rbind, lapply(names(data), function(sheet){
  d <- data[[sheet]]
  stacked_data <- do.call(rbind, lapply(seq(0,1,by=0.05), function(i){
    subset <- d$Prey[d$SaintScore >= i]
    count <- length(subset)
    known_count_hitproteins <- sum(subset %in% hitproteins)
    known_count_inweb <- sum(subset %in% inweb$gene[inweb$significant])
    known_count_irefindex <- sum(subset %in% irefindex$gene[irefindex$significant])
    
    tpr_hitproteins <-  known_count_hitproteins/count 
    tpr_inweb <-  known_count_inweb/count 
    tpr_irefindex <-  known_count_irefindex/count 
      
    df <- data.frame(threshold = i, count = count, tpr_hitproteins = tpr_hitproteins, tpr_inweb = tpr_inweb, tpr_irefindex = tpr_irefindex)
    return(df)
  }))
  
  stacked_data$experiment <- sheet
  return(stacked_data)
  
}))


plt <- ggplot(result, aes(y = threshold, x = tpr_hitproteins, color = experiment)) +
  geom_point() + 
  ylab('SaintScore Threshold') + 
  xlab('TPR') +
  facet_wrap(~experiment)
plt

ggsave('210317_TPR_vs_SaintScore.pdf', plt, width = 10, height = 10)

