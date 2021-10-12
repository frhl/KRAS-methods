devtools::load_all()
library(igraph)
library(qgraph)
library(RColorBrewer)
library(genoppi) # can be downloaded from github

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# what proteins are shared between all of the 



