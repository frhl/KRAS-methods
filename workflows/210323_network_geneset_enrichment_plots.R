# libraries (all are strictly required to run the following code)
devtools::load_all()
library(RColorBrewer)
library(genoppi) # can be downloaded from github. Ask frederik.

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get MS/MS data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data, SaintScore >= 0.6 & LogFC > 4)

# get go enrichment genes
files <- list.files(path = 'derived/tables/', pattern = '210319', full.names = T)
files <- files[grepl('\\.xlsx',files)]

# get 

name = names(data)[1]

for (name in names(data)){
  
  # read go enrichment
  go <- read_excel(files[1], sheet = name)[1:10,]
  geneset <- strsplit(go$overlap_genes, split = ';')
  names(geneset) <- go$dataset
  
  # read ppi data
  df <- data[[name]][,c(1,2)]
  #make_pathway_graph(df, geneset, remove_unassigned = T, collapse_overlap = '\n')
  the_graph(df, geneset, collapse_overlap = '\n')
  
}
graphics.off()

