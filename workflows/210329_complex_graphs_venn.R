# libraries (all are strictly required to run the following code)
devtools::load_all()

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- filter_workbook(data, SaintScore >= 0.8 & LogFC > 4)

# protein complexes
data("all_corum_complexes")
corum_human <- all_corum_complexes[all_corum_complexes$organism == 'Human',]
db <- corum_human 

#pdf('derived/plots/210318_corum_graph_ss06_lf4.pdf', width = 14, height = 14)
params <- 'Corum Complexes (SaintScore >= 0.8 & LogFC > 4)'





#for (name in names(data)){
  
  # generate plots
  df <- data[[name]][,c(1,2)]
  complexes <- annotate_prey_complexes(df$Prey, db, n = 3)
  #make_complex_graph(df, complexes, main = paste(name, '\n', params))
  make_pathway_qgraph(head(df,10), NULL)
  
  # generate 
  outfile = paste0('derived/tables/210318_corum_table_',name,'_ss06_lf4.txt')
  complexes_table <- get_prey_complexes(df$Prey, db)
  write.table(complexes_table, outfile, row.names = F, quote = F, sep = '\t')
  
  
#}
#graphics.off()

