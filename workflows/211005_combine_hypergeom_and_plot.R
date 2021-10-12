# Explore what corum complexes are enriched in the data through 

# load these libraries
library(genoppi)
devtools::load_all()
library(writexl)

# files and databases
files <- list.files('derived/KRAS/tables/211001/', full.names = TRUE, pattern = '.csv')
dbs <- c("msigdb_h" , "go_bp", "go_cc" , "go_mf", "kegg")

# create an excel file
list_of_dbs <- lapply(dbs, function(d){
  selected <- files[grepl(d, files)]
  newnames <- gsub('_enrichment','',tools::file_path_sans_ext(basename(selected)))
  list_of_data <- lapply(selected, function(f){
    fread(f)
  })
  names(list_of_data) <- newnames
  outfile = paste0('derived/KRAS/tables/211001_', d,'_hypergeom_analysis.xlsx')
  #write_xlsx(list_of_data, path = outfile)
  return(list_of_data)
})
names(list_of_dbs) <- dbs

# Create heatmap
list_of_dbs


