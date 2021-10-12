
library(data.table)

data <- fread('inst/extdata/allCorumComplexes.txt')
data$`subunits(Gene name)`

all_corum_complexes <- (do.call(rbind,lapply(1:nrow(data), function(i){
  entry = data[i,]
  complex = entry$ComplexName
  organism = entry$Organism
  genes = unique(unlist(strsplit(entry$`subunits(Gene name)`, '\\ *;\\ *')))
  d <- data.frame(genes = genes, complex = complex, organism = organism)
  return(d)
})))

all_corum_complexes <- all_corum_complexes[!duplicated(all_corum_complexes),]
save(all_corum_complexes, file = 'data/all_corum_complexes.rda', compress = 'xz')






