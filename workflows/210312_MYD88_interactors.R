# created by fhl
devtools::load_all('~/Projects/07_genoppi/Genoppi/')

## load in the data
d <- read.delim('~/Projects/06_kessler_rotation/data/SAINT MyD88 final analysisb.txt')
d <- d[d$BFDR < 0.1 & d$SaintScore >= 0.8, ]
head(d)

# InWeb
inweb_interactors <- get_inweb_list('MYD88')
inweb_overlap <- list(inweb = inweb_interactors$gene[inweb_interactors$significant],experiment = d$Prey)
venn <- draw_genoppi_venn(inweb_overlap, colors = c('yellow','red'), main = 'InWeb Overlap')
intersect(inweb_interactors$gene[inweb_interactors$significant], d$Prey)
plot_venn(venn)

# IrefIndex
irefindex_interactors <- get_irefindex_list('MYD88')
iref_overlap <- list(irefindex = irefindex_interactors$gene[irefindex_interactors$significant],experiment = d$Prey)
venn <- draw_genoppi_venn(iref_overlap, colors = c('blue','red'), main = 'irefindex Overlap')
intersect(irefindex_interactors$gene[irefindex_interactors$significant], d$Prey)
plot_venn(venn)


