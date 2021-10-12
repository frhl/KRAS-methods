# libraries (all are strictly required to run the following code)
devtools::load_all()
library(igraph)
library(qgraph)
library(RColorBrewer)
library(genoppi) # can be downloaded from github. Ask frederik.


# analysis variables
threshold <- "SaintScore >= 0.8, LogFC > 4"
date = '210321'

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]

# get data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
data <- lapply(data, function(x){
  x$significant <- x$SaintScore >= 0.8 & x$LogFC > 4
  return(x)
})


make_compatible <- function(x){
  colnames(x)[2] <- 'genes'
  x <- x[,c('genes','significant')]
  x <- as.data.frame(x)
  return(x)
}

# check how many protein protein interactions can
# be recapitulated across experiments
experiments <- names(data)
result <- do.call(rbind, lapply(experiments, function(e1){
  d1 <- make_compatible(data[[e1]])
  do.call(rbind, lapply(experiments, function(e2){
    d2 <- make_compatible(data[[e2]])
    
    # get overlap
    genes_origin <- unique(d1$genes[d1$significant])
    genes_target <- unique(d2$genes[d2$significant])
    
    # set operations
    genes_union <- unique(c(genes_origin, genes_target))
    genes_intersect <- intersect(genes_target, genes_origin)
    genes_not_intersect <- genes_union[! genes_union %in% genes_intersect]
    overlap_pct <- length(genes_intersect) / length(genes_union)
    
    #if (overlap_pct < 0.70) browser()
    # calculate hypergeoemtric enrichment
    hypergeom <- calc_hyper(d1, d2, bait = 'KRAS', intersectDf = data.frame(intersectN = F))
    
    # format output
    outdf <- data.frame(experiment_origin = e1, experiment_target = e2, threshold = threshold)
    outdf <- cbind(outdf, hypergeom$statistics)
    outdf$overlap_pct <- round(overlap_pct,3)*100
    outdf$intersect <- paste(sort(unique(genes_intersect)), collapse = ';')
    outdf$non_intersect <- paste(sort(unique(genes_not_intersect)), collapse = ';')
    outdf$label <- paste0(outdf$overlap_pct, '%\n(', length(genes_intersect),'/', length(genes_union),')')
      
    return(outdf)
  }))
}))

#write.table(result, 'derived/tables/210322_experiments_hypergeometric_overlap.txt', sep = '\t', quote = F, row.names = F)



# make heatmap
mat_df <- reshape(result[,c(1,2,14)], idvar = 'experiment_origin', timevar = 'experiment_target',  direction = 'wide')
row_names <- mat_df$experiment_origin
mat_df$experiment_origin <- NULL
col_names <- gsub('label.','',colnames(mat_df))
dim(mat_df)
str(mat_df)
mat <- as.matrix(mat_df)
rownames(mat) <- row_names
colnames(mat) <- col_names

# generate numerical matrix for clustering
nummat <- apply(as.matrix(gsub('%\n.+$','',mat)),2, as.numeric)
rownames(nummat) <- row_names
colnames(nummat) <- col_names

# generate matrix for labels
labelmat <- mat
#labels <- round(as.matrix(mat), 2)

# color palettte
palette = colorRampPalette(c('white', 'orange'))(n = 200)
outfile = paste0('derived/plots/',date,'_experiments_hypergeometric_overlap_heatmap_logfc_4.pdf')
pdf(outfile, width = 12, height = 10)
gplots::heatmap.2(nummat,
                  main = threshold,
                  cellnote = labelmat,    # same data set for cell labels
                  notecol="black",      # change font color of cell labels to black
                  density.info="none",  # turns off density plot inside color legend
                  trace="none",         # turns off trace lines inside the heat map
                  margins =c(12,9),     # widens margins around plot
                  col=palette)       # use on color palette defined earlier

graphics.off()















