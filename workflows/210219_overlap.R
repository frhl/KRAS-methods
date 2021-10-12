# load these libraries
library(readxl)
library(ggplot2)
library(genoppi)
library(ggrepel)
library(RColorBrewer)
devtools::load_all()

# data path
path <- 'inst/extdata/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G12)|(G13)|(Q61H)',sheets)]
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
#in_sheets <- in_sheets[grepl('(starv)', tolower(in_sheets))]


# helper functions
to_df <- function(df){
  df1 <- df[,c('Prey','SaintScore')]
  df1$significant = df1$SaintScore >= 0.8
  colnames(df1) <- c('gene','saintscore','significant')
  return(as.data.frame(df1))
}

# function for finding duplicated names
lower_triangle_duplicated <- function(x, split = '\\.'){
  return(duplicated(unlist(lapply(strsplit(x, split = split), function(x) paste(sort(x), collapse = '.')))))
}

# for quickly indexing
index <- function(name, i){unlist(lapply(strsplit(name, split = '_'), function(x) x[i]))}

# how similar are the IPs really?
dnames <- names(data)
stats <- list()
for (dname1 in dnames){
  stats[[dname1]] <- list()
  for (dname2 in dnames){
    
    # load data
    d1 <- to_df(data[[dname1]])
    d2 <- to_df(data[[dname2]])
    
    # calculate hypergeometric overlap
    statistics <- calc_hyper(d1, d2, intersectDf = data.frame(intersectN = T), bait = 'KRAS')
    
    # get overlap between interacors
    overlap_significant <- intersect(d1$gene[d1$significant], 
                                      d2$gene[d2$significant])
    d1_only <- d1$gene[d1$significant & (! d1$gene %in% overlap_significant)]
    d2_only <- d2$gene[d2$significant & (! d2$gene %in% overlap_significant)]
    overlap_pct <- length(overlap_significant) / 
      (length(d1_only) + length(d2_only) + length(overlap_significant))
    
    # add to statistics df
    statistics$statistics$pvalue_lt_005 <- statistics$statistics$pvalue < 0.05
    statistics$statistics$d1_n_only <- length(d1_only)
    statistics$statistics$d2_n_only <- length(d2_only)
    statistics$statistics$overlap_n<- length(overlap_significant)
    statistics$statistics$overlap_pct <- overlap_pct
    
    # combine data
    namedf <- data.frame(
      bait = 'KRAS',
      cell1 = index(dname1, 1),
      condition1 = index(dname1, 2),
      cell2 = index(dname2, 1),
      condition2 = index(dname2, 2),
      SaintScore = 0.8,
      logfc = 0
    )
    
    stats[[dname1]][[dname2]] <- cbind(namedf, statistics$statistics)
  }
}


# load external dataests
path <- 'inst/extdata/K-Ras known interactor and putative effectors  new.xlsx'
effector_df <- read_excel(path, 'Strong interactors')
effector_lst <- lapply(colnames(effector_df), function(x) data.frame(gene = unique(na.omit(effector_df[[x]])), significant = TRUE))
names(effector_lst) <- colnames(effector_df)
inweb <- get_inweb_list('KRAS')
effector_lst$InWeb <- inweb[inweb$significant,]
effectors <- fread('inst/extdata/effector_genes.txt')$gene
effectors <- effectors[! effectors %in% 'LZTR1']
effectors <- effectors[! effectors %in% '']
effector_lst$effectors <- data.frame(gene = effectors, significant = TRUE)


# calculate enrichment of various PPI genesets
databases <- list()
for (dname1 in dnames){
  databases[[dname1]] <- list()
  for (dbname in names(effector_lst)){
    
    # combine data
    namedf <- data.frame(
      bait = 'KRAS',
      cell = index(dname1, 1),
      condition = index(dname1, 2),
      database = dbname,
      SaintScore = 0.8,
      logfc = 0
    )
    
    # load data
    d1 <- to_df(data[[dname1]])
    db <- effector_lst[[dbname]]
    
    # calculate statistics
    statistics <- calc_hyper(d1, db, intersectDf = data.frame(intersectN = F), bait = 'KRAS')
    statistics$statistics$pvalue_lt_005 <- statistics$statistics$pvalue < 0.05
    statistics$statistics$overlap_genes <- paste(statistics$genes$mylist$successInSample_genes, collapse = '; ')
    databases[[dname1]][[dbname]] <- cbind(namedf, statistics$statistics)
  }
}


# combine lists and write excel file
stats_df <- do.call(rbind, lapply(stats, function(x) do.call(rbind, x)))
stats_df <- stats_df[!lower_triangle_duplicated(rownames(stats_df)),]
stats_df <- stats_df[order(stats_df$cell1, stats_df$condition1, stats_df$cell2, stats_df$condition2),]
database_df <- do.call(rbind, lapply(databases, function(x) do.call(rbind, x)))
database_df <- database_df[order(database_df$bait, database_df$cell, database_df$pvalue),]

outlist <- list(overlap = stats_df, database_overlap = database_df)
write_xlsx(outlist, path = 'derived/tables/210419_overlap_statistics.xlsx')

