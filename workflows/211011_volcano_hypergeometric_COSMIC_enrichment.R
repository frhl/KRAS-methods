# Explore what corum complexes are enriched in the data through 

# load these libraries
library(genoppi)
devtools::load_all()

# setup file paths and read in files
in_dir <- 'inst/extdata/KRAS/210927/'
files <- list.files(in_dir, full.names = TRUE)
newnames <- gsub(' ','_',tools::file_path_sans_ext(basename(files)))
newnames <- gsub('Starvation','STV',newnames)
newnames <- gsub('_Vs_','_',newnames)
indata <- lapply(files, fread)
names(indata) <- newnames

indata$G12D_WT_FCS

# Iterate over data 
data <- list()
directions <- c('P','N','B') # positive/negative/both
for (d_name in names(indata)){
  for (direction in directions){
    
    # annotate data
    d <- indata[[d_name]]
    d$comparison <- d_name
    ncols <- ncol(d)
    d$how <- paste0(unlist(strsplit(d$comparison, split = '_'))[c(1,3)], collapse = '_')
    d$condition <- unlist(strsplit(d$comparison, split = '_'))[4]
    
    # setup new variables with better naming
    d$gene <- d[['Gene names']]
    d$logp <- d[["-Log(P-value)"]]
    d$logFC <- d[["Difference"]]
    
    # remove columns we don't need
    d <- d[,-(1:ncols), with = F]
    
    # set significance thresholds
    if (direction == 'P'){
      d$significant <- d$logp >= 0.8 & d$logFC >= 0.5
    } else if (direction == 'N'){
      d$significant <- d$logp >= 0.8 & d$logFC <= -0.5
    } else if (direction == 'B') {
      d$significant <- d$logp >= 0.8 & abs(d$logFC) >= 0.5
    } else {
      stop('direction input is not valid!')
    }
    lstname <- paste0(d_name,'_',direction)
    data[[lstname]] <- d
  }
}




# kegg database is a subset


# load COSMIC
cosmic <- fread('/well/lindgren/flassen/ressources/genesets/genesets/data/cosmic/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv')
#cosmic <- fread('~/Projects/08_genesets/genesets/data/cosmic/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv')
hallmarks <- sort(table(cosmic$HALLMARK))
cell_types <- sort(table(cosmic$CELL_TYPE))

# gene-hallmarks
d_hallmarks <- cosmic[,c('GENE_NAME','HALLMARK')]
colnames(d_hallmarks) <- c('gene','geneset')
d_hallmarks <- d_hallmarks[!duplicated(d_hallmarks),]

# cells
d_cell_type <- cosmic[,c('GENE_NAME','CELL_TYPE')]
colnames(d_cell_type) <- c('gene','geneset')
d_cell_type <- d_cell_type[d_cell_type$geneset != '',]
d_cell_type <- d_cell_type[!duplicated(d_cell_type),]


# get databases
dblist <- list(
  hallmarks = d_hallmarks,
  celltypes = d_cell_type
)

#sum(grepl('kegg',reactome))

dblist <- lapply(dblist, function(x) {
  colnames(x) <- c('gene','geneset')
  return(x)
})

# load background genes
human <- fread('/well/lindgren/flassen//ressources/genesets/genesets/data/biomart/protein_coding_genes.tsv')
#human <- fread('~/Projects/08_genesets/genesets/data/biomart/protein_coding_genes.tsv')

# HYPERGEOMETRIC OVERLAP TESTING
outlist <- list()
for (dbname in names(dblist)){
  
  # get database
  print(dbname)
  signature_pathways <- dblist[[dbname]]
  colnames(signature_pathways) <- c('gene','geneset')
  
  enrichment <- lapply(names(data), function(cur_sheet){
    
    print(cur_sheet)
    
    # get current data
    cur_data <- data[[cur_sheet]] 
    cur_data <- as.data.frame(cur_data[,colnames(cur_data) %in% c('gene','significant'), with = FALSE])
    
    # if current peptide is both significant and insignificant (e.g. duplicate), set to significant
    dup_gene <- cur_data$gene[duplicated(cur_data$gene)]
    cur_data <- cur_data[!duplicated(cur_data),]
    cur_data$significant[cur_data$gene %in% dup_gene ] <- TRUE
    cur_data <- cur_data[!duplicated(cur_data),]
    
    # only do analysis if enriched proteins are present
    if (any(cur_data$significant)){
      
      pathways <- unique(signature_pathways$geneset)
      background_genes <- unique(human$hgnc_symbol[!human$hgnc_symbol %in% cur_data$gene])
      cur_data_global <- rbind(cur_data, data.frame(gene = background_genes, significant = FALSE))
      
      # iterate over pathways
      my_enrichment <- lapply(pathways , function(pathway_name){
        
        # conditional analysis
        cur_pathway <- signature_pathways
        cur_pathway$significant <- cur_pathway$geneset %in% pathway_name
        cur_pathway <- data.frame(gene = cur_pathway$gene, significant = cur_pathway$significant)
        hyper_cond <- suppressWarnings(calc_hyper(cur_data, cur_pathway, intersectDf = data.frame(intersectN = F)))
        hyper_cond$statistics$list_name <- pathway_name
        hyper_cond$statistics$analysis <- 'conditional'
        
        # global analysis
        hyper_global <- suppressWarnings(calc_hyper(cur_data_global, cur_pathway, intersectDf = data.frame(intersectN = F)))
        hyper_global$statistics$list_name <- pathway_name
        hyper_global$statistics$analysis <- 'global'
        
        # combine
        combined <- rbind(hyper_cond$statistics, hyper_global$statistics)
        
        return(combined)
      })
      
      enrichment <- do.call(rbind, my_enrichment)
      enrichment$FDR[enrichment$analysis == 'conditional'] <- p.adjust(enrichment$pvalue[enrichment$analysis == 'conditional'])
      enrichment$FDR[enrichment$analysis == 'global'] <- p.adjust(enrichment$pvalue[enrichment$analysis == 'global'])
      enrichment_sorted <- enrichment[(order(enrichment$pvalue)),]
      
      #return(enrichment_sorted)
      
      #write.csv(enrichment_sorted, )
      outfile = paste0('derived/KRAS/tables/211011/',cur_sheet,'_',dbname,'_enrichment.csv')
      write.csv(enrichment_sorted, outfile)
      
    }
  })
  
  
  
  # write to a table
  #names(enrichment) <- names(data)
  #enrichment <- null_omit(enrichment)
  #outfile = paste0('derived/KRAS/tables/211001/',dbname,'_enrichment.xlsx')
  #write_xlsx(enrichment, outfile)
  #outlist[[dbname]] <- enrichment
  
}






