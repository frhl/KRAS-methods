devtools::load_all()
library(data.table)

# setup file paths
in_dir <- 'inst/extdata/KRAS/210927/'
files <- list.files(in_dir, full.names = TRUE)
newnames <- gsub(' ','_',tools::file_path_sans_ext(basename(files)))

# read files
dts <- lapply(files, fread)
names(dts) <- newnames

# select proteins to use
dts <- lapply(names(dts), function(d_name){
  
  # annotate data
  d <- dts[[d_name]]
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

  return(d)
  
})

# get all proteins/genes that pass Andreas' threshold
# and subsequently subset the data for these genes
dt <- do.call(rbind, dts)
genes <- unique(dt$gene[abs(dt$logFC) >= 0.5 & dt$logp >= 0.8])
dt <- dt[dt$gene %in% genes]

# let's group the genes using a categorical metric
keep <- c(
  'activator','Classical Ras Effector',
  'Proximitome','receptor','repressor'
)

# merge with subsetted data.frame
mapping <- read_xlsx('inst/extdata/KRAS/20210826 KRAS interaction network-final.xlsx')
colnames(mapping) <- c('uniprot','gene', 'class1', 'class2', 'class3')
mapping <- mapping[mapping$class1 %in% keep,]
dt <- merge(dt, mapping, all.x = T)
dt$class1[is.na(dt$class1)] <- 'Novel'
dt$class1 <- as.factor(dt$class1)

# setup current condition for comparison
#condition <- 'Starvation'
conditions <- c('FCS','Starvation')

# For each group, we are doing horizontal clustering
out_table <- list()
max_obs_order <- c(0)
limits <- c(min(dt$logFC), max(dt$logFC))

# for each condition
query_classes <- levels(dt$class1)
clustered_dts <- lapply(query_classes, function(class){
  
  # Get genes to keep across both Starvartion and FCS
  bool_class <- dt$class1 %in% class
  dt_keep <- dt[bool_class,]
  n_max <- min(50, nrow(dt_keep))
  genes_keep <- dt_keep$gene[rev(order(abs(dt_keep$logFC)))[1:n_max]]
  
  # Get starvation and FCS condition
  clust_dt <- dt[bool_class,]
  clust_dt_stv <- clust_dt[clust_dt$condition == 'Starvation']
  clust_dt_fcs <- clust_dt[clust_dt$condition == 'FCS']
    
  # convert to "wide" format so that we can perform hclustering
  clust_dt_stv <- clust_dt_stv[,c('how','gene','logFC')]
  clust_dt_stv <- reshape(clust_dt_stv, idvar = "gene",timevar = "how", direction = "wide")
  clust_dt_fcs <- clust_dt_fcs[,c('how','gene','logFC')]
  clust_dt_fcs <- reshape(clust_dt_fcs, idvar = "gene",timevar = "how", direction = "wide")
  
  # save to list
  clust_dt_stv_out <- clust_dt_stv
  clust_dt_stv_out$in_plot <- clust_dt_stv$gene %in% genes_keep
  out_table[[paste0('Stv',class)]] <<- clust_dt_stv_out
  clust_dt_fcs_out <- clust_dt_fcs
  clust_dt_fcs_out$in_plot <- clust_dt_fcs$gene %in% genes_keep
  out_table[[paste('FCS',class)]] <<- clust_dt_fcs_out
  
  # subset by genes
  clust_dt_stv <-  clust_dt_stv[ clust_dt_stv$gene %in% genes_keep, ]
  clust_dt_fcs <-  clust_dt_fcs[ clust_dt_fcs$gene %in% genes_keep, ]
  enough_genes <- (nrow(clust_dt_stv) >= 2) | (nrow(clust_dt_fcs) >= 2)

  if (enough_genes){
    
    # remove missing values
    clust_dt_narm <- clust_dt_stv
    clust_dt_narm[is.na(clust_dt_narm)] <- 0
    clust_dt_narm <- as.matrix(clust_dt_narm[,2:4])
    
    # clustering
    dist_matrix <- dist(clust_dt_narm)
    clustering <- hclust(dist_matrix, method = 'complete')
    
    # for every class we increment order numbers
    order_matrix <- data.frame(gene = clust_dt_stv$gene, order = clustering$order)
    max_obs_order <<- c(max_obs_order, max(order_matrix$order))
    
    # combine with order_matrix
    melt_mrg_stv <- merge(melt(clust_dt_stv), order_matrix)
    melt_mrg_stv$condition <- 'Starvation'
    melt_mrg_fcs <- merge(melt(clust_dt_fcs), order_matrix)
    melt_mrg_fcs$condition <- 'FCS'
    
    # Combined FCS and STV for plotting purposes
    combined <- rbind(melt_mrg_stv, melt_mrg_fcs)
    combined$gene <- factor(combined$gene, levels = order_matrix$gene[order_matrix$order])
    colnames(combined) <- c('gene','how','logFC','order','condition')
    combined$how <- gsub('logFC\\.','',combined$how)
    
    # get p-values from dt
    dt_pval <- dt[,c('gene','how','condition','logp')]
    combined <- merge(combined, dt_pval, by = c('gene','how','condition'), all.x = TRUE)
    
    # deal with NAs and labels
    combined$logFC[is.na(combined$logFC)] <- 0
    combined$label <- ifelse(abs(combined$logFC) > 0.5 & combined$logp >= 0.8, '*', '')
    combined$condition <- factor(combined$condition, c("Starvation", "FCS"))
    
    p <- ggplot(combined, aes(x = how, y = gene, fill = logFC, label = label)) +
      geom_tile() +
      geom_text() +
      theme_bw() +
      xlab('Experiment') + 
      ylab('Prey protein') +
      ggtitle(paste0(class)) +
      scale_fill_gradient2(low = 'navyblue', mid = 'white', high = 'firebrick1', midpoint = 0, limits = limits) +
      labs(fill = expression(Log[2]*FC[mt/wt])) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) + facet_grid(~condition, scale="free")
    
    outplt <- paste0('derived/KRAS/tables/210930_volcano_comparison_heatmap_plot_',class,'.pdf')
    #l2 <- log2(n_max)*2
    ggsave(p, filename = outplt, width = 4, height = 2+(n_max)/7)
  }
})

# Write xlsx table with plotting data
out_xlsx <- 'derived/KRAS/tables/210930_volcano_comparison_heatmap_data.xlsx'
write_xlsx(x = out_table, path = out_xlsx)

