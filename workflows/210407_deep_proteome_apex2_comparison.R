
# a script for comparing kras intensity data with deep proteome

library(readxl)
library(ggplot2)
library(genoppi)
library(cowplot)
devtools::load_all()


#pdf('derived/plots/210407_KRAS_deep_proteome_versus_apex2_minrep3.pdf', width = 12, height = 10)
#for (minrep in 1:3){

minrep = 2

# helpers
clean_lfq <- function(x, func = median, log = TRUE){
  keep <- apply(x[,1:3], 1, function(x) sum(x!=0) >= minrep)
  if (log){
    median <- apply(log10(x[,1:3]), 1, function(x) median(x[is.finite(x)]))
  } else {
    median <- apply(x[,1:3], 1, function(x) median(x[is.finite(x)]))
  }
  df <- data.frame(gene = x$gene, LFQ.median = func)
  df <- df[keep,]
  return(df)
}

# read master sheets which contains experiment well IDs
master <- read.csv('inst/extdata/KRAS/kras_wells_master.csv')
master <- master[!is.na(master$experiment),]
master$experiment <- gsub(' ', '_', master$experiment)

# read MS data for KRAS intensity values
d <- fread('../data/proteinGroups Kras.txt')
kras <- grep_ms_data(d, what = 'LFQ intensity', master)
kras$gene <- unlist(lapply(strsplit(kras$`Gene names`, split = ';'), function(x) x[1]))

# combine each experiment in list
experiments <- unique(master$experiment)
d_kras <- lapply(experiments, function(x) { as.data.frame(kras[, grepl(x, colnames(kras)) | colnames(kras) %in% 'gene', with = F])})
names(d_kras) <- experiments

# kras
x <- d_kras[[1]]
d_x <- lapply(d_kras, clean_lfq)
d_x <- lapply(d_x, function(x){colnames(x) = c('gene','kras.LFQ'); return(x)})


# load deeep proteome
dp <- read_xlsx('inst/extdata/KRAS/11 Cell line proteome data mcp.M111.014050-3.xlsx', skip = 1)
dp$gene <- unlist(lapply(strsplit(dp$`Gene Names`, split = ';'), function(x) x[1]))
dp[dp == 'NaN'] <- NA

# subset deep proteoome by LFQ
dp_hek <- dp[,c('LFQ Intensity HEK293_1','LFQ Intensity HEK293_2','LFQ Intensity HEK293_3')]
dp_hek <- as.data.frame(apply(dp_hek, 2, as.numeric))
dp_hek[is.na(dp_hek)] <- 0
dp_hek$gene <- dp$gene

d_y <- clean_lfq(dp_hek, log = F)
colnames(d_y) <- c('gene','dp.LFQ')
write.csv(d_y, 'derived/KRAS/tables/210408_deep_proteome_min2reps.csv', quote = F, row.names = F)

## combine and merge deep proteome data with KRAS
combi <- do.call(rbind, lapply(names(d_x), function(x){
  mrg <- merge(d_x[[x]], d_y)
  mrg$name <- x
  return(mrg)
}))

# make combi
effectors <- fread('inst/extdata/KRAS/effector_genes.txt')$gene
combi$what <- "NaN"
combi$color <- 'black'
combi$what[combi$gene %in% effectors] <- 'Effector'
combi$color[combi$gene %in% effectors] <- 'blue'
combi$what[combi$gene %in% 'KRAS'] <- 'KRAS'
combi$color[combi$gene %in% 'KRAS'] <- 'red'

# subset combi
combi_s <- combi[combi$gene %in% effectors | combi$gene %in% 'KRAS',]

plt <- ggplot(combi, aes(x = kras.LFQ, y = dp.LFQ, label = gene)) + 
  geom_point(alpha = 0.2, color = 'lightblue') + 
  theme_bw() + 
  facet_wrap(~name)
plt <- plt + geom_point(data = combi_s, color = combi_s$color)
plt <- plt + ggrepel::geom_text_repel(data = combi_s)

# add text
plt <- plt + 
  ggtitle('Deep proteome experiment (LFQ) vs KRAS apex 2 experiment (LFQ)',
          paste('Minimum non-zero replicates:', minrep)) +
  xlab('KRAS apex-2 log10(LFQ)') +
  ylab('KRAS Deep proteome log10(LFQ)')

print(plt)
#ggsave('derived/plots/KRAS_deep_proteome_versus_apex2_minrep3.pdf', width = 12, height = 10)
#}
#graphics.off()
