#4-6: Beads Control Starvation 
#7-9: Beads Control FCS
#13-15: WT Starvation
#16-18: WT FCS
#22-24: G12D starvation
#25-27: G12D FCS
#31-33: G13D Starvation
#34-36: G13 FCS
#40-42: Q61H Starvation
#43-45: Q61H FCS

library(ggplot2)
library(ggrepel)

# read master sheets which contains experiment well IDs
master <- read.csv('inst/extdata/kras_wells_master.csv')
master <- master[!is.na(master$experiment),]
master$experiment <- gsub(' ', '_', master$experiment)

# read MS data
d <- fread('../data/proteinGroups Kras.txt')
str(d)
ncol(d)


ggplot_pca <- function(pca){
  
  # setup data
  pca 
  pc_var <- (pca$sdev^2)/sum((pca$sdev)^2)
  pca_rot <- pca$rotation[,c(1,2)]
  pca_d <- data.frame(pca_rot)
  pca_d$group <- rownames(pca_rot)
  pca_d$replicate <- rep(1:3, length.out = nrow(pca_d))
  
  # pca contributions
  xlab_new <- paste('PC1:', round(pc_var[1],4)*100,'%')
  ylab_new <- paste('PC2:', round(pc_var[2],4)*100,'%')
  
  # setup plot
  p1 <- ggplot(pca_d, aes(PC1, PC2, color = group, label = group)) + 
    geom_point(size = 3) +
    xlab(xlab_new) +
    ylab(ylab_new) +
    geom_text_repel() +
    theme_bw() 
  
  return(p1)
  
}


# generate various PCA plots
pdf('derived/210331_kras_pca_ms_data_qc.pdf', width = 10, height = 10)
  ggplot_pca(prcomp(grep_ms_data(d, what = 'Peptides', master))) + ggtitle('Peptide Count (PCA)')
  ggplot_pca(prcomp(grep_ms_data(d, what = 'Razor \\+ unique peptides', master))) + ggtitle('Razor + Unique Peptides (PCA)')
  ggplot_pca(prcomp(grep_ms_data(d, what = 'Intensity', master))) + ggtitle('Intensity (PCA)')
  ggplot_pca(prcomp(grep_ms_data(d, what = 'LFQ intensity', master))) + ggtitle('LFQ Intensity (PCA)')
  ggplot_pca(prcomp(grep_ms_data(d, what = 'MS/MS count', master))) + ggtitle('MS/MS count (PCA)')
  ggplot_pca(prcomp(grep_ms_data(d, what = 'Sequence coverage', master))) + ggtitle('Sequence coverage % (PCA)')
graphics.off()

pdf('derived/210331_kras_zeros_nonzeros_qc.pdf', width = 10, height = 10)

## intensity
d_int <- grep_ms_data(d, what = 'Intensity', master)
d_int_zeros <- apply(d_int, 2, function(x) sum(x == 0))
plot(d_int_zeros, main = 'Count of intensities equal to zero across samples')
text(d_int_zeros, names(d_int_zeros))
par(mar = c(15.1, 4.1, 4.1, 4.1))
d_int_not_zeros <- apply(sapply(d_int, as.integer), 2, function(x) log10(as.numeric(x[x!=0 & !is.na(x)])))
boxplot(d_int_not_zeros, las = 2, main = 'log10(Intensities) (that are not zero) across samples')

## intensity
d_int <- grep_ms_data(d, what = 'Peptides', master)
d_int_zeros <- apply(d_int, 2, function(x) sum(x == 0))
plot(d_int_zeros, main = 'Count of Peptides equal to zero across samples')
text(d_int_zeros, names(d_int_zeros))
par(mar = c(15.1, 4.1, 4.1, 4.1))
d_int_not_zeros <- apply(sapply(d_int, as.integer), 2, function(x) (as.numeric(x[x!=0 & !is.na(x)])))
boxplot(d_int_not_zeros, las = 2, main = 'Peptides (that are not zero) across samples')

graphics.off()






