# compare WT FCS versus WT starvation in scatter plot fashion

# load these libraries
library(readxl)
library(ggplot2)
library(genoppi)
library(ggrepel)
library(cowplot)
devtools::load_all()

# setup paths
path <- 'inst/extdata/KRAS/210216_SAINTexpress_KRAS.xlsx'
sheets <- excel_sheets(path)
in_sheets <- sheets[grepl('(WT)|(G|Q)',sheets)]
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
data <- set_saintscore_workbook(data)
saintscore <- 0.6

# effectors and interactors
effector_df <- read_excel('inst/extdata/KRAS/K-Ras known interactor and putative effectors  new.xlsx', 'Strong interactors')
effector_lst <- lapply(colnames(effector_df), function(x) unique(na.omit(effector_df[[x]])))
names(effector_lst) <- colnames(effector_df)
gene_effectors <- fread('inst/extdata/KRAS/effector_genes.txt')$gene
gene_effectors <- gene_effectors[! gene_effectors %in% 'LZTR1']
gene_effectors <- gene_effectors[! gene_effectors %in% '']
effector_lst$effectors <- gene_effectors
effectors <- data.frame(Prey = effector_lst$effectors, Type = 'effectors')
unique_interactors <- data.frame(Prey = effector_lst$`UNIQUE INTERACTORS  Kovalski JR (2019)`, Type = 'Unique Interactors')
hc_interactors <- fread('inst/extdata/KRAS/hc_ehc_kras_interactors.csv')
hc_interactors <- data.frame(Prey = hc_interactors$EHC[hc_interactors$EHC != "" ], Type ="High Confidence Interactors")


# throughput
int <- as.data.frame(read_xlsx('inst/extdata/KRAS/BIOGRID-GENE-110043-4.3.196.New.xlsx'))
int_ht <- int[int$`Experimental System` %in% 'Proximity Label-MS' & int$Throughput == 'High Throughput',] # not labelled
int_ht <- data.frame(Prey = unique(c(int_ht$`Official Symbol Interactor A`, int_ht$`Official Symbol Interactor B`)), Type = 'High Throughput')
int_lt <- int[int$Throughput == 'Low Throughput',] 
int_lt <- data.frame(Prey = unique(c(int_lt$`Official Symbol Interactor A`, int_lt$`Official Symbol Interactor B`)), Type = 'Low Throughput')

# get overlay df and remove duplicates
overlay_df <- rbind(effectors, int_lt, int_ht)
overlay_df <- overlay_df[! (overlay_df$Prey %in% intersect(int_ht$Prey, int_lt$Prey) & overlay_df$Type == 'High Throughput'),]
overlay_df <- overlay_df[! (overlay_df$Prey %in% intersect(int_lt$Prey, effectors$Prey) & overlay_df$Type == 'Low Throughput'),]

#overlay_df[overlay_df$Prey %in% 'LZTR1',]

#overlay_df <- rbind(effectors, hc_interactors, unique_interactors)


lst <- list()

for (x1 in names(data)){
  for (x2 in names(data)){
    
    # load data
    x1_data <- data[[x1]]
    x2_data <- data[[x2]]
    prey_keep <- unique(c(x1_data$Prey[x1_data$SaintScore >= saintscore],
                          x2_data$Prey[x2_data$SaintScore >= saintscore]))
    
    # merge data
    mrg <- merge(x1_data, x2_data, all.x = T, all.y = T, by.x = 'Prey',by.y = 'Prey', suffixes = c('.x1','.x2'))
    mrg$LogFC.x1[is.na(mrg$LogFC.x1)] <- 0
    mrg$LogFC.x2[is.na(mrg$LogFC.x2)] <- 0
    mrg <- mrg[mrg$Prey %in% prey_keep,]
    
    # merge with effectors
    plot_df <- merge(mrg, overlay_df, all.x = T)
    plot_df <- plot_df[!duplicated(plot_df$Prey),]
    
    # color
    plot_df$Type[is.na(plot_df$Type)] <- 'N/A'
    plot_df$Type <- as.factor(plot_df$Type)
    colors_levels <- c("effectors", "High Throughput", "Low Throughput", "N/A")
    my_colors <-c('blue','orange','red','black')
    names(my_colors) <- colors_levels
    color_scale <- scale_colour_manual(name = "Type",values = my_colors)
    
    # scatter plot
    plt <- ggplot(plot_df, aes(x = LogFC.x1, y = LogFC.x2, color = Type, label = Prey)) +
      geom_point(data = plot_df[plot_df$Type == 'N/A',], size = 1.5, alpha = 1, color = 'black', show.legend = F) + 
      geom_point(data = plot_df[plot_df$Type == 'High Throughput',], size = 1.5, show.legend = F) +
      geom_point(data = plot_df[plot_df$Type == 'Low Throughput',], size = 1.5, show.legend = F) + 
      geom_point(data = plot_df[plot_df$Type == 'effectors',], size = 1.5, show.legend = F) + 
      #geom_label_repel(data = plot_df[!is.na(plot_df$Type) & (plot_df$Type != 'High Throughput'),], show.legend  = F, 
      #                 max.overlaps = Inf, box.padding = 1, size = 3) +
      geom_label_repel(data = plot_df[plot_df$Prey %in% "LZTR1", ], show.legend  = F, 
                       max.overlaps = Inf, box.padding = 1, size = 4) +
      xlab(paste('Log2', x1)) +
      ylab(paste('Log2', x2)) +
      #ggtitle(paste(x1,'-',x2)) +
      color_scale + 
      geom_abline(linetype = 'dashed') +
      theme_classic()
    
    outname <- paste0(x1,'-',x2)
    lst[[outname]] <- plt
    
  }
}

length(lst)

pdf('derived/KRAS/plots/210528_all_scatterplots_nolabel.pdf', width = 24, height = 24)
plot_grid(plotlist = lst, nrow = 8, ncol = 8, labels = 1:64)
graphics.off()

# 
p1 <- lst$`G12_FCS-G12_Starvation`
p2 <- lst$`G12_FCS-WT_Starvation`

ggsave('derived/plots/210519_G12_FCS-G12_starvation.pdf', p1, width = 10, height = 8)
ggsave('derived/plots/210519_G12_FCS-WT_Starvation.pdf', p2, width = 10, height = 8)

# merge data



#ggsave('derived/plots/210423_wtstarv_vs_wtfcs_scatterplot_selected_ss04.pdf', plt, width = 10, height = 8)


