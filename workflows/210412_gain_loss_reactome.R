# what genesets are enriched in the WT versus MT.

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
in_sheets <- in_sheets[grepl('(starv)', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
#data <- set_saintscore_workbook(data)
#data <- filter_workbook(data, SaintScore > 0.8)
wt <- data[["WT_Starvation"]]
data[["WT_Starvation"]] <- NULL

# load up/down regulaton data
reg <- as.data.frame(fread('derived/tables/210406_kras_log_wt_mt.csv'))
colnames(reg)[2:4] <- c("G12_Starvation", "G13_Starvation", "Q61H_Starvation")
reg$cluster_complete_order <- NULL
reg$Prey <- reg$wt_prey
reg$wt_prey <- NULL

# setup pathway enrichment
#signature_pathways <- all_corum_complexes
#signature_pathways <- signature_pathways[signature_pathways$organism == 'Human',]
#signature_pathways$organism <- NULL
#signature_pathways <- msigdb_h_table
#signature_pathways <- goa_cc_table
signature_pathways <- all_reactome_complexes[,c(1,4)]
signature_pathways <- signature_pathways[duplicated(signature_pathways),]
#signature_pathways$GO.ID <- NULL
colnames(signature_pathways) <- c('gene','geneset')
bait = 'KRAS'
pathway_enrichment <- calc_nested_enrichment(data, signature_pathways)
names(pathway_enrichment) <- names(data)

#writexl::write_xlsx(pathway_enrichment, 'derived/210412_regulation_pathway_enrichment.xlsx')
extracted_pathways <- lapply(pathway_enrichment, function(x){ head(x$list_name, 30)})
names(extracted_pathways) <- names(pathway_enrichment)

# generate complex data
pdf('derived/plots/210413_regulation_corum_enrichment.pdf', width = 12, height = 14)
for (name in names(pathway_enrichment)){
  
  # subset data we want to plot
  cur_data <- data[[name]]
  cur_data <- as.data.frame(cur_data[cur_data$SaintScore >= 0.8,])
  
  # get complexes data
  cur_pathway <- extracted_pathways[[name]]
  cur_data_pathways <- as.data.frame(signature_pathways[signature_pathways$geneset%in% cur_pathway,])
  cur_pathways_list <- lapply(unique(cur_data_pathways$geneset), function(x) cur_data_pathways$gene[cur_data_pathways$geneset %in% x])
  names(cur_pathways_list) <- unique(cur_data_pathways$geneset)
  
  # get regulation (colors)
  cur_reg <- reg[,colnames(reg) %in% c(name,"Prey")]
  colnames(cur_reg) <- c('val', 'Prey')
  node_colors <- assign_color_scale(cur_reg, colors = c('yellow','white','red'))
  
  # convert to list 
  #plot(node_colors$val, col = node_colors$color)
  prey_unique <- unique(node_colors$Prey)
  cur_col_list <- lapply( prey_unique , function(x) node_colors$color[node_colors$Prey %in% x] )
  names(cur_col_list) <- prey_unique
  
  title = paste(name, '\n Reactome Curated Complexes - Enrichment Threshold: SaintScore > 0.8 ONLY')
  
  area_force = 1.95 #2
  repulse = 3.05 #2.9
  
  # integrate data
  setup_group_graph(cur_data, cur_pathways_list, simple_plot = T,
                    node_colors = cur_col_list,
                    size_group = 8, size_not_group =  4, 
                    area_force =  area_force, repulse_rad = repulse, 
                    plot_labels = "NONE", main = title)
  
  setup_group_graph(cur_data, cur_pathways_list, simple_plot = T,
                    node_colors = cur_col_list,
                    size_group = 8, size_not_group =  4, 
                    area_force =  area_force, repulse_rad = repulse, 
                    plot_labels = 'NODES', main = title)
  
  # get P-value for enrichment
  selected_pathways <- pathway_enrichment[[name]]
  selected_pathways_stats <- selected_pathways[selected_pathways$list_name %in% cur_pathway,]
  cur_pathways_list_sig <- cur_pathways_list
  
  # count complexes
  count <- lapply(cur_pathways_list_sig,function(x) paste0(sum(x %in% cur_data$Prey),'/',length(x)))
  names(cur_pathways_list_sig) <- paste(names(cur_pathways_list_sig), count)
  
  setup_group_graph(cur_data, cur_pathways_list_sig, simple_plot = T,
                    node_colors = cur_col_list,
                    size_group = 8, size_not_group =  4, 
                    area_force =  area_force, repulse_rad = repulse, 
                    plot_labels = T, main = title)
  
}
graphics.off()










