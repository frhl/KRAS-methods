# Explore what corum complexes are enriched in the data through 
# networks that are colored by heatmap expression. Compares WT versus MT
# starvartion condition only.

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
in_sheets <- in_sheets[grepl('(starv)', tolower(in_sheets))]

# get proteomics data
data <- read_workbook(path, in_sheets, rename_saint_sheets(in_sheets))
#data <- set_saintscore_workbook(data)

# ext data
effector_df <- read_excel('inst/extdata/K-Ras known interactor and putative effectors  new.xlsx', 'Strong interactors')
effector_lst <- lapply(colnames(effector_df), function(x) unique(na.omit(effector_df[[x]])))
names(effector_lst) <- colnames(effector_df)
gene_effectors <- fread('inst/extdata/effector_genes.txt')$gene
gene_effectors <- gene_effectors[! gene_effectors %in% 'LZTR1']
gene_effectors <- gene_effectors[! gene_effectors %in% '']
effector_lst$effectors <- gene_effectors
effectors <- unique(unlist(effector_lst))

pdf('210419_kras_graph_effectors_based_on_logfc.pdf', width = 10, height = 10)
for (name in names(data)){
  
  # subset data we want to plot
  #name <- names(data)[1]
  cur_data <- as.data.frame(data[[name]])
  cur_data$val <- log2(cur_data$FoldChange)
  
  # remove non-effectors
  cur_data <- cur_data[cur_data$Prey %in% effectors,]
  cur_data <- cur_data[!cur_data$Prey %in% 'KRAS', ]

  
  # get widths for edges
  cur_data_ss <- cur_data
  cur_data_ss$width <- exp(cur_data_ss$SaintScore)
  edge_width_df <- cur_data_ss[,c('Bait','Prey','width', 'SaintScore')]
  
  # set color for edges
  cur_data$width <- exp(cur_data$SaintScore) / max(exp(cur_data$SaintScore))
  cur_data_edges <- cur_data
  cur_data_edges <- cur_data_edges$Prey[cur_data_edges$SaintScore >= 0.8]
  
  
  #edge_weight_lst <- lapply(cur_data_ss$Prey, function(x) cur_data_ss[cur_data_ss$Prey %in% x, c('Bait','Prey','weight')])
  
  #edge_color_lst <- lapply(cur_data_ss$Prey, function(x) data.frame(source='KRAS', target=x, color = 'blue'))
  
  #names(edge_color_lst) <- cur_data_ss$Prey

  # get regulation (colors) for nodes
  min_color <- -0.05
  max_color <- 10 #max(cur_data$val)*1.05
  colorize <- function(x) assign_color_scale(x, colors = c('white','red'), force_min =  min_color, force_max = max_color, force_length_out = 100)
  node_colors <- colorize(cur_data)
  node_colors <-  node_colors[!is.na(node_colors$color),]
  node_colors_lst <- as.list(node_colors$color)
  names(node_colors_lst) <- node_colors$Prey
  
  # get regulation (colors) for hulls
  convex_groups <- lapply(names(effector_lst), function(xname) {
    x <- effector_lst[[xname]]
    prey_in_data <- x[x %in% cur_data$Prey]
    group <- colorize(data.frame(Prey = xname, val = mean(node_colors$val[node_colors$Prey %in% prey_in_data])))
    return(group$color)
  })
  names(convex_groups) <- names(effector_lst)
  
  # titles
  area_force = 2.5 #2.3 #2
  repulse_rad =2.9 #2.9
  
  # setup legend
  legend <- suppressWarnings(setup_graph_legend(min_color, max_color, colorize, draw_every = 1, text = paste(name,'Log2FC')))
  print(legend)

  pdf('derived/plots/repulse_params.pdf', width = 10, height = 10)

  for (area_force in seq(1,5, by = 0.2)){
    
    for(repulse_rad in seq(1,5, by = 0.2)){
    
      NAME = paste('area_force =',area_force,', repulse_rad =',repulse_rad)
      print(NAME)
      # setup plot with node labels
      make_group_graph(cur_data, effector_lst, 
                       node_colors = node_colors_lst, 
                       convex_group_hull = convex_groups,
                       known_interactors =  cur_data_edges,
                       edge_df = edge_width_df,
                       remove_virtual_nodes = T,
                       size_group = 12, size_not_group =  4, 
                       area_force = area_force, repulse_rad = repulse_rad , 
                       plot_labels = 'NONE', main = NAME)  
      
    }
    
  }

  
  graphics.off()
  
  make_group_graph(cur_data, effector_lst, 
                   node_colors = node_colors_lst, 
                   convex_group_hull = convex_groups,
                   known_interactors =  cur_data_edges,
                   edge_df = edge_width_df,
                   remove_virtual_nodes = T,
                   size_group = 12, size_not_group =  4, 
                   area_force = area_force, repulse_rad = repulse_rad , 
                   plot_labels = 'NODES', main = name)
  
  make_group_graph(cur_data, effector_lst, 
                   node_colors = node_colors_lst, 
                   convex_group_hull = convex_groups,
                   known_interactors =  cur_data_edges,
                   edge_df = edge_width_df,
                   size_group = 12, size_not_group =  4, 
                   area_force = area_force, repulse_rad = repulse_rad , 
                   plot_labels = 'GROUPS', main = name)
  
  make_group_graph(cur_data, effector_lst, 
                   node_colors = node_colors_lst, 
                   convex_group_hull = convex_groups,
                   known_interactors =  cur_data_edges,
                   edge_df = edge_width_df,
                   size_group = 12, size_not_group =  4, 
                   area_force = area_force, repulse_rad = repulse_rad , 
                   plot_labels = T, main = name)
  
  
  
  
}
graphics.off()






