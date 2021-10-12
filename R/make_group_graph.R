#' @title plot a group graph
#' @param df a data.frame containg two columns with 'Bait' and 'Prey'. If these are un-named
#' the function will assume that the first column is the 'Bait' and the second columns is 'Prey'.
#' @param geneset A named list. The list name correspond to each group and
#' each list item should be a string or vector of strings also present in df$prey.
#' Using a geneset will create a virtual node for each item of the list corresponding
#' to a group memebership.
#' @param remove_virtual_nodes boolean. Should virtual nodes from geneset be displayed?
#' @param convex_group_hull A named list. The list name correspond to each group and
#' each list item should be a string or vector of strings also present in df$prey. Each item 
#' in the list will be overlayed with a convex hull.
#' @param convex_hull_alpha numeric between 0 and 1. Opacity for the hull.
#' @param convex_hull_borders string. Borders for convex hull.
#' 
#' @param baits baits string or vector of strings. What are the baits in the data?
#' @param known_interactors vector of known gene-names of known interacors.
#' @param main title of the plot
#' @param node_colors list of preys. Names correspond to colors and entries correspond to nodes.
#' 
#' @param size_bait numeric. Size of bait.
#' @param size_group numeric. Size of group.
#' @param size_not_group numeric. Size of things not in groups.
#' @param simple_plot boolean. Removes bait and simplifies the plot.
#' 
#' @param area_force numeric. Force for generating the graph.
#' @param repulse_rad numeric. Repulse values. 
#' 
#' @param plot_labels string or boolean. Can be either TRUE, FALSE, "GROUPS", "NODES" or "NONE".
#' 
#' @family plotting
#' @export

make_group_graph <- function(df, geneset, convex_group_hull = NULL, baits = 'KRAS', 
                              edge_df = NULL, known_interactors = NULL, main = '',
                              node_colors = NULL, 
                              convex_hull_alpha = 0.5, convex_hull_borders = NA,
                              size_bait = 10, size_group = 10, size_not_group = 2, simple_plot = F,
                              area_force = 3.8, repulse_rad = 8.2, 
                              plot_labels = 'ALL', remove_virtual_nodes = F, verbose = T){
  
  
  # check inpuy geneset
  stopifnot(ncol(df) > 1)
  if (is.null(names(geneset))) {
    stop('param "geneset" must be a named list!')
    genesets_genes <- unlist(geneset) %in% unlist(df[,c(1,2)])
    if (verbose) write(paste('Note:',sum(genesets_genes), '/', length(genesets_genes), 'genes from param "geneset" was found in df. The rest will be discarded.'),stderr())
  }
  
  # check input convex hull
  if (!is.null(convex_group_hull)) {
    if (is.null(names(convex_group_hull))) stop('param "convex_group_hull" must be a named list!')
    if (!all(names(convex_group_hull) %in% names(geneset))) stop('names of named list "convex_group_hull" must all be names in "geneset"')
  }
  
  # check node colors
  if (is.null(node_colors)){
    if (is.null(names(node_colors))) stop('param "node_colors" must be a named list!')
    node_genes <- names(node_colors) %in% unlist(df[,c(1,2)])
    if (verbose) write(paste('Note:',sum(node_genes ), '/', length(node_genes ), 'genes from param "node_colors" was found in df. The rest will be discarded.'),stderr())
  }
  
  # check edge df
  if (!is.null(edge_df)){
    if (!inherits(edge_df, "data.frame")) stop('param "edge_df" must be a data.frame!')
    if (ncol(edge_df) < 2) stop('param "edge_df" must be a data.frame with at least two columns (source, target)!')
    if (nrow(edge_df) < 1) stop('param "edge_df" should at least contain one row!')
  }
  
  # check input df for unusual
  #if (length())
  #
  
  
  # setup nopes and verticies
  int_lst <- setup_int_group_df(df[,c(1,2)], geneset)
  
  # remove bait and non group items
  if (simple_plot){
    int_lst$edges <- int_lst$edges[!int_lst$edges$source %in% baits,]
    int_lst$nodes <- int_lst$nodes[int_lst$nodes$name %in% c(int_lst$edges$target, int_lst$virt_edges$source),]
  }
  
  # remove verticies not in geneset
  nodes = int_lst$nodes
  edges = int_lst$edges
  int_network <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  # set standard paramaters
  V(int_network)$bait <- V(int_network)$name %in% baits
  groups <- names(geneset)
  max_num_groups <- length(groups)
  
  # group colors
  colors <- brewer.pal(6, "Purples")
  group_node_colors <- as.list(rep(colors[2:length(colors)], length.out = length(groups)))
  names(group_node_colors) <- groups
  group_node_colors[['X']] <- 'cyan'
  
  # check input for user colors
  if (!is.null(node_colors)){
    stopifnot(!is.data.frame(node_colors))
    stopifnot(is.list(node_colors))
    na_color <- 'grey'
  }

  # replace node size / colors
  V(int_network)$size <- size_not_group
  V(int_network)$color <- group_node_colors[['X']]
  V(int_network)$shape <- 'circle'
  V(int_network)$label.cex <- 0.5
  
  for (x in V(int_network)$name) {
    
    # set genset nodes
    x_group <- int_lst$edges[int_lst$edges$source %in% names(geneset) & int_lst$edges$target %in% x,]
    
    if (nrow(x_group) > 0){
      V(int_network)$size[V(int_network)$name==x] <- size_group
      V(int_network)$color[V(int_network)$name==x] <- group_node_colors[[x_group$source[1]]]
    }
    
    # get bait groups
    if (x %in% baits) {
      V(int_network)$size[V(int_network)$name==x] <- size_bait
      V(int_network)$color[V(int_network)$name==x] <- "red"
    } 
    
    # if in geneset
    if (x %in% names(geneset)){
      V(int_network)$size[V(int_network)$name==x] <- 0
      #V(int_network)$shape[V(int_network)$name==x] <- 'rectangle'
      V(int_network)$color[V(int_network)$name==x] <- group_node_colors[[x]]
      V(int_network)$label.cex[V(int_network)$name==x] <- 0.75
    }
    
    # user defined colors
    if (!is.null(node_colors)){
      if (x %in% names(node_colors)){
        V(int_network)$color[V(int_network)$name==x] <- node_colors[[x]]
      } else {
        V(int_network)$color[V(int_network)$name==x] <- 'grey'
      }
    }
    
    
  }
  
  # label sizes
  V(int_network)$label <- NA
  V(int_network)$label[V(int_network)$size>=3 | V(int_network)$size==0] <- V(int_network)$name[V(int_network)$size>=3 | V(int_network)$size==0 ]
  
  # edges (colors and weights)
  E(int_network)$color <- 'grey'
  E(int_network)$color[!edges$source %in% names(geneset)] <- 'black'
  
  if (!is.null(known_interactors)){
    E(int_network)$color[edges$target %in% known_interactors & ! edges$source %in% names(geneset)] <- 'blue'
  }
  
  
  # a sub routine that sets up the width of the vectors
  width_vec <- 0.5
  if (!is.null(edge_df)){
    edgelist <- E(int_network)
    width_vec <- c()
    for (i in 1:length(edgelist)){
      source <-  V(int_network)$name[as.integer(tail_of(int_network, i))]
      target <-  V(int_network)$name[as.integer(head_of(int_network, i))]
      bool <- edge_df[[1]] == source & edge_df[[2]] == target
      width_vec <- c(width_vec, ifelse(any(bool), edge_df$width[bool], 0.5))
    }
  }
  
  # set up
  e <- get.edgelist(int_network,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(int_network),niter=7000,
                                         area=(vcount(int_network)^area_force),repulse.rad=(vcount(int_network)^repulse_rad))
  
  
  # remove virtual nodes
  if (remove_virtual_nodes){
    l <- l[-int_lst$virt_nodes$id,]
  }
  
  # setup plotting labels (groups or nodes)
  if (toupper(plot_labels) == 'GROUPS'){
    V(int_network)$label[! V(int_network)$label %in% groups] <- ''
  } else if (toupper(plot_labels) == 'NODES'){
    V(int_network)$label[  V(int_network)$label %in% groups] <- ''
  } else if (toupper(plot_labels) == 'NONE'){
    V(int_network)$label <- ''
  }
  
  
  # setup convex hulls
  convex_hull_list = NULL
  convex_hull_colors = NULL
  if (!is.null(convex_group_hull)){
    
    # get list of integeres corresponding to baits
    ids = data.frame(id = 1:length(V(int_network)), Prey = names(V(int_network)))
    ids = ids[!ids$Prey %in% groups,]
    stacked_geneset = stack(geneset)
    colnames(stacked_geneset) = c('Prey','group')
    stacked_geneset = stacked_geneset[stacked_geneset$group %in% names(convex_group_hull),]
    merged_groups = merge(ids, stacked_geneset)
    
    # setup colors and list of inndexes
    convex_hull_colors <- unlist(convex_group_hull)
    convex_hull_colors <- adjustcolor(convex_hull_colors, alpha.f = convex_hull_alpha)
    convex_group_names <- unique(names(convex_group_hull))
    convex_hull_list <- lapply(convex_group_names, function(x){merged_groups$id[merged_groups$group == x]})
    names(convex_hull_list) <- convex_group_names
    
  }

  # finally plot the data
  plot(int_network,
       vertex.color=V(int_network)$color,
       #V(int_network)$size+1, 
       vertex.frame.color="grey30",
       vertex.label= V(int_network)$label,
       vertex.label.color="black",
       vertex.label.cex=V(int_network)$label.cex,
       mark.groups = convex_hull_list, 
       mark.col = convex_hull_colors,
       mark.border = convex_hull_borders,
       edge.width = width_vec, #V(int_network)$width,
       layout=l, main = main)
  
  
}

























