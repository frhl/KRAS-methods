
make_complex_graph <- function(df, complexes, baits = 'KRAS', geneset = NULL, main = NULL){
  
  # read in data
  colnames(df) <- c('int1','int2')
  df$interaction <- 'Bait-Prey'
  inweb_list <- get_inweb_list(baits)
  df$known_inweb_interactor <- df$int2 %in% inweb_list$gene[inweb_list$significant]
  
  # read in complexes
  df_complexes <- do.call(rbind, complexes)
  df_complexes <- df_complexes[df_complexes$int1 != df_complexes$int2,]
  df_complexes$interaction <- gsub('\\.[0-9]+$','',rownames(df_complexes))
  rownames(df_complexes) <- NULL
  df_complexes$known_inweb_interactor <- F
  
  # combinen data/complexes
  df_int <- unique(rbind(df, df_complexes))
  
  # define groups
  groups <- unique(df_int$interaction)
  groups <- groups[!groups %in% 'Bait-Prey']
  
  # setup nodes and edges
  tmp = unique(na.omit(c(as.character(df_int$int1), as.character(df_int$int2))))
  nodes = data.frame(name=tmp, id = 1:length(tmp))
  edges = data.frame(source=df_int$int1, target = df_int$int2, weight = 1)
  edges = edges[complete.cases(edges),]
  
  # rename nodes / edges
  nodes$name <- as.character(nodes$name)
  edges$source <- as.character(edges$source)
  edges$target <- as.character(edges$target)
  int_network <- graph_from_data_frame(d=edges, vertices=nodes, directed=FALSE)
  
  # set colors
  V(int_network)$bait <- V(int_network)$name %in% unique(df_int$int1)
  max_num_groups <- length(unique(df_int$interaction))
  node_colors <- brewer.pal(9,"Purples")
  
  # group colors
  group_node_colors <- as.list(rep(node_colors[2:6], 10, length.out = length(groups)))
  names(group_node_colors) <- groups
  group_node_colors[['Bait-Prey']] <- node_colors[1]
  
  # replace node size / colors
  V(int_network)$size <- NA
  V(int_network)$color <- NA
  for (x in V(int_network)$name) {
    
    group = unique(df_int[df_int$int1 %in% x | df_int$int2 %in% x, ]$interaction)
    
    if (length(group) == 1){
      V(int_network)$size[V(int_network)$name==x] <- 2
      V(int_network)$color[V(int_network)$name==x] <- group_node_colors[[group]]
    } else {
      V(int_network)$size[V(int_network)$name==x] <- 1+length(group) 
      V(int_network)$color[V(int_network)$name==x] <- group_node_colors[[group[2]]]  #node_colors[2]
    }
    
    if (x %in% baits) {
      V(int_network)$size[V(int_network)$name==x] <- 5
      V(int_network)$color[V(int_network)$name==x] <- "red"
    } 
    
  }
  
  # label sizes
  V(int_network)$label <- NA
  V(int_network)$label[V(int_network)$size>=3] <- V(int_network)$name[V(int_network)$size>=3]
  
  # set edge attributes (for edge color)
  E(int_network)$color <- 'grey'
  E(int_network)$color[df_int$known_inweb_interactor] <- 'blue'
  E(int_network)$color[!df_int$interaction %in% 'Bait-Prey'] <- 'black'
  
  # set up
  e <- get.edgelist(int_network,names=F)
  l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(int_network),niter=5000,
                                         area=(vcount(int_network)^1.7),repulse.rad=(vcount(int_network)^2.1))
  
  plot(int_network, edge.width=0.5,
       vertex.color=V(int_network)$color,
       vertex.size=V(int_network)$size+1, vertex.frame.color="grey30",
       vertex.label=V(int_network)$label, vertex.label.color="black",
       vertex.label.cex=0.5,
       layout=l, main = main)
  
  # plot legend
  x <- sort(table(df_int$interaction[!df_int$interaction %in% 'Bait-Prey']))
  labels <- names(x)
  labels <- paste0(labels, ' [',x,' edges]')
  labels <- gsub("\\s*\\([^\\)]+\\)","", labels)
  labels <- tail(labels, n = 10)
  legend("bottomright",paste(rev(labels), collapse = '\n'), bty="n", cex = 0.7)
  
}






















