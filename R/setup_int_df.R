#' @title setup interactor group data.frame
#' @description This function will take a data.frame of two two columns (interactor1, interactor2)
#' and a named list and create virtual nodes for each interactor in the group. The returned
#' list will contain the input data.frame + another data.frame with the virtual edges and nodes
#' @param df a data.frame with two columns (interactor1, interactor2)
#' @param groups a named list, where each list name corresponds to a group, and items 
#' in the group are present in the interactor1 or interactor2 column.
#' @return a list of two data.frames
#' @export


setup_int_group_df <- function(df, groups){
  
  stopifnot(ncol(df) == 2)
  stopifnot(!is.data.frame(groups)) # groups must be a list
  stopifnot(is.list(groups))
  
  # prepare data.frame
  colnames(df) <- c('int1','int2')
  df <- df[!duplicated(df),]
  df_int <- df
  
  # setup nodes and edges
  tmp = unique(na.omit(c(as.character(df_int$int1), as.character(df_int$int2))))
  nodes = data.frame(name=tmp, id = 1:length(tmp))
  edges = data.frame(source=df_int$int1, target = df_int$int2, weight = 1)
  edges = edges[complete.cases(edges),]
  
  # setup virtual nodes and edges
  virt_group <- names(groups)
  id_start = max(nodes$id)
  id_end = max(nodes$id)+length(virt_group)-1
  virt_nodes <- data.frame(
    name = virt_group, 
    id = id_start:id_end+1
  )
  
  # iterate through each group and get the edges
  virt_group_edges_lst <- lapply(virt_group, function(x){
    
    groups_genes <- groups[[x]]
    nodes_involved <- groups_genes[groups_genes %in% df$int2]
    edges_new <- NULL
    if (length(nodes_involved) > 0){
      edges_new <- data.frame(source = x, target = nodes_involved, weight = 10)
    } 
    return(edges_new)
  })
  
  # omit NULLs and combine into data.frame
  virt_group_edges_lst <- null_omit(virt_group_edges_lst)
  virt_group_edges <- do.call(rbind, virt_group_edges_lst)
  virt_group_edges <- virt_group_edges[complete.cases(virt_group_edges),]
  
  # combine with normal nodes and edges
  edges <- rbind(edges, virt_group_edges)
  nodes <- rbind(nodes, virt_nodes)
  
  # rename nodes / edges
  nodes$name <- as.character(nodes$name)
  edges$source <- as.character(edges$source)
  edges$target <- as.character(edges$target)
  return(list(edges = edges, nodes = nodes, virt_edges = virt_group_edges, virt_nodes = virt_nodes))
}