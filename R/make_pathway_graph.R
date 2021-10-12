
make_pathway_graph <- function(df, geneset, collapse_overlap = '\n', remove_unassigned = F){
  
  # prepare geneset
  #geneset <- list(X = c('PMM2', 'SPECC1','DIMT1'), Y = c('RNF114', 'OPA1', 'TSG101','DIMT1'))
  geneset <- stack(geneset)
  colnames(geneset) <- c('int2','group')
  
  # prepare data.frame
  colnames(df) <- c('int1','int2')
  df <- df[!duplicated(df),]
  df_merge <- merge(df, geneset, all.x = T)
  df_merge$group <- as.character(df_merge$group)
  df_merge$group[is.na(df_merge$group)] <- 'Not assigned'
  
  # remove unassigned
  if (remove_unassigned){
    df_merge <- df_merge[!df_merge$group %in% 'Not assigned',]
  }
  
  # deal with interactors in multiple groups
  int_count <- table(df_merge$int2)
  df_int <- do.call(rbind, lapply(names(int_count), function(x){
    
    tmp_df <- df_merge[df_merge$int2 %in% x,]
    int1 = unique(tmp_df$int1)
    int2 = unique(tmp_df$int2)
    groups = paste(unique(tmp_df$group), collapse = collapse_overlap)
    return(data.frame(int1, int2, groups))
    
  }))
  
  # make list count 
  unique_groups <- unique(df_int$groups)
  names(unique_groups) <- 1:length(unique_groups)
  
  # prepare groups
  groups <- lapply(unique_groups, function(x){
    return(as.numeric(rownames(df_int[df_int$groups %in% x,])) + 1)
  })
  names(groups) <- unique_groups
  
  # make plot
  qgraph(df_int[,c(1,2)], 
         groups = groups, 
         layout = 'spring', 
         borders = FALSE,
         cut = 0.1,
         legend = T)
  
}






















