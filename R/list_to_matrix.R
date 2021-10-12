# borrowed from here: https://rdrr.io/github/jokergoo/ComplexHeatmap/src/R/Upset.R
# since R-package is not available

list_to_matrix = function(lt, universal_set = NULL) {
  if(!is.null(universal_set)) {
    lt = lapply(lt, function(x) intersect(x, universal_set))
  } else {
    universal_set = unique(unlist(lt))
  }
  
  mat = matrix(0, nrow = length(universal_set), ncol = length(lt))
  rownames(mat) = sort(universal_set)
  colnames(mat) = names(lt)
  for(i in seq_along(lt)) {
    mat[as.character(unique(lt[[i]])), i] = 1
  }
  return(mat)
}