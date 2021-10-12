#' @title duplicated interaction
#' @description Assuming an uni-directed graph, this function will return TRUE
#' for all entries that are duplicated.
#' @param x a data.frame or matrix with at least two columns representing edges,
#' e.g 'Bait','Prey'
#' @family graphs
#' @export

duplicated_interaction <- function(x){
  return(duplicated(t((apply(x, 1, sort)))))
}

#' @title lower triangle dupllicated
#' @description used to deal with a post-merge data.frame, 
#' and find unidirectional comparisons.
#' @param x a vector of strings
#' @param split how the string should be splitted
#' @export

# function for finding duplicated names
lower_triangle_duplicated <- function(x, split = '\\.'){
  return(duplicated(unlist(lapply(strsplit(x, split = split), function(x) paste(sort(x), collapse = '.')))))
}

#' @title index 
#' @description index a vector by a splitted string.
#' @param name string or vector of strings.
#' @param i index
#' @param split how should the strings be splitted.
#' @export
# for quickly indexing
index <- function(name, i, split = '_'){unlist(lapply(strsplit(name, split = split), function(x) x[i]))}


#' @title as interactor complex
#' @description Converts a vector of N strings
#' to a 2 x (N ^ 2) data.frame that indicates 
#' a complex. I.e. all proteins interactor 
#' with all proteins in a directed manner.
#' @note duplicated (directed) entries are returned.
#' @export 
as_interactor_complex <- function(x){
  df <- do.call(rbind, lapply(x, function(i){ do.call(rbind, lapply(x, function(j){ data.frame(int1 = i, int2 = j) })) }))
  return(df)
}

