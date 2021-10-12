#' @title Impute values
#' @description Replacing missing values with randomly sampled values from normal distribution,
#' with width SD x width and down-shifted Median-Sd x shift compared to observed sample distribution.
#' This is building upon the assumption that missing values have arisen due to low expression that 
#' can't be quantified. Therfore, shifting the median to lower expression levels will provide a 
#' proxy of this. In contrast to, locf imputation, this ensures that the variance is not reduced
#' which would consequently impact the moderated t.test.
#' @param df a data.frame with numeric columns
#' @param width numeric. change the factor of the standard deviation.
#' @param shift numeric. Negative values will shift the median distribution downwards.
#' @param seed numeric. Random number seed to be used in the analysis.
#' @note No down-shifting and stdwith of 0.5 do not simualte low abudant missing values.
#' down-shifting of 0.8 and stdwidth of 0.5 simulates low abundant missing values. 
#' down-shifting of 3.6 and stdwith of 0.5 results in (usually undesired) bi-modal distribution.
#' 
#' @references (Perseus, Tyanova et al. 2016)
#' @return data.frame with missing values imputed.
#' @family processing
#' @export

impute_by_col <- function(df, width = 0.3, shift = -1.8, seed = 4295, verbose = T){
  
  set.seed(seed)
  
  numerics <- as.logical(sapply(df, function(x) is.numeric(x)))
  cols <- as.vector(unlist(lapply(df, function(x) is.numeric(x) & any(is.na(x)))))
  #df$imputed <- as.logical(apply(df[numerics], 1, function(x) any(is.na(x)))) 
  
  impute <- function(x){
    std <- sd(x, na.rm = T) * width # adjsuted/down-shifted mean (sample mean - SD * shift)
    me <- median(x, na.rm = T) + sd(x, na.rm = T) * shift  # adjsuted SD (sample SD * width)
    x[is.na(x)] <- rnorm(sum(as.numeric(is.na(x))), mean = me, sd = std)
    return(x)
  }
  
  #impute(df$)
  
  nimputed <- sum(is.na(df[,numerics]))
  if (verbose & nimputed > 0) write(paste('[impute] imputed',nimputed, 'value(s).'),stderr())
  if (sum(cols) == 1) df[, cols] <- impute(df[,cols & numerics])
  if (sum(cols) > 1) df[, cols] <- lapply(df[,cols & numerics], impute)
  return(df)
}
