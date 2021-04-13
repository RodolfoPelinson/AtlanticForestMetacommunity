#' @title Coefficient of Variation
#'
#' @description Computes de coefficient of variation of a vector
#' @param x vector from which coefficient of variation will be computed

#'
#' @export
#'
#'
#'
#'

coef_var <- function(x){
  result <- sd(x) / mean(x)
  return(result)
}
