#' @title Automating the exclusion of multicolinear variables
#'
#' @description This function takes a matrix of variables and excludes those who exhibit
#' a variance inflation factor higher than 3 for giving matrix of response variables. It uses the function vif.cca() from package vegan.
#' @param Y A matrix of response variables
#' @param X A matrix of predictor variables from which multicolinear variables will be excluded.
#' @export
#'

VIF_selection <- function(Y,X){
  New_X <- X
  for(i in 1:length(colnames(X))){
    XVIF <- vif.cca(rda(Y,New_X))
    delete <-NA
    higher_vif <- max(XVIF)

    if(higher_vif>3){
      delete <- match(higher_vif, XVIF)
    }

    if(is.na(delete)){break}else{
      New_X <- New_X[,-delete]
    }
  }
  return(New_X)
}
