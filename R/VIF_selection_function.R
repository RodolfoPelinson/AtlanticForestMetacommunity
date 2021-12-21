#' @title Automating the exclusion of multicolinear variables
#'
#' @description This function takes a matrix of variables and excludes those who exhibit
#' a variance inflation factor higher than 3 for giving matrix of response variables. It uses the function vif.cca() from package vegan.
#' @param Y A matrix of response variables
#' @param X A matrix of predictor variables from which multicolinear variables will be excluded.
#' @export
#'


VIF_selection <- function(Y,X){

New_X <- X[is.na(rowSums(X)) == FALSE,]
New_Y <- Y[is.na(rowSums(X)) == FALSE,]
vifs <- list()

  for(i in 1:length(colnames(X))){
    XVIF <- vif.cca(rda(New_Y,New_X, na.rm = TRUE))
    vifs[[i]] <- XVIF
    delete <-NA
    higher_vif <- max(XVIF, na.rm = TRUE)

    if(higher_vif>3){
      delete <- match(higher_vif, XVIF)
    }

    if(is.na(delete)){break}else{
      New_X <- New_X[,-delete]
      X <- X[,-delete]

    }
    vifs[[i]] <- XVIF
  }
  return(list(variables = X,VIFs = vifs))
}
