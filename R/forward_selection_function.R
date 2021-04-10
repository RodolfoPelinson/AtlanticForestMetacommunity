#' @title Automating the forward selection of variables
#'
#' @description This function takes a matrix of variables and excludes those who are not important in explaining species occurrence.
#' If no variables are important, it does not exclude any variables. It uses the function forward.sel() from package adespatial.
#' @param Y A matrix of response variables
#' @param X A matrix of predictor variables from  variables will be selected.
#' @export
#'

forward_selection <- function(Y,X){
  X_rda <- rda(Y,X)
  p <- anova.cca(X_rda, permutations = 9999)
  if(na.omit(p$`Pr(>F)`) <= 0.05){
    X.R2 <- RsquareAdj(X_rda)$adj.r.squared

    res <- try(X_forw <- forward.sel(Y, as.matrix(X), adjR2thresh = X.R2, nperm = 9999))
    if(inherits(res, "try-error")){
      message("No variables selected")
      NEW_X <- data.frame(X)
    }else{
      NEW_X <- data.frame(X[,match(X_forw$variables, colnames(X))])
      colnames(NEW_X) <- colnames(X)[1:nrow(X_forw)]
    }

  }
  else{
    message("Forward selection NOT performed. p > 0.05")
    NEW_X <- data.frame(X)
  }
  if(is.null(NEW_X)){NEW_X <- data.frame(X)}
  return(NEW_X)
}
