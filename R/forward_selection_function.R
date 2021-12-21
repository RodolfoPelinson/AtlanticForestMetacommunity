#' @title Automating the forward selection of variables
#'
#' @description This function takes a matrix of variables and excludes those who are not important in explaining species occurrence.
#' If no variables are important, it does not exclude any variables. It uses the function forward.sel() from package adespatial.
#' @param Y A matrix of response variables
#' @param X A matrix of predictor variables from  variables will be selected.
#' @export
#'

forward_selection <- function(Y,X){

  New_X <- X[is.na(rowSums(X)) == FALSE,]
  New_Y <- Y[is.na(rowSums(X)) == FALSE,]

  X_rda <- rda(New_Y,New_X)
  p <- anova.cca(X_rda, permutations = 9999)
  if(na.omit(p$`Pr(>F)`) <= 0.05){
    X.R2 <- RsquareAdj(X_rda)$adj.r.squared

    res <- try(X_forw <- forward.sel(New_Y, as.matrix(New_X), adjR2thresh = X.R2, nperm = 9999))
    if(inherits(res, "try-error")){
      message("No variables selected")
      New_X_2 <- data.frame(X)
    }else{
      New_X_2 <- data.frame(X[,match(X_forw$variables, colnames(X))])
      colnames(New_X_2) <- colnames(X)[1:nrow(X_forw)]
    }

  }
  else{
    message("Forward selection NOT performed. p > 0.05")
    New_X_2 <- data.frame(X)
    X_forw <- p
  }
  if(is.null(New_X_2)){New_X_2 <- data.frame(X)}
  result <- list(forward_results = X_forw,selected_variables = New_X_2)
  return(result)
}
