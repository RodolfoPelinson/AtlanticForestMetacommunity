#' @title Performs Spearman rank correlations or Kruskal-Wallis Rank Sum Test on the First CA Axis.
#'
#' @description This function performs Spearman rank correlations between the first axis of a CA computed on a incidencematrix and multiple gradient variables.
#' . If the gradient variable is a factor, it performs Kruskal-Wallis Rank Sum Tests. P values are corrected for multiple comparisons.
#' @param Y An incidence matrix.
#' @param X A matrix of gradient variables from, either continuous gradients or a factor.
#' @param p.adjust.method Method of p-value correction. Same as for the "p.adjust" function.
#' @export
#'


My_spearman <- function(Y, X, p.adjust.method = "fdr"){
  Y_scores=OrderMatrix(Y, outputScores=TRUE, scores = 1)$sitescores
  rho <- rep(NA, ncol(X))
  chisq <- rep(NA, ncol(X))

  p <- rep(NA, ncol(X))
  chisq_p <- rep(NA, ncol(X))

  for(i in 1:ncol(X)){
    if(is.factor(X[,i])){
      kruskal <- kruskal.test(Y_scores, X[,i])
      chisq[i] <- kruskal$statistic
      chisq_p[i] <- kruskal$p.value
    }else{
      cor <- cor.test(Y_scores, X[,i], method = "spearman", exact = NULL)
      rho[i] <- cor$estimate
      p[i] <- cor$p.value
    }
  }
  new_p <- rep(NA, ncol(X))
  for(i in 1:ncol(X)){
    if(is.na(p[i])){new_p[i] <- chisq_p[i]}
    else{new_p[i] <- p[i]}
  }

  result <- cbind(rho, `chi-squared` = chisq, p = new_p, adj_p = p.adjust(new_p, method = p.adjust.method))
  rownames(result)<- colnames(X)
  return(result)
}
