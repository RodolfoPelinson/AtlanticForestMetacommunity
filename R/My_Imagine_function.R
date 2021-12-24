#' @title Plot presence absence matirx with up to five associated environmental gradients.
#'
#' @description This function plots the incidence matrix, ordinated or not, with up to five associated environmental gradients associates.
#' @param comm An incidence matrix.
#' @param col The colors of absences, true presences and filled presences, respectively.
#' @param order Logical. Should the matrix be order by reciprocal averaging (CA)? Defaut set is to TRUE.
#' @param scores Numeric. What axis of the CA should be used for ordination? Defaut is 1.
#' @param fill Logical. Should embeded absences be filled?
#' @param xlab X axis Label. Defaut is to "Species".
#' @param ylab Y axis Label. Defaut is to "Sites".
#' @param yline Line position of y label. Defaut is 2.
#' @param xline Line position of x label. Defaut is 2.
#' @param sitenames A vector of site Names. Defaut is row names.
#' @param speciesnames A vector of species Names. Defaut is column names.
#' @param Env1 The first environmental gradient
#' @param Env.label_1 Name of the first environmental gradient.
#' @param Env.col_1 The colors of the first Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.
#' @param Env2 The second environmental gradient.
#' @param Env.label_2 Name of the second environmental gradient.
#' @param Env.col_2 The colors of the second Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.
#' @param Env3 The third environmental gradient.
#' @param Env.label_3 Name of the third environmental gradient.
#' @param Env.col_3 The colors of the third Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.
#' @param Env4 The fourth environmental gradient.
#' @param Env.label_4 Name of the fourth environmental gradient.
#' @param Env.col_4 The colors of the fourth Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.
#' @param Env5 The fifth environmental gradient.
#' @param Env.label_5 Name of the fifth environmental gradient.
#' @param Env.col_5 The colors of the fifth Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.
#' @param cex.site The size of the site labels.
#' @param cex.species The size of the species labels.
#' @param cex.envlab The size of the environmental gradient labels
#' @param top_margin The space left in the top of the plot. Useful to fit long species names.
#' @param left_margin The space left in the left side of the plot. Useful to fit long site names.
#' @param bottom_margin The space left in the bottom side of the plot.
#' @param right_margin The space left in the right side of the plot.
#' @param EigenVal Should the relative eigenvalue of the ordination axis used be plotted?
#' @param xlab Label of the x axis
#' @param ylab Label of the y axis
#' @param xlab_line Line position of the x label
#' @param ylab_line Line position of the y label
#' @param cex.xlab Size of the x label
#' @param cex.ylab Size of the y label
#' @param main Main Title
#' @param main_line Line poisition of the main title.
#' @param cex.main Size of the main title
#' @param box.lwd The width of the plot box.
#' @export
#'


My_Imagine <- function (comm, col = c(0,1,2), order = TRUE, scores = 1, fill = TRUE,
                        xlab = "", ylab = "", yline = 0, xline = 0, sitenames = rownames(comm), cex.envlab = 1,
                        speciesnames = colnames(comm),   xlab_line = 3, cex.xlab = 2, ylab_line = 3, cex.ylab = 2, main = "", main_line = 3, cex.main = 2,
                        Env1 = NULL, Env2 = NULL, Env.col_1 = NULL, Env.label_1 = NULL, Env.col_2 = NULL, Env.label_2 = NULL,
                        cex.site = 1, cex.species = 1, top_margin = 2, left_margin = 2, bottom_margin = 3, right_margin = 1,EigenVal = F, box.lwd = 1,
                        Env.col_3 = NULL, Env.label_3 = NULL, Env3 = NULL,
                        Env.col_4 = NULL, Env.label_4 = NULL, Env4 = NULL,
                        Env.col_5 = NULL, Env.label_5 = NULL, Env5 = NULL)



{
  if(isFALSE(is.null(Env1))  & is.null(Env2) & is.null(Env3)& is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths = c(0.95,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & is.null(Env3)& is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), widths = c(0.9,0.05,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & is.null(Env4)& is.null(Env5)){
    layout(matrix(c(1,2,3,4), 1, 4, byrow = TRUE), widths = c(0.85,0.05,0.05,0.05))
  }
  if(isFALSE(is.null(Env1)) & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & isFALSE(is.null(Env4)) & is.null(Env5)){
    layout(matrix(c(1,2,3,4,5), 1, 5, byrow = TRUE), widths = c(0.868,0.033,0.033,0.033,0.033))
  }
  if(isFALSE(is.null(Env1))  & isFALSE(is.null(Env2)) & isFALSE(is.null(Env3)) & isFALSE(is.null(Env4)) & isFALSE(is.null(Env5))){
    layout(matrix(c(1,2,3,4,5,6), 1, 6, byrow = TRUE), widths = c(0.835,0.033,0.033,0.033,0.033,0.033))
  }

  require(metacom)
  if (order == TRUE) {
    new_scores <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = T)
    if(is.null(sitenames)==FALSE){
      sitenames <- sitenames[order(new_scores$sitescores)]
    }
    if(is.null(speciesnames)==FALSE){
      speciesnames <- speciesnames[order(new_scores$speciesscores)]
    }
    if(is.null(Env1) == FALSE){
      Env1 <- Env1[order(new_scores$sitescores)]
    }
    if(is.null(Env2) == FALSE){
      Env2 <- Env2[order(new_scores$sitescores)]
    }

    if(is.null(Env3) == FALSE){
      Env3 <- Env3[order(new_scores$sitescores)]
    }

    if(is.null(Env4) == FALSE){
      Env4 <- Env4[order(new_scores$sitescores)]
    }

    if(is.null(Env5) == FALSE){
      Env5 <- Env5[order(new_scores$sitescores)]
    }



    EigenValues <- decorana(comm, ira = 1)$evals
    Relative_EigenValues <- EigenValues/sum(EigenValues)
    comm <- OrderMatrix(comm, binary = TRUE, scores = scores, outputScores = FALSE)

  }


  #Filling embedded absences
  if (fill == TRUE) {
    temp1 = comm
    for (i in 1:dim(comm)[2]) {
      temp2 = comm[, i]
      if (sum(temp2) < 2) {
        comm[, i] = temp2
      }
      else {
        first <- min(which(temp2 > 0))
        last <- max(which(temp2 > 0))
        temp1[first:last, i] <- 1
      }
    }
    for (i in 1:dim(comm)[1]){
      for (j in 1:dim(comm)[2]){
        if(temp1[i,j] == 1 & comm[i,j] == 0){
          comm[i,j] <- 2
        }
      }
    }
  }

  #Reversing the order of rows?
  reverse <- nrow(comm):1
  comm <- comm[reverse, ]

  #ploting
  par(mar = c(bottom_margin,left_margin, top_margin, right_margin), cex = 1)
  image(1:dim(comm)[2], 1:dim(comm)[1], t(comm),  col = col,
        xlab = "", ylab = "", axes = FALSE)
  title(xlab = xlab, line = xlab_line, cex.lab = cex.xlab, adj = 1)
  title(ylab = ylab, line = ylab_line, cex.lab = cex.ylab)
  title(main = main, line = main_line, cex.main = cex.main, adj = 0, font.main = 1)


  box(lwd = box.lwd)
  if (length(sitenames) > 1) {
    axis(2, at = 1:dim(comm)[1], labels = sitenames, las = 1,
         cex.axis = cex.site, tick = FALSE, line = yline)
  }
  if (length(speciesnames) > 1) {
    axis(3, at = 1:dim(comm)[2], labels = speciesnames, las = 2,
         cex.axis = cex.species, tick = FALSE, line = xline)
  }

  if(is.null(Env1)==F){

    if(is.numeric(Env1)){
      pal <- col_numeric(
        palette = c(Env.col_1[1], Env.col_1[2]),
        domain = Env1,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_1 <-pal(Env1)
    }

    if(is.factor(Env1)){
      Env1<-as.numeric(Env1)
    }
    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env1),x = 1, z = t(as.matrix(1:length(Env1))),
          col = Env.col_1,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_1, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)

  }

  if(is.null(Env2)==F){

    if(is.numeric(Env2)){
      pal <- col_numeric(
        palette = c(Env.col_2[1], Env.col_2[2]),
        domain = Env2,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_2 <-pal(Env2)
    }

    if(is.factor(Env2)){
      Env2<-as.numeric(Env2)
    }

    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env2),x = 1, z = t(as.matrix(1:length(Env2))),
          col = Env.col_2,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_2, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }

  if(is.null(Env3)==F){

    if(is.numeric(Env3)){
      pal <- col_numeric(
        palette = c(Env.col_3[1], Env.col_3[2]),
        domain = Env3,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_3 <-pal(Env3)
    }

    if(is.factor(Env3)){
      Env3<-as.numeric(Env3)
    }

    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env3),x = 1, z = t(as.matrix(1:length(Env3))),
          col = Env.col_3,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_3, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }


  if(is.null(Env4)==F){

    if(is.numeric(Env4)){
      pal <- col_numeric(
        palette = c(Env.col_4[1], Env.col_4[2]),
        domain = Env4,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE
      )
      Env.col_4 <- pal(Env4)
    }

    if(is.factor(Env4)){
      Env4<-as.numeric(Env4)
    }

    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env4),x = 1, z = t(as.matrix(1:length(Env4))),
          col = Env.col_4,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_4, las = 2,
         cex.axis = cex.envlab, tick = FALSE, line = xline)
    box(lwd = box.lwd)
  }


  if(is.null(Env5)==F){

    if(is.numeric(Env5)){
      pal <- col_numeric(
        palette = c(Env.col_5[1], Env.col_5[2]),
        domain = Env5,
        na.color = "grey50",
        alpha = FALSE,
        reverse = FALSE)
      Env.col_5 <-pal(Env5)
    }

    if(is.factor(Env5)){
      Env5<-as.numeric(Env5)
    }

    par(mar = c(bottom_margin,0, top_margin, right_margin))
    image(y = 1:length(Env5),x = 1, z = t(as.matrix(1:length(Env5))),
          col = Env.col_5,
          axes = FALSE, xlab = "", ylab = "", xlim = c(0,1))
    axis(3, at = 0.4, labels = Env.label_5, las = 2,
         cex.axis = cex.envlab, tick = FALSE,  line = xline)
    box(lwd = box.lwd)
  }



  if(isTRUE(EigenVal)){
    y <- -((dim(comm)[2]^2)/10000)
    text(x=(dim(comm)[2]),y=y, labels = paste(round(Relative_EigenValues[scores]*100,2),"%", sep = " "), xpd=NA, adj=c(1,1))
  }

  if(isTRUE(EigenVal)){
    y <- -((dim(comm)[2]^2)/10000)
    text(x=(dim(comm)[2]),y=y, labels = paste(round(Relative_EigenValues[scores]*100,2),"%", sep = " "), xpd=NA, adj=c(1,1))
  }

  layout(matrix(c(1), 1, 1))
  par(mar = c(4,4,4,4))
}
