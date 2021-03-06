% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/My_Imagine_function.R
\name{My_Imagine}
\alias{My_Imagine}
\title{Plot presence absence matirx with up to five associated environmental gradients.}
\usage{
My_Imagine(
  comm,
  col = c(0, 1, 2),
  order = TRUE,
  scores = 1,
  fill = TRUE,
  xlab = "",
  ylab = "",
  yline = 0,
  xline = 0,
  sitenames = rownames(comm),
  cex.envlab = 1,
  speciesnames = colnames(comm),
  xlab_line = 3,
  cex.xlab = 2,
  ylab_line = 3,
  cex.ylab = 2,
  main = "",
  main_line = 3,
  cex.main = 2,
  Env1 = NULL,
  Env2 = NULL,
  Env.col_1 = NULL,
  Env.label_1 = NULL,
  Env.col_2 = NULL,
  Env.label_2 = NULL,
  cex.site = 1,
  cex.species = 1,
  top_margin = 2,
  left_margin = 2,
  bottom_margin = 3,
  right_margin = 1,
  EigenVal = F,
  box.lwd = 1,
  Env.col_3 = NULL,
  Env.label_3 = NULL,
  Env3 = NULL,
  Env.col_4 = NULL,
  Env.label_4 = NULL,
  Env4 = NULL,
  Env.col_5 = NULL,
  Env.label_5 = NULL,
  Env5 = NULL
)
}
\arguments{
\item{comm}{An incidence matrix.}

\item{col}{The colors of absences, true presences and filled presences, respectively.}

\item{order}{Logical. Should the matrix be order by reciprocal averaging (CA)? Defaut set is to TRUE.}

\item{scores}{Numeric. What axis of the CA should be used for ordination? Defaut is 1.}

\item{fill}{Logical. Should embeded absences be filled?}

\item{xlab}{Label of the x axis}

\item{ylab}{Label of the y axis}

\item{yline}{Line position of y label. Defaut is 2.}

\item{xline}{Line position of x label. Defaut is 2.}

\item{sitenames}{A vector of site Names. Defaut is row names.}

\item{cex.envlab}{The size of the environmental gradient labels}

\item{speciesnames}{A vector of species Names. Defaut is column names.}

\item{xlab_line}{Line position of the x label}

\item{cex.xlab}{Size of the x label}

\item{ylab_line}{Line position of the y label}

\item{cex.ylab}{Size of the y label}

\item{main}{Main Title}

\item{main_line}{Line poisition of the main title.}

\item{cex.main}{Size of the main title}

\item{Env1}{The first environmental gradient}

\item{Env2}{The second environmental gradient.}

\item{Env.col_1}{The colors of the first Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.}

\item{Env.label_1}{Name of the first environmental gradient.}

\item{Env.col_2}{The colors of the second Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.}

\item{Env.label_2}{Name of the second environmental gradient.}

\item{cex.site}{The size of the site labels.}

\item{cex.species}{The size of the species labels.}

\item{top_margin}{The space left in the top of the plot. Useful to fit long species names.}

\item{left_margin}{The space left in the left side of the plot. Useful to fit long site names.}

\item{bottom_margin}{The space left in the bottom side of the plot.}

\item{right_margin}{The space left in the right side of the plot.}

\item{EigenVal}{Should the relative eigenvalue of the ordination axis used be plotted?}

\item{box.lwd}{The width of the plot box.}

\item{Env.col_3}{The colors of the third Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.}

\item{Env.label_3}{Name of the third environmental gradient.}

\item{Env3}{The third environmental gradient.}

\item{Env.col_4}{The colors of the fourth Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.}

\item{Env.label_4}{Name of the fourth environmental gradient.}

\item{Env4}{The fourth environmental gradient.}

\item{Env.col_5}{The colors of the fifth Environmental Gradient. If it is a continuous gradient, this should be a vector with the two most extreme colors. If it is a factor, it should be a vector with one color for each of the levels of the factor.}

\item{Env.label_5}{Name of the fifth environmental gradient.}

\item{Env5}{The fifth environmental gradient.}
}
\description{
This function plots the incidence matrix, ordinated or not, with up to five associated environmental gradients associates.
}
