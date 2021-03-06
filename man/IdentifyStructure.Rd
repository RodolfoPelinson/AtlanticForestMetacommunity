% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IdentifyStructure_function.R
\name{IdentifyStructure}
\alias{IdentifyStructure}
\title{Run the EMS analysis and identify the idealized structure}
\usage{
IdentifyStructure(
  comm,
  names = NULL,
  scores = 1,
  CoherenceMethod = "curveball",
  turnoverMethod = "EMS",
  sims = 1000,
  order = T,
  orderNulls = F,
  seed = NULL,
  fill = T,
  round = NULL,
  elapsed_time = TRUE
)
}
\arguments{
\item{comm}{A list of incidence matrices (i.e. metacommunities).}

\item{names}{The names of the metacommunities. A vector with the same length as the list of metacommunities.}

\item{scores}{Numeric. What axis of the CA should be used for ordination? (Defaut is 1)}

\item{CoherenceMethod}{null model randomization method used by 'nullmaker' to compute Coherence. See the Coherence function from package metacom. (default is "curveball")}

\item{turnoverMethod}{null model randomization method used by 'nullmaker' or 'EMS' to use the approach outlined in Leibold and Mikkelson 2002. See the Turnover function from package metacom. (default is "EMS")}

\item{sims}{Number of randomizations (default is 1000)}

\item{order}{Should the original matrix be ordered by reciprocal averaging?}

\item{orderNulls}{Should the null communities be ordered? (default is TRUE)}

\item{seed}{seed for simulating the null model. Null matrices should be repeatable.}

\item{fill}{should embedded absences be filled before the statistics are calculated? (default is TRUE)}

\item{round}{Numeric. Should numeric results be rounded? If so, how many digits? Defaut is set to NULL.}

\item{elapsed_time}{Logical. Should a message with the elapsed time be returned?}
}
\description{
This function simply runs the functions Coherence(), Turnover() and BoundaryClump() from package metacom and automaticaly identify the idealizes structure according to Presley et al. 2010 "A comprehensive framework for the evaluation of metacommunity structure". Oikos 119 (6).
}
