#' @title Run the EMS analysis and identify the idealized structure
#'
#' @description This function simply runs the functions Coherence(), Turnover() and BoundaryClump() from package metacom and automaticaly identify the idealizes structure according to Presley et al. 2010 "A comprehensive framework for the evaluation of metacommunity structure". Oikos 119 (6).
#' @param comm A list of incidence matrices (i.e. metacommunities).
#' @param names The names of the metacommunities. A vector with the same length as the list of metacommunities.
#' @param scores Numeric. What axis of the CA should be used for ordination? (Defaut is 1)
#' @param CoherenceMethod null model randomization method used by 'nullmaker' to compute Coherence. See the Coherence function from package metacom. (default is "curveball")
#' @param turnoverMethod null model randomization method used by 'nullmaker' or 'EMS' to use the approach outlined in Leibold and Mikkelson 2002. See the Turnover function from package metacom. (default is "EMS")
#' @param sims Number of randomizations (default is 1000)
#' @param order Should the original matrix be ordered by reciprocal averaging?
#' @param orderNulls Should the null communities be ordered? (default is TRUE)
#' @param seed seed for simulating the null model. Null matrices should be repeatable.
#' @param fill should embedded absences be filled before the statistics are calculated? (default is TRUE)
#' @param round Numeric. Should numeric results be rounded? If so, how many digits? Defaut is set to NULL.
#' @param elapsed_time Logical. Should a message with the elapsed time be returned?

#' @export
#'

IdentifyStructure <- function(comm, names = NULL, scores = 1,CoherenceMethod = "curveball", turnoverMethod = "EMS", sims = 1000, order = T, orderNulls = F,  seed = NULL, fill = T, round = NULL, elapsed_time = TRUE){

  if(isTRUE(elapsed_time)){start_time <- Sys.time()}

  if(is.null(names)){names <- c(1:length(comm))}

  Embeded_Absences <- rep(NA, length(comm))
  Simulated_Embeded_Absences <- rep(NA, length(comm))
  percent_difference_EmbAbs <- rep(NA, length(comm))
  z_Coherence <- rep(NA, length(comm))
  p_Coherence <- rep(NA, length(comm))

  Turnover <- rep(NA, length(comm))
  Simulated_Turnover <- rep(NA, length(comm))
  percent_difference_Turn <- rep(NA, length(comm))
  z_Turnover = rep(NA, length(comm))
  p_Turnover <- rep(NA, length(comm))

  I_Index <- rep(NA, length(comm))
  p_I_Index <- rep(NA, length(comm))
  N_sites <- rep(NA, length(comm))
  N_species = rep(NA, length(comm))
  Structure <- rep(NA, length(comm))

  pb <- txtProgressBar(min = 0, max = length(comm), style = 3)


  for (i in 1:length(comm)){
    #  Sys.sleep(0.1)

    tryCatch({
      metacom_cohe <- metacom::Coherence(comm=comm[[i]], scores=scores, method=CoherenceMethod, sims=sims, order=order, allowEmpty=FALSE, verbose = F,  orderNulls = orderNulls, seed = seed)
      metacom_turn <- metacom::Turnover(comm=comm[[i]], scores=scores, method= turnoverMethod, sims=sims, order=order, allowEmpty=FALSE, verbose = F,  orderNulls = orderNulls, seed = seed, fill = fill)
      metacom_bound <- metacom::BoundaryClump(comm=comm[[i]], scores=scores, order=order)


      Quasi <- "Quasi"
      Loss <- NULL

      embAbs <- metacom_cohe$stat[1]
      z_embAbs <- metacom_cohe$stat[2]
      null_embAbs <- metacom_cohe$stat[4]
      Coe_p <- metacom_cohe$stat[3]

      turnover <- metacom_turn$stat[1]
      z_turn <-metacom_turn$stat[2]
      null_Turnover <- metacom_turn$stat[4]
      Turnover_p <- metacom_turn$stat[3]

      Boundary <- metacom_bound$stat[1]
      Boundary_p <- metacom_bound$stat[2]

      if(is.nan(Coe_p)){
        if(turnover<null_Turnover){structure <- "Nested"}
        if(turnover>null_Turnover ){
          if(Boundary_p > 0.05){structure <- "Gleasonian"}
          if(Boundary_p <= 0.05 & Boundary > 1){structure <- "Clementsian"}
          if(Boundary_p <= 0.05 & Boundary < 1){structure <- "Evenly-Spaced"}
        }
      }else{
        if(Coe_p > 0.05){structure <- "Random"}
        if(Coe_p <= 0.05){
          if(embAbs < null_embAbs){
            if(turnover<null_Turnover){structure <- "Nested"}
            if(turnover>null_Turnover ){
              if(Boundary_p > 0.05){structure <- "Gleasonian"}
              if(Boundary_p <= 0.05 & Boundary > 1){structure <- "Clementsian"}
              if(Boundary_p <= 0.05 & Boundary < 1){structure <- "Evenly-Spaced"}
            }
          }
        }
      }


      if(is.nan(Turnover_p)==F){
        if(is.nan(Coe_p)){
          if(Turnover_p>0.05){structure <- paste(Quasi, structure, sep = "-")}
        }else{
          if(Turnover_p>0.05 & Coe_p < 0.05){structure <- paste(Quasi, structure, sep = "-")}
        }
      }

      if(is.nan(Turnover_p) & is.nan(Coe_p)){structure <- "Nested"}

      if(Boundary_p > 0.05){Loss <- "with random species loss"}
      if(Boundary_p <= 0.05 & Boundary > 1){Loss <- "with clumped species loss"}
      if(Boundary_p <= 0.05 & Boundary < 1){Loss <- "with hyperdispersed species loss"}

      if(is.nan(Coe_p)){
        if(turnover <= null_Turnover){structure <- paste(structure, Loss, sep = " ")}
      }else{
        if(turnover <= null_Turnover & Coe_p <= 0.05){structure <- paste(structure, Loss, sep = " ")}
      }

      if(Coe_p <= 0.05){
        if(embAbs > null_embAbs){
          if(embAbs > null_embAbs){structure <- "Checkerboard"}
        }
      }


      #Now we compute other useful information
      percent_emb_abs <- 1-(embAbs/null_embAbs)
      percent_turnover <- 1-(turnover/null_Turnover)



      Embeded_Absences[i] = embAbs
      Simulated_Embeded_Absences[i] = null_embAbs
      percent_difference_EmbAbs[i] = percent_emb_abs
      z_Coherence[i] = z_embAbs
      p_Coherence[i] = Coe_p

      Turnover[i] = turnover
      Simulated_Turnover[i] = null_Turnover
      percent_difference_Turn[i] = percent_turnover
      z_Turnover[i] = z_turn
      p_Turnover[i] = Turnover_p

      I_Index[i] = Boundary
      p_I_Index[i] = Boundary_p

      N_sites[i] <- nrow(comm[[i]])
      N_species[i] = ncol(comm[[i]])

      Structure[i] = structure

      # update progress bar
      setTxtProgressBar(pb, i)

    }, error=function(e){})



  }
  close(pb)


  results <- data.frame(Embeded_Absences,
                        Simulated_Embeded_Absences,
                        percent_difference_EmbAbs,
                        z_Coherence,
                        p_Coherence,
                        Turnover,
                        Simulated_Turnover,
                        percent_difference_Turn,
                        z_Turnover,
                        p_Turnover,
                        I_Index,
                        p_I_Index,
                        N_sites,
                        N_species)
  if(is.null(round) == FALSE){
    results <- round(results, round)
  }

  results <- data.frame(results, Structure)
  rownames(results) <- names

  if(isTRUE(elapsed_time)){
    end_time <- Sys.time()
    message(paste("Elapsed time: ", round(end_time - start_time, 3) ) )
    }



  return(results)
}
