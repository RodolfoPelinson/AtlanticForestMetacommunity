#' @title Run the EMS analysis and identify the idealized structure
#'
#' @description This function is the same as IdentifyStructure but uses apply functions intead of for loops. It enables multicore processing, which may be useful if you have too many metacommunties (> 100). It simply runs the functions Coherence(), Turnover() and BoundaryClump() from package metacom and automaticaly identify the idealizes structure according to Presley et al. 2010 "A comprehensive framework for the evaluation of metacommunity structure". Oikos 119 (6).
#' @param comm A list of incidence matrices (i.e. metacommunities).
#' @param names The names of the metacommunities. A vector with the same length as the list of metacommunities.
#' @param scores Numeric. What axis of the CA should be used for ordination? (Defaut is 1)
#' @param CoherenceMethod null model randomization method used by 'nullmaker' to compute Coherence. See the Coherence function from package metacom. (default is "curveball")
#' @param turnoverMethod null model randomization method used by 'nullmaker' or 'EMS' to use the approach outlined in Leibold and Mikkelson 2002. See the Turnover function from package metacom. (default is "EMS")
#' @param allow_Checkerboard Logical. Should Checkerboard structures be Identified? Presley et al. 2019 advocate that checkerboards should be restricted to small sets of ecologically similar species for which interspecific interactions may lead to mutual exclusion.
#' @param sims Number of randomizations (default is 1000)
#' @param order Should the original matrix be ordered by reciprocal averaging?
#' @param orderNulls Should the null communities be ordered? (default is TRUE)
#' @param seed seed for simulating the null model. Null matrices should be repeatable.
#' @param fill should embedded absences be filled before the statistics are calculated? (default is TRUE)
#' @param round Numeric. Should numeric results be rounded? If so, how many digits? Defaut is set to NULL.
#' @param elapsed_time Logical. Should a message with the elapsed time be returned?
#' @param multicore Logical. Should multicore processing be enabled? The default number of cores is the total number of cores minus 1. Alternative number of cores can be set by the n_cores argument.
#' @param n_cores Number of cores to be used with multicore processing. It will be ignored if multicore is set to FALSE.


#' @export
#'


IdentifyStructure2 <- function(comm, names = NULL, scores = 1,
                               CoherenceMethod = "curveball",
                               turnoverMethod = "EMS",
                               sims = 1000,
                               order = T,
                               orderNulls = F,
                               seed = NULL,
                               fill = T,
                               round = NULL,
                               elapsed_time = TRUE,
                               multicore = FALSE,
                               n_cores = NULL,
                               allow_Checkerboard = FALSE){

  require(pbapply)

  if(isTRUE(elapsed_time)){start_time <- Sys.time()}

  if(is.null(names)){names <- c(1:length(comm))}


  Identify <- function(comm, ...){

    results <- c(Embeded_Absences = NA,
                 Simulated_Embeded_Absences = NA,
                 percent_difference_EmbAbs = NA,
                 z_Coherence = NA,
                 p_Coherence = NA,
                 Turnover = NA,
                 Simulated_Turnover = NA,
                 percent_difference_Turn = NA,
                 z_Turnover = NA,
                 p_Turnover = NA,
                 I_Index = NA,
                 p_I_Index = NA,
                 N_sites = NA,
                 N_species = NA,
                 Structure = NA)

    tryCatch({
      metacom_cohe <- metacom::Coherence(comm=comm, scores=scores, method=CoherenceMethod, sims=sims, order=order, allowEmpty=FALSE, verbose = F,  orderNulls = orderNulls, seed = seed)
      metacom_turn <- metacom::Turnover(comm=comm, scores=scores, method= turnoverMethod, sims=sims, order=order, allowEmpty=FALSE, verbose = F,  orderNulls = orderNulls, seed = seed, fill = fill)
      metacom_bound <- metacom::BoundaryClump(comm=comm, scores=scores, order=order)


      Quasi <- "Quasi"
      Loss <- NULL

      if(is.null(metacom_cohe) | is.null(metacom_turn) | is.null(metacom_bound)){stop()}

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

      if(isTRUE(allow_Checkerboard)){
        if(Coe_p <= 0.05){
          if(embAbs > null_embAbs){
            if(embAbs > null_embAbs){structure <- "Checkerboard"}
          }
        }
      }else{
        if(Coe_p <= 0.05){
          if(embAbs > null_embAbs){
            if(embAbs > null_embAbs){structure <- "Random"}
          }
        }
      }



      #Now we compute other useful information
      percent_emb_abs <- 1-(embAbs/null_embAbs)
      percent_turnover <- 1-(turnover/null_Turnover)



      Embeded_Absences = embAbs
      Simulated_Embeded_Absences = null_embAbs
      percent_difference_EmbAbs = percent_emb_abs
      z_Coherence = z_embAbs
      p_Coherence = Coe_p
      Turnover = turnover
      Simulated_Turnover = null_Turnover
      percent_difference_Turn = percent_turnover
      z_Turnover = z_turn
      p_Turnover = Turnover_p
      I_Index = Boundary
      p_I_Index = Boundary_p
      N_sites <- nrow(comm)
      N_species = ncol(comm)
      Structure = structure

      results <- c(Embeded_Absences = Embeded_Absences,
                   Simulated_Embeded_Absences = Simulated_Embeded_Absences,
                   percent_difference_EmbAbs = percent_difference_EmbAbs,
                   z_Coherence = z_Coherence,
                   p_Coherence = p_Coherence,
                   Turnover = Turnover,
                   Simulated_Turnover = Simulated_Turnover,
                   percent_difference_Turn = percent_difference_Turn,
                   z_Turnover = z_Turnover,
                   p_Turnover = p_Turnover,
                   I_Index = I_Index,
                   p_I_Index = p_I_Index,
                   N_sites = N_sites,
                   N_species = N_species)

      if(is.null(round) == FALSE){
        results <- round(results, round)
      }
      results <- c(results, Structure = Structure)

    }, error=function(e){})

    return(results)
  }


  if(isTRUE(multicore)){
    require(doParallel)
    if(is.null(n_cores)){
      n_cores <- detectCores() - 1
      message(paste("Using ", n_cores, " cores"))
      registerDoParallel(cores=n_cores)
      cl <- makeCluster(n_cores, type="PSOCK")
      results <- pblapply(comm, Identify, cl = cl)
      #results <- parLapply(cl, comm, Identify)
      stopCluster(cl)
    }else{
      if(detectCores() >= n_cores){
        message(paste("Using ", n_cores, " cores"))
        registerDoParallel(cores=n_cores)
        cl <- makeCluster(n_cores, type="PSOCK")
        results <- pblapply(comm, Identify, cl = cl)
        #results <- parLapply(cl, comm, Identify)
        stopCluster(cl)
      }else{
        message(paste("n_cores > total number of cores. Proceeding with single core processing"))
        results <- pblapply(comm, Identify)
      }
    }
  }else{
    results <- pblapply(comm, Identify)
  }




  results <- data.frame(do.call(rbind,results))

  results <- data.frame(apply(results[,1:(ncol(results)-1)], 2, as.numeric),  Structure = results[,ncol(results)])
  rownames(results) <- names


  if(isTRUE(elapsed_time)){
    end_time <- Sys.time()
    message(paste("Elapsed time: ", round(end_time - start_time, 3) ) )
  }


  return(results)
}


