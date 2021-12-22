#' @title Performing variation partitioning with p values
#'
#' @description This function uses the `varpart()` function from vegan to compute all fractions and their respective p-values in a single table
#' @param Y Species occurence matrix
#' @param env Matrix of environmental variables
#' @param clim Matrix of climate variables
#' @param spa Matrix of spatial variables
#' @param percent_r2 Logical. Should the R2 values be converted in percentages?
#' @param allow_negative_r2 Logical. Should negative values be used to compute the R2 values?
#'
#' @export
#'


var_partitioning <- function(Y, env, clim = NULL, spa, percent_r2 = F, allow_negative_r2 = TRUE, round = 3){

  #Without climate variables
  if((is.null(clim) || nrow(clim) == 0) & is.null(spa) == FALSE){
    comp.varpart <- varpart(Y, env, spa)

    independent.variables <- data.frame(env,spa)
    #spa_clim <- data.frame(clim,spa)
    #env_clim <- data.frame(clim,env)
    env_spa <- data.frame(spa,env)

    anova_pure_env <- anova.cca(rda(Y, env, spa),permutations=9999)
    anova_pure_spa <- anova.cca(rda(Y, spa, env),permutations=9999)
    #anova_pure_clim <- anova.cca(rda(Y, clim, env_spa),permutations=9999)

    anova_env <- anova.cca(rda(Y, env),permutations=9999)
    anova_spa <- anova.cca(rda(Y, spa),permutations=9999)
    #anova_clim <- anova.cca(rda(Y, clim),permutations=9999)
    anova_all <- anova.cca(rda(Y, independent.variables),permutations=9999)

    p_all <- anova_all$`Pr(>F)`[1]
    p_env <- anova_env$`Pr(>F)`[1]
    p_spa <- anova_spa$`Pr(>F)`[1]
    p_clim <- NA

    p_pure_env <- anova_pure_env$`Pr(>F)`[1]
    p_pure_spa <- anova_pure_spa$`Pr(>F)`[1]
    p_pure_clim <- NA

    F_all <- anova_all$`F`[1]
    F_env <- anova_env$`F`[1]
    F_spa <- anova_spa$`F`[1]
    F_clim <- NA

    F_pure_env <- anova_pure_env$`F`[1]
    F_pure_spa <- anova_pure_spa$`F`[1]
    F_pure_clim <- NA


    Df_all <- anova_all$`Df`[1]
    Df_env <- anova_env$`Df`[1]
    Df_spa <- anova_spa$`Df`[1]
    Df_clim <- NA

    Df_pure_env <- anova_pure_env$`Df`[1]
    Df_pure_spa <- anova_pure_spa$`Df`[1]
    Df_pure_clim <- NA

    r2_env <- comp.varpart$part$fract$Adj.R.square[1]
    r2_clim <- NA
    r2_spa <- comp.varpart$part$fract$Adj.R.square[2]
    r2_all <- comp.varpart$part$fract$Adj.R.square[3]
    r2_pure_env <- comp.varpart$part$indfract$Adj.R.square[1]
    r2_pure_clim <- NA
    r2_pure_spa <- comp.varpart$part$indfract$Adj.R.square[3]
    r2_share_env_spa <- comp.varpart$part$indfract$Adj.R.square[2]
    r2_share_clim_spa <- NA
    r2_share_clim_env <- NA
    r2_share_clim_env_spa <- NA
    r2_resid <- comp.varpart$part$indfract$Adj.R.square[4]


##################################################

    if(allow_negative_r2 == FALSE){

      if(r2_spa < 0){
        #If negative spatial fraction
        r2_pure_env <- r2_env
        p_pure_env <- p_env
        Df_pure_env <- Df_env
        F_pure_env <- F_env

        r2_all <- r2_env
        p_all <- p_env
        Df_all <- Df_env
        F_all <- F_env

        r2_spa <- 0
        p_spa <-NA
        Df_spa <- NA
        F_spa <-NA

        r2_pure_spa <- 0
        p_pure_spa <-NA
        Df_pure_spa <- NA
        F_pure_spa <-NA

        r2_share_env_spa <- 0

        r2_resid <- 1 - r2_all
      }

      if(r2_env < 0){
        #If negative environmental fraction
        r2_pure_spa <- r2_spa
        p_pure_spa <- p_spa
        Df_pure_spa <- Df_spa
        F_pure_spa <- F_spa

        r2_all <- r2_spa
        p_all <- p_spa
        Df_all <- Df_spa
        F_all <- F_spa

        r2_env <- 0
        p_env <-NA
        Df_env <- NA
        F_env <-NA

        r2_pure_env <- 0
        p_pure_env <-NA
        Df_pure_env <- NA
        F_pure_env <-NA

        r2_share_env_spa <- 0

        r2_resid <- 1 - r2_all
      }

      if(r2_share_env_spa < 0) {
        #If negative environmental space
        r2_pure_spa <- r2_spa
        p_pure_spa <- p_spa
        Df_pure_spa <- Df_spa
        F_pure_spa <- F_spa

        r2_all <- r2_spa + r2_spa
        p_all <- NA
        Df_all <- NA
        F_all <- NA

        r2_pure_env <- r2_env
        p_pure_env <- p_env
        Df_pure_env <- Df_env
        F_pure_env <- F_env

        r2_share_env_spa <- 0

        r2_resid <- 1 - r2_all

      }

      if(r2_spa >= 0 & r2_pure_spa < 0 & r2_env >=0 & r2_pure_env < 0){
        #If negative spatial fraction
        r2_pure_env <- 0
        p_pure_env <- NA
        Df_pure_env <- NA
        F_pure_env <- NA

        r2_pure_spa <- 0
        p_pure_spa <-NA
        Df_pure_spa <- NA
        F_pure_spa <-NA

        if(r2_spa > r2_env){

          r2_pure_spa <- r2_spa - r2_env
          p_pure_spa <-NA
          Df_pure_spa <- NA
          F_pure_spa <-NA

          r2_share_env_spa <- r2_spa - r2_pure_spa

          r2_all <- r2_spa
          p_all <- p_spa
          Df_all <- Df_spa
          F_all <- F_spa


        }else{
          r2_pure_env <- r2_env - r2_spa
          p_pure_env <-NA
          Df_pure_env <- NA
          F_pure_env <-NA

          r2_share_env_spa <- r2_env - r2_pure_env

          r2_all <- r2_env
          p_all <- p_env
          Df_all <- Df_env
          F_all <- F_env
        }

        r2_resid <- 1 - r2_all
      }


    }

    Env <- c(Adj_R2 = r2_env,
             Df = Df_env,
             `F` = F_env,
             p = p_env)

    Clim <- c(Adj_R2 = r2_clim,
              Df = Df_clim,
              `F` = F_clim,
              p = p_clim)

    Spa <- c(Adj_R2 = r2_spa,
             Df = Df_spa,
             `F` = F_spa,
             p = p_spa)

    All <- c(Adj_R2 = r2_all,
             Df = Df_all,
             `F` = F_all,
             p = p_all)

    Pure_Env <- c(Adj_R2 = r2_pure_env,
                  Df = Df_pure_env,
                  `F` = F_pure_env,
                  p = p_pure_env)

    Pure_Clim <- c(Adj_R2 = r2_pure_clim,
                   Df = Df_pure_clim,
                   `F` = F_pure_clim,
                   p = p_pure_clim)

    Pure_Spa <- c(Adj_R2 = r2_pure_spa,
                  Df = Df_pure_spa,
                  `F` = F_pure_spa,
                  p = p_pure_spa)

    Env_Spa <- c(Adj_R2 = r2_share_env_spa, NA, NA, NA)
    Env_Clim <- c(Adj_R2 = r2_share_clim_env, NA, NA, NA)
    Spa_Clim <- c(Adj_R2 = r2_share_clim_spa, NA, NA, NA)
    Spa_Clim_Env <- c(Adj_R2 = r2_share_clim_env_spa, NA, NA, NA)

    Resid <- c(Adj_R2 = r2_resid, NA, NA, NA)



    result <- rbind(All, Env, Clim, Spa, Pure_Env, Pure_Clim, Pure_Spa, Env_Spa, Env_Clim, Spa_Clim, Spa_Clim_Env, Resid)
    if(isTRUE(percent_r2)){
      result[,1] <- round(result[,1]*100,1)
    }else{result[,1] <- round(result[,1],round)}

    return(result)

  }





  #With Environment, Climate and Space
  if((is.null(clim) == FALSE | nrow(clim) == 0) & (is.null(spa) == FALSE | nrow(spa) == 0)){



    comp.varpart <- varpart(Y, env, clim,spa)

    independent.variables <- data.frame(env,clim,spa)
    spa_clim <- data.frame(clim,spa)
    env_clim <- data.frame(clim,env)
    env_spa <- data.frame(spa,env)

    anova_pure_env <- anova.cca(rda(Y, env, spa_clim),permutations=9999)
    anova_pure_spa <- anova.cca(rda(Y, spa, env_clim),permutations=9999)
    anova_pure_clim <- anova.cca(rda(Y, clim, env_spa),permutations=9999)

    anova_env <- anova.cca(rda(Y, env),permutations=9999)
    anova_spa <- anova.cca(rda(Y, spa),permutations=9999)
    anova_clim <- anova.cca(rda(Y, clim),permutations=9999)
    anova_all <- anova.cca(rda(Y, independent.variables),permutations=9999)

    p_all <- anova_all$`Pr(>F)`[1]
    p_env <- anova_env$`Pr(>F)`[1]
    p_spa <- anova_spa$`Pr(>F)`[1]
    p_clim <- anova_clim$`Pr(>F)`[1]

    p_pure_env <- anova_pure_env$`Pr(>F)`[1]
    p_pure_spa <- anova_pure_spa$`Pr(>F)`[1]
    p_pure_clim <- anova_pure_clim$`Pr(>F)`[1]

    F_all <- anova_all$`F`[1]
    F_env <- anova_env$`F`[1]
    F_spa <- anova_spa$`F`[1]
    F_clim <- anova_clim$`F`[1]

    F_pure_env <- anova_pure_env$`F`[1]
    F_pure_spa <- anova_pure_spa$`F`[1]
    F_pure_clim <- anova_pure_clim$`F`[1]


    Df_all <- anova_all$`Df`[1]
    Df_env <- anova_env$`Df`[1]
    Df_spa <- anova_spa$`Df`[1]
    Df_clim <- anova_clim$`Df`[1]

    Df_pure_env <- anova_pure_env$`Df`[1]
    Df_pure_spa <- anova_pure_spa$`Df`[1]
    Df_pure_clim <- anova_pure_clim$`Df`[1]


    r2_env <- comp.varpart$part$fract$Adj.R.square[1]
    r2_clim <- comp.varpart$part$fract$Adj.R.square[2]
    r2_spa <- comp.varpart$part$fract$Adj.R.square[3]
    r2_all <- comp.varpart$part$fract$Adj.R.square[7]
    r2_pure_env <- comp.varpart$part$indfract$Adj.R.square[1]
    r2_pure_clim <- comp.varpart$part$indfract$Adj.R.square[2]
    r2_pure_spa <- comp.varpart$part$indfract$Adj.R.square[3]
    r2_share_env_spa <- comp.varpart$part$indfract$Adj.R.square[6]
    r2_share_clim_spa <- comp.varpart$part$indfract$Adj.R.square[5]
    r2_share_clim_env <- comp.varpart$part$indfract$Adj.R.square[4]
    r2_share_clim_env_spa <- comp.varpart$part$indfract$Adj.R.square[7]
    r2_resid <- comp.varpart$part$indfract$Adj.R.square[8]



    if(allow_negative_r2 == FALSE){

      if(r2_share_clim_env < 0){
        r2_pure_env <- r2_pure_env + r2_share_clim_env
        r2_pure_clim <- r2_pure_clim + r2_share_clim_env
        r2_share_clim_env <- 0
      }
      if(r2_share_clim_spa < 0){
        r2_pure_spa <- r2_pure_spa + r2_share_clim_spa
        r2_pure_clim <- r2_pure_clim + r2_share_clim_spa
        r2_share_clim_spa <- 0
      }
      if(r2_share_env_spa < 0){
        r2_pure_spa <- r2_pure_spa + r2_share_env_spa
        r2_pure_env <- r2_pure_env + r2_share_env_spa
        r2_share_env_spa <- 0
      }
      if(r2_share_clim_env_spa < 0){
        r2_pure_spa <- r2_pure_spa + r2_share_clim_env_spa
        r2_pure_env <- r2_pure_env + r2_share_clim_env_spa
        r2_pure_clim <- r2_pure_clim + r2_share_clim_env_spa
        r2_share_clim_env_spa <- 0
      }
    }

    Env <- c(Adj_R2 = r2_env,
             Df = Df_env,
             `F` = F_env,
             p = p_env)

    Clim <- c(Adj_R2 = r2_clim,
              Df = Df_clim,
              `F` = F_clim,
              p = p_clim)

    Spa <- c(Adj_R2 = r2_spa,
             Df = Df_spa,
             `F` = F_spa,
             p = p_spa)

    All <- c(Adj_R2 = r2_all,
             Df = Df_all,
             `F` = F_all,
             p = p_all)

    Pure_Env <- c(Adj_R2 = r2_pure_env,
                  Df = Df_pure_env,
                  `F` = F_pure_env,
                  p = p_pure_env)

    Pure_Clim <- c(Adj_R2 = r2_pure_clim,
                   Df = Df_pure_clim,
                   `F` = F_pure_clim,
                   p = p_pure_clim)

    Pure_Spa <- c(Adj_R2 = r2_pure_spa,
                  Df = Df_pure_spa,
                  `F` = F_pure_spa,
                  p = p_pure_spa)

    Env_Spa <- c(Adj_R2 = r2_share_env_spa, NA, NA, NA)
    Env_Clim <- c(Adj_R2 = r2_share_clim_env, NA, NA, NA)
    Spa_Clim <- c(Adj_R2 = r2_share_clim_spa, NA, NA, NA)
    Spa_Clim_Env <- c(Adj_R2 = r2_share_clim_env_spa, NA, NA, NA)

    Resid <- c(Adj_R2 = r2_resid, NA, NA, NA)

    result <- rbind(All, Env, Clim, Spa, Pure_Env, Pure_Clim, Pure_Spa, Env_Spa, Env_Clim, Spa_Clim, Spa_Clim_Env, Resid)
    if(isTRUE(percent_r2)){
      result[,1] <- round(result[,1]*100,1)
    }else{result[,1] <- round(result[,1],round)}

    return(result)

  }
}
