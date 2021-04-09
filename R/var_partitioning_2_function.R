#' @title Performing variation partitioning with p values
#'
#' @description This function uses the `varpart()` function from vegan to compute all fractions and their respective p-values in a single table. It differs from the `var_partitioning_1` because it can or cannot omit negative values in the calculations of R2.
#' @param Y Species occurence matrix
#' @param env Matrix of environmental variables
#' @param clim Matrix of climate variables
#' @param spa Matrix of spatial variables
#' @param percent_r2 Logical. Should the R2 values be converted in percentages?
#' @param allow_negative_r2 Logical. Should negative values be used to compute the R2 values?

#'
#' @export
#'


var_partitioning_2 <- function(Y, env, spa, allow_negative_r2 = FALSE, percent_r2 = FALSE){

  comp.varpart <- varpart(Y, env, spa)

  independent.variables <- data.frame(env,spa)
  anova_pure_env <- anova.cca(rda(Y, env, spa),permutations=9999)
  anova_pure_spa <- anova.cca(rda(Y, spa, env),permutations=9999)
  anova_env <- anova.cca(rda(Y, env),permutations=9999)
  anova_spa <- anova.cca(rda(Y, spa),permutations=9999)
  anova_all <- anova.cca(rda(Y, independent.variables),permutations=9999)

  p_all <- anova_all$`Pr(>F)`[1]
  p_env <- anova_env$`Pr(>F)`[1]
  p_spa <- anova_spa$`Pr(>F)`[1]
  p_pure_env <- anova_pure_env$`Pr(>F)`[1]
  p_pure_spa <- anova_pure_spa$`Pr(>F)`[1]

  F_all <- anova_all$`F`[1]
  F_env <- anova_env$`F`[1]
  F_spa <- anova_spa$`F`[1]
  F_pure_env <- anova_pure_env$`F`[1]
  F_pure_spa <- anova_pure_spa$`F`[1]

  Df_all <- anova_all$`Df`[1]
  Df_env <- anova_env$`Df`[1]
  Df_spa <- anova_spa$`Df`[1]
  Df_pure_env <- anova_pure_env$`Df`[1]
  Df_pure_spa <- anova_pure_spa$`Df`[1]


  r2_env <- comp.varpart$part$fract$Adj.R.squared[1]
  r2_spa <- comp.varpart$part$fract$Adj.R.squared[2]
  r2_all <- comp.varpart$part$fract$Adj.R.squared[3]
  r2_pure_env <- comp.varpart$part$indfract$Adj.R.squared[1]
  r2_pure_spa <- comp.varpart$part$indfract$Adj.R.squared[3]
  r2_env_spatially_structured <- comp.varpart$part$indfract$Adj.R.squared[2]
  r2_resid <- comp.varpart$part$indfract$Adj.R.squared[4]

  Env <- c(Adj_R2 = r2_env,
           Df = Df_env,
           `F` = F_env,
           p = p_env)

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

  Pure_Spa <- c(Adj_R2 = r2_pure_spa,
                Df = Df_pure_spa,
                `F` = F_pure_spa,
                p = p_pure_spa)

  Env_Spa <- c(Adj_R2 = r2_env_spatially_structured, NA, NA, NA)

  Resid <- c(Adj_R2 = r2_resid, NA, NA, NA)

  if(isFALSE(allow_negative_r2)){
    if(r2_pure_spa < 0 | r2_spa < 0){
      Resid <- c(Adj_R2 = 1-r2_env, NA, NA, NA)
      result <- rbind(Pure_Env = Env, Resid)
    }else{

      if(r2_pure_env < 0 | r2_env < 0){
        Resid <- c(Adj_R2 = 1-r2_spa, NA, NA, NA)
        result <- rbind(Pure_Spa = Spa, Resid)
      }else{
        if(r2_env_spatially_structured < 0){
          Resid2 <- c(Adj_R2 = 1-(r2_env + r2_spa), NA, NA, NA)
          All2 <- c(Adj_R2 = (r2_env + r2_spa), NA, NA, NA)
          result <- rbind(All_original = All, All = All2, Pure_Env = Env, Pure_Spa = Spa, Resid_original = Resid, Resid = Resid2)
        }else{
          result <- rbind(All, Env, Spa, Pure_Env, Pure_Spa, Env_Spa, Resid)
        }
      }

    }
  }else{
    result <- rbind(All, Env, Spa, Pure_Env, Pure_Spa, Env_Spa, Resid)
  }

  if(isTRUE(percent_r2)){
    result[,1] <- round(result[,1]*100,1)
  }else{result[,1] <- round(result[,1],3)}

  return(result)

}
