Variation Partitioning
================
Rodolfo Pelinson
07/04/2021

First we need preparare all data matrices from the main dataset. These
are prepared sourcing the “Loading\_data.R” file in the Auxiliary
Scripts folder.

``` r
library(AtlanticForestMetacommunity)
source("Loading_data.R")
```

    ## Error in get(genname, envir = envir) : 
    ##   objeto 'testthat_print' não encontrado

The used packages to run this analysis are:

`vegan` version 2.5-6  
`ade4` version 1.7-15  
`adespatial` version 0.3-8

Some packages used to generate plots and tables

`dplyr` version 0.8.5  
`plotrix` version 3.7-8

       

## Reducing Multicolinearity

To remove multicolinear variables from our Environmental Matrix I built
a function that removes variables with the variance inflation factor
(VIF) higher than 3. The function `VIF_selection()` simply uses the
`vif.cca()` from the `vegan` package.

``` r
Broad_env_VIF <- VIF_selection(Broad_pa,  Broad_env_st[,-7])
Broad_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover            nvt 
    ##       1.237565       1.118756       1.340180       1.488527       1.617198 
    ## dist_to_forest 
    ##       1.319102

``` r
SSF_env_VIF <- VIF_selection(SSF_pa,  SSF_env_st)
SSF_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover            nvt 
    ##       1.287250       1.287689       1.443912       1.293051       1.492950 
    ## dist_to_forest 
    ##       1.136891

``` r
ST_env_VIF <- VIF_selection(ST_pa,  ST_env_st)
ST_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth            nvt dist_to_forest 
    ##      27.479904      11.917983       5.386203      20.148713       1.807948 
    ## 
    ## [[2]]
    ##           area          depth            nvt dist_to_forest 
    ##       3.533032       2.360081       2.639484       1.668725 
    ## 
    ## [[3]]
    ##          depth            nvt dist_to_forest 
    ##       1.444899       1.301229       1.367747

``` r
IC_env_VIF <- VIF_selection(IC_pa,  IC_env_st)
IC_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth            nvt dist_to_forest 
    ##       5.342364       5.518649       3.465826       3.219143       3.331453 
    ## 
    ## [[2]]
    ##    hydroperiod          depth            nvt dist_to_forest 
    ##       2.895120       1.516899       1.606162       2.715728

``` r
NI_env_VIF <- VIF_selection(NI_pa,  NI_env_st)
NI_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth            nvt dist_to_forest 
    ##       1.502173       4.868245       4.697667       1.834039       2.256836 
    ## 
    ## [[2]]
    ##    hydroperiod          depth            nvt dist_to_forest 
    ##       1.487856       1.157613       1.493676       1.441587

``` r
MD_env_VIF <- VIF_selection(MD_pa,  MD_env_st)
MD_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover            nvt 
    ##       6.123821      15.263901       3.075865      36.928117       3.436439 
    ## dist_to_forest 
    ##      43.153516 
    ## 
    ## [[2]]
    ##  hydroperiod         area        depth canopy_cover          nvt 
    ##     3.843659     5.706331     2.499002     1.637214     3.434896 
    ## 
    ## [[3]]
    ##  hydroperiod        depth canopy_cover          nvt 
    ##     1.770317     1.879809     1.410612     1.508743

``` r
JA_env_VIF <- VIF_selection(JA_pa,  JA_env_st)
JA_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover            nvt 
    ##      14.962311       3.081842       5.988699       7.691347       3.453394 
    ## dist_to_forest 
    ##      12.622399 
    ## 
    ## [[2]]
    ##           area          depth   canopy_cover            nvt dist_to_forest 
    ##       1.961897       4.322304       6.871158       3.356278       7.533870 
    ## 
    ## [[3]]
    ##         area        depth canopy_cover          nvt 
    ##     1.450483     2.089961     2.081341     2.864143

``` r
DRF_env_VIF <- VIF_selection(DRF_pa,  DRF_env_st)
DRF_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover            nvt 
    ##       1.338500       1.296349       1.356465       1.272615       1.048673 
    ## dist_to_forest 
    ##       1.368314

``` r
UBA_env_VIF <- VIF_selection(UBA_pa,  UBA_env_st)
UBA_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover dist_to_forest 
    ##       1.436660       7.001901       1.750901       1.459906       6.699443 
    ## 
    ## [[2]]
    ##    hydroperiod          depth   canopy_cover dist_to_forest 
    ##       1.341421       1.599294       1.415973       1.643447

``` r
BER_env_VIF <- VIF_selection(BER_pa,  BER_env_st)
BER_env_VIF$VIFs
```

    ## [[1]]
    ##    hydroperiod           area          depth   canopy_cover dist_to_forest 
    ##       3.117432       8.428997       2.304574       1.663244       7.972957 
    ## 
    ## [[2]]
    ##    hydroperiod          depth   canopy_cover dist_to_forest 
    ##       2.192425       2.253493       1.561949       4.070420 
    ## 
    ## [[3]]
    ##  hydroperiod        depth canopy_cover 
    ##     1.541888     1.154773     1.384145

``` r
ITA_env_VIF <- VIF_selection(ITA_pa,  ITA_env_st)
ITA_env_VIF$VIFs
```

    ## [[1]]
    ##  hydroperiod         area        depth canopy_cover          nvt 
    ##     2.342853     1.671319     1.996469     1.828709     1.458903

   

Same thing for climate variables

``` r
Broad_clim_VIF <- VIF_selection(Broad_pa,  Broad_clim_st[,-6])
Broad_clim_VIF$VIFs
```

    ## [[1]]
    ## temp_Season  range_temp  total_prec prec_season 
    ##    16.76472    11.05404    37.92818    37.69109 
    ## 
    ## [[2]]
    ## temp_Season  range_temp prec_season 
    ##    3.035253    2.145741    4.869951 
    ## 
    ## [[3]]
    ## temp_Season  range_temp 
    ##    1.167625    1.167625

``` r
SSF_clim_VIF <- VIF_selection(SSF_pa,  SSF_clim_st)
SSF_clim_VIF$VIFs
```

    ## [[1]]
    ## temp_Season  range_temp  total_prec prec_season 
    ##   45.128386    5.922860    3.172604   52.669397 
    ## 
    ## [[2]]
    ## temp_Season  range_temp  total_prec 
    ##    2.035689    4.718471    3.046979 
    ## 
    ## [[3]]
    ## temp_Season  total_prec 
    ##     1.15459     1.15459

``` r
DRF_clim_VIF <- VIF_selection(DRF_pa,  DRF_clim_st)
DRF_clim_VIF$VIFs
```

    ## [[1]]
    ## temp_Season  range_temp  total_prec prec_season 
    ##    6.512557    3.413672   11.946292    6.515921 
    ## 
    ## [[2]]
    ## temp_Season  range_temp prec_season 
    ##    1.881731    2.874642    2.327331

   

## Spatial Variables

Constructing matrices of spatial filters

``` r
set.seed(5)

candidates_Broad <- listw.candidates(Broad_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup"), y_fdown = 5, y_fup = 0.5)
```

    ## 
    ##      PLEASE NOTE:  The components "delsgs" and "summary" of the
    ##  object returned by deldir() are now DATA FRAMES rather than
    ##  matrices (as they were prior to release 0.0-18).
    ##  See help("deldir").
    ##  
    ##      PLEASE NOTE: The process that deldir() uses for determining
    ##  duplicated points has changed from that used in version
    ##  0.0-9 of this package (and previously). See help("deldir").

``` r
Broad_MEM <- listw.select(Broad_pa, candidates = candidates_Broad, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 11 is 0.054295 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 18 is 0.059994 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 17 is 0.070493 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 20 is 0.056794 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 21 is 0.052895 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 23 is 0.065293 (> 0.050000)

``` r
Broad_MEM_FS <- Broad_MEM$best$MEM.select

adegraphics::s.label(Broad_coord, nb = candidates_Broad[[Broad_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("Broad - ",names(Broad_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-1.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
################################################################################


candidates_DRF <- listw.candidates(DRF_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup"), y_fdown = 5, y_fup = 0.25)

DRF_MEM <- listw.select(DRF_pa, candidates = candidates_DRF, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.057465 with 2 variables (> 0.050471)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.086391 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.106289 (> 0.050000)

``` r
DRF_MEM_FS <- DRF_MEM$best$MEM.select


adegraphics::s.label(DRF_coord, nb = candidates_DRF[[DRF_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("DRF - ",names(DRF_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-2.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
################################################################################


candidates_SSF <- listw.candidates(SSF_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup"), y_fdown = 5, y_fup = 0.25)

SSF_MEM <- listw.select(SSF_pa, candidates = candidates_SSF, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.079183 with 2 variables (> 0.079183)
    ## Procedure stopped (alpha criteria): pvalue for variable 7 is 0.083592 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.304147 with 9 variables (> 0.302779)
    ## Procedure stopped (alpha criteria): pvalue for variable 7 is 0.077392 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.069193 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 9 is 0.071693 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 10 is 0.057494 (> 0.050000)

``` r
SSF_MEM_FS <- SSF_MEM$best$MEM.select


adegraphics::s.label(SSF_coord, nb = candidates_SSF[[SSF_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("SSF - ",names(SSF_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-3.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
################################################################################


#We restricted our options for optimization because of the low statistical power (low replication) that we had at our small spatial scale. We only used linear weights as we do not believe that a exponential decay would make sense at this scale (even relative large distances between sites are actually small). Also We restricted our graph conectivity matrix to only three that yields relatively different scenarios of conectance. Those are the Delauney triangulation, the Relative neighbour and the PCNM

################################################################################



candidates_ITA <- listw.candidates(ITA_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

ITA_MEM <- listw.select(ITA_pa, candidates = candidates_ITA, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)


ITA_MEM_FS <- dbmem(ITA_coord, MEM.autocor = c("positive"), silent = TRUE)
################################################################################


candidates_BER <- listw.candidates(BER_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

BER_MEM <- listw.select(BER_pa, candidates = candidates_BER, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)

BER_MEM_FS <- dbmem(BER_coord, MEM.autocor = c("positive"), silent = TRUE)
################################################################################

candidates_UBA <- listw.candidates(UBA_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))
```

    ## Warning in nb2listw(nb.object, style = style, glist = lapply(nb.dist, f1, : zero
    ## sum general weights

``` r
UBA_MEM <- listw.select(UBA_pa, candidates = candidates_UBA, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.091291 (> 0.050000)

``` r
UBA_MEM_FS <- UBA_MEM$best$MEM.select

adegraphics::s.label(UBA_coord, nb = candidates_UBA[[UBA_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("UBA - ",names(UBA_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-4.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
################################################################################


candidates_ST <- listw.candidates(ST_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))
```

    ## Warning in nb2listw(nb.object, style = style, glist = lapply(nb.dist, f1, : zero
    ## sum general weights

``` r
ST_MEM <- listw.select(ST_pa, candidates = candidates_ST, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)

ST_MEM_FS <- dbmem(ST_coord, MEM.autocor = c("positive"), silent = TRUE)

################################################################################


candidates_IC <- listw.candidates(IC_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

IC_MEM <- listw.select(IC_pa, candidates = candidates_IC, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)

IC_MEM_FS <- dbmem(IC_coord, MEM.autocor = c("positive"), silent = TRUE)

################################################################################

candidates_NI <- listw.candidates(NI_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

NI_MEM <- listw.select(NI_pa, candidates = candidates_NI, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)


adegraphics::s.label(NI_coord, nb = candidates_NI[[NI_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("NI - ",names(NI_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-5.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
NI_MEM_FS <- dbmem(NI_coord, MEM.autocor = c("positive"), silent = TRUE)


################################################################################

candidates_MD <- listw.candidates(MD_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

MD_MEM <- listw.select(MD_pa, candidates = candidates_MD, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)

MD_MEM_FS <- dbmem(MD_coord, MEM.autocor = c("positive"), silent = TRUE)

################################################################################


candidates_JA <- listw.candidates(JA_coord, style = "B", nb = c("del", "rel", "pcnm"),
                                   weights = c("flin"))

JA_MEM <- listw.select(JA_pa, candidates = candidates_JA, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = TRUE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.082792 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.103890 (> 0.050000)

``` r
JA_MEM_FS <- JA_MEM$best$MEM.select


adegraphics::s.label(JA_coord, nb = candidates_JA[[JA_MEM$best.id]],
                     pnb.edge.col = "red", main = paste("JA - ",names(JA_MEM$best.id)), plot = TRUE, labels = NULL)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-6-6.png" width="400" height="400" style="display: block; margin: auto;" />

``` r
################################################################################
```

``` r
###PLOTING WEIGHTS ####
nb <- chooseCN(coordinates(UBA_coord), type = 1, plot.nb = FALSE)

#Delaunay triangulation (type 1)
#Gabriel graph (type 2)
#Relative neighbours (type 3)
#Minimum spanning tree (type 4)
#Neighbourhood by distance (type 5)
#K nearests neighbours (type 6)
#Inverse distances (type 7)

dist_UBA <- dist(UBA_coord)

distnb <- spdep::nbdists(nb, as.matrix(UBA_coord))

a_fdown <- 5
a_f_up <- 0.5

fdist <- lapply(distnb, function(x) 1 - x/max(dist_UBA)) #linear
fdist <- lapply(distnb, function(x) 1 - (x/max(dist_UBA))^a_fdown) #fdown
fdist <- lapply(distnb, function(x) 1 / x^a_f_up) #fup
    
lw <- nb2listw(nb, style = 'B', zero.policy = TRUE, glist = fdist) 

m1 <- as.matrix(dist_UBA)
m2 <- spdep::listw2mat(lw)

plot(m1[!diag(ncol(m1))], m2[!diag(ncol(m2))], pch = 20, xlab = "distance", ylab = "spatial weights")
################
```

       

## Forward Selection

Similarly to removing variables with VIF higher than 3, I built a
function to select the most important variables, when they could
significantly explain species occurrences. This function mostly rely on
the function `fs()` from package `adespatial`. Forward selection is only
performed when the whole predictor matrix can significantly explain (p
\< 0.05) the variation in species occurrences.

### Broad Spatial Extent

``` r
set.seed(5)
Broad_env_Forward <- forward_selection(Broad_pa, Broad_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Testing variable 4
    ## Testing variable 5
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.176872 with 5 variables (> 0.175018)

``` r
Broad_env_Forward$forward_results
```

    ##        variables order         R2     R2Cum  AdjR2Cum         F pvalue
    ## 1            nvt     5 0.12280383 0.1228038 0.1134720 13.159611 0.0001
    ## 2 dist_to_forest     6 0.03553797 0.1583418 0.1402416  3.926809 0.0001
    ## 3   canopy_cover     4 0.02254986 0.1808917 0.1541816  2.532738 0.0018
    ## 4    hydroperiod     1 0.02048938 0.2013810 0.1662769  2.334697 0.0066

``` r
Broad_env_FS <- Broad_env_Forward$selected_variables

Broad_clim_Forward <- forward_selection(Broad_pa, Broad_clim_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.211063 with 2 variables (> 0.211063)

``` r
Broad_clim_Forward$forward_results
```

    ##    variables order        R2     R2Cum AdjR2Cum        F pvalue
    ## 1 range_temp     2 0.1858759 0.1858759 0.177215 21.46152  1e-04

``` r
Broad_clim_FS <- Broad_clim_Forward$selected_variables

#_______________________________________________________________
SSF_env_Forward <- forward_selection(SSF_pa, SSF_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Testing variable 4
    ## Testing variable 5
    ## Procedure stopped (alpha criteria): pvalue for variable 5 is 0.122000 (> 0.050000)

``` r
SSF_env_Forward$forward_results
```

    ##        variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1 dist_to_forest     6 0.06279225 0.06279225 0.04149207 2.947968 0.0013
    ## 2    hydroperiod     1 0.05740783 0.12020008 0.07927915 2.805793 0.0017
    ## 3          depth     3 0.05098863 0.17118871 0.11198790 2.583848 0.0043
    ## 4   canopy_cover     4 0.04138347 0.21257218 0.13574995 2.154765 0.0132

``` r
SSF_env_FS <- SSF_env_Forward$selected_variables

SSF_clim_Forward <- forward_selection(SSF_pa, SSF_clim_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2

``` r
SSF_clim_Forward$forward_results
```

    ##     variables order        R2     R2Cum   AdjR2Cum        F pvalue
    ## 1  total_prec     2 0.1035181 0.1035181 0.08314354 5.080747  2e-04
    ## 2 temp_Season     1 0.1046860 0.2082041 0.17137639 5.685174  1e-04

``` r
SSF_clim_FS <- SSF_clim_Forward$selected_variables

#_______________________________________________________________
ST_env_Forward <- forward_selection(ST_pa, ST_env_VIF$variables)
```

    ## Forward selection NOT performed. p > 0.05

``` r
ST_env_Forward$forward_results
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: rda(X = New_Y, Y = New_X)
    ##          Df Variance     F Pr(>F)
    ## Model     3   1.5233 1.632  0.103
    ## Residual  4   1.2446

``` r
ST_env_FS <- ST_env_Forward$selected_variables

IC_env_Forward <- forward_selection(IC_pa, IC_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.414300 (> 0.050000)

``` r
IC_env_Forward$forward_results
```

    ##        variables order        R2     R2Cum  AdjR2Cum        F pvalue
    ## 1 dist_to_forest     4 0.2762401 0.2762401 0.2038641 3.816737 0.0005
    ## 2    hydroperiod     1 0.1448824 0.4211225 0.2924830 2.252534 0.0105

``` r
IC_env_FS <- IC_env_Forward$selected_variables

NI_env_Forward <- forward_selection(NI_pa, NI_env_VIF$variables)
```

    ## Forward selection NOT performed. p > 0.05

``` r
NI_env_Forward$forward_results
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: rda(X = New_Y, Y = New_X)
    ##          Df Variance      F Pr(>F)
    ## Model     4   2.2203 1.3581 0.1508
    ## Residual  3   1.2261

``` r
NI_env_FS <- NI_env_Forward$selected_variables

MD_env_Forward <- forward_selection(MD_pa, MD_env_VIF$variables)
```

    ## Forward selection NOT performed. p > 0.05

``` r
MD_env_Forward$forward_results
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: rda(X = New_Y, Y = New_X)
    ##          Df Variance      F Pr(>F)
    ## Model     4   1.9377 1.2819 0.2664
    ## Residual  3   1.1337

``` r
MD_env_FS <- MD_env_Forward$selected_variables

JA_env_Forward <- forward_selection(JA_pa, JA_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.058600 (> 0.050000)

``` r
JA_env_Forward$forward_results
```

    ##   variables order        R2     R2Cum  AdjR2Cum        F pvalue
    ## 1     depth     2 0.2978698 0.2978698 0.2101035 3.393898 0.0016

``` r
JA_env_FS <- JA_env_Forward$selected_variables

#______________________________________________________________
DRF_env_Forward <- forward_selection(DRF_pa, DRF_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Testing variable 3
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.082422 with 3 variables (> 0.081064)

``` r
DRF_env_Forward$forward_results
```

    ##      variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1 canopy_cover     4 0.06947213 0.06947213 0.05008614 3.583624 0.0006
    ## 2        depth     3 0.04152619 0.11099833 0.07316847 2.195419 0.0191

``` r
DRF_env_FS <- DRF_env_Forward$selected_variables

DRF_clim_Forward <- forward_selection(DRF_pa, DRF_clim_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.059286 with 2 variables (> 0.054878)

``` r
DRF_clim_Forward$forward_results
```

    ##    variables order         R2      R2Cum   AdjR2Cum        F pvalue
    ## 1 range_temp     2 0.06813262 0.06813262 0.04871872 3.509476  7e-04

``` r
DRF_clim_FS <- DRF_clim_Forward$selected_variables

#______________________________________________________________
UBA_env_Forward <- forward_selection(UBA_pa, UBA_env_VIF$variables)
```

    ## Testing variable 1
    ## Testing variable 2
    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.279600 (> 0.050000)

``` r
UBA_env_Forward$forward_results
```

    ##      variables order        R2     R2Cum   AdjR2Cum        F pvalue
    ## 1 canopy_cover     3 0.1360444 0.1360444 0.09490361 3.306803 0.0015

``` r
UBA_env_FS <- UBA_env_Forward$selected_variables

BER_env_Forward <- forward_selection(BER_pa, BER_env_VIF$variables)
```

    ## Forward selection NOT performed. p > 0.05

``` r
BER_env_Forward$forward_results
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: rda(X = New_Y, Y = New_X)
    ##          Df Variance      F Pr(>F)
    ## Model     3  0.61185 0.8438 0.6938
    ## Residual  8  1.93360

``` r
BER_env_FS <- BER_env_Forward$selected_variables

ITA_env_Forward <- forward_selection(ITA_pa, ITA_env_VIF$variables)
```

    ## Forward selection NOT performed. p > 0.05

``` r
ITA_env_Forward$forward_results
```

    ## Permutation test for rda under reduced model
    ## Permutation: free
    ## Number of permutations: 9999
    ## 
    ## Model: rda(X = New_Y, Y = New_X)
    ##          Df Variance      F Pr(>F)
    ## Model     5   1.0943 1.1421 0.3307
    ## Residual  9   1.7247

``` r
ITA_env_FS <- ITA_env_Forward$selected_variables
```

We can plot the important spatial variables to better understand what
spatial patterns they describe.

We are only going to plot the four most important ones

``` r
par(mfrow = c(2,2))

sr_value(Broad_coord, data.frame(Broad_MEM$best$MEM.select)[,1], ylim = c(-24.49270,-20.17833), xlim = c(-52.64536, -44.51492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("LARGE SCALE", colnames(data.frame(Broad_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)

sr_value(Broad_coord, data.frame(Broad_MEM$best$MEM.select)[,2], ylim = c(-24.49270,-20.17833), xlim = c(-52.64536, -44.51492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("LARGE SCALE", colnames(data.frame(Broad_MEM$best$MEM.select))[2]), line = 3, outer = F, adj = 1)

sr_value(Broad_coord, data.frame(Broad_MEM$best$MEM.select)[,3], ylim = c(-24.49270,-20.17833), xlim = c(-52.64536, -44.51492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("LARGE SCALE", colnames(data.frame(Broad_MEM$best$MEM.select))[3]), line = 3, outer = F, adj = 1)

sr_value(Broad_coord, data.frame(Broad_MEM$best$MEM.select)[,4], ylim = c(-24.49270,-20.17833), xlim = c(-52.64536, -44.51492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("LARGE SCALE", colnames(data.frame(Broad_MEM$best$MEM.select))[4]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-9-1.png" width="800" height="800" style="display: block; margin: auto;" />

       

Intermediate Scale DRF.

``` r
par(mfrow = c(2,2))

sr_value(DRF_coord, data.frame(DRF_MEM$best$MEM.select)[,1], ylim = c(-24.59270,-23.33123), xlim = c(-47.61181, -44.61492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("DRF", colnames(data.frame(DRF_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)

sr_value(DRF_coord, data.frame(DRF_MEM$best$MEM.select)[,2], ylim = c(-24.59270,-23.33123), xlim = c(-47.61181, -44.61492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("DRF", colnames(data.frame(DRF_MEM$best$MEM.select))[2]), line = 3, outer = F, adj = 1)

sr_value(DRF_coord, data.frame(DRF_MEM$best$MEM.select)[,3], ylim = c(-24.59270,-23.33123), xlim = c(-47.61181, -44.61492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("DRF", colnames(data.frame(DRF_MEM$best$MEM.select))[3]), line = 3, outer = F, adj = 1)

sr_value(DRF_coord, data.frame(DRF_MEM$best$MEM.select)[,4], ylim = c(-24.59270,-23.33123), xlim = c(-47.61181, -44.61492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("DRF", colnames(data.frame(DRF_MEM$best$MEM.select))[4]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-10-1.png" width="800" height="800" style="display: block; margin: auto;" />
       

Intermediate Scale SSF.

``` r
par(mfrow = c(2,2))

sr_value(SSF_coord, data.frame(SSF_MEM$best$MEM.select)[,1], ylim = c(-22.61958,-20.17833), xlim = c(-52.54536, -47.52589), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("SSF", colnames(data.frame(SSF_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)

sr_value(SSF_coord, data.frame(SSF_MEM$best$MEM.select)[,2], ylim = c(-22.61958,-20.17833), xlim = c(-52.54536, -47.52589), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("SSF", colnames(data.frame(SSF_MEM$best$MEM.select))[2]), line = 3, outer = F, adj = 1)

sr_value(SSF_coord, data.frame(SSF_MEM$best$MEM.select)[,3], ylim = c(-22.61958,-20.17833), xlim = c(-52.54536, -47.52589), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("SSF", colnames(data.frame(SSF_MEM$best$MEM.select))[3]), line = 3, outer = F, adj = 1)

sr_value(SSF_coord, data.frame(SSF_MEM$best$MEM.select)[,4], ylim = c(-22.61958,-20.17833), xlim = c(-52.54536, -47.52589), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("SSF", colnames(data.frame(SSF_MEM$best$MEM.select))[4]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-11-1.png" width="800" height="800" style="display: block; margin: auto;" />

       

Small Extent DRF - Ubatuba

``` r
par(mfrow = c(1,2))

sr_value(UBA_coord, data.frame(UBA_MEM$best$MEM.select)[,1], ylim = c(-23.37694,-23.33123), xlim = c(-44.95004, -44.80492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("UBA", colnames(data.frame(UBA_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)

sr_value(UBA_coord, data.frame(UBA_MEM$best$MEM.select)[,2], ylim = c(-23.37694,-23.33123), xlim = c(-44.95004, -44.80492), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("UBA", colnames(data.frame(UBA_MEM$best$MEM.select))[2]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-12-1.png" width="800" height="400" style="display: block; margin: auto;" />

   

Small Extent SSF - Nova Itapirema.

``` r
par(mfrow = c(1,1))

sr_value(NI_coord, data.frame(NI_MEM$best$MEM.select)[,1], ylim = c(-21.08111,-21.07333), xlim = c(-49.54072, -49.51689), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("NI", colnames(data.frame(NI_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-13-1.png" width="400" height="400" style="display: block; margin: auto;" />

Small Extent SSF - Jataí.

``` r
par(mfrow = c(1,2))

sr_value(JA_coord, data.frame(JA_MEM$best$MEM.select)[,1], ylim = c(-21.58556, -21.56376), xlim = c(-47.79086, -47.72289), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("JA", colnames(data.frame(JA_MEM$best$MEM.select))[1]), line = 3, outer = F, adj = 1)

sr_value(JA_coord, data.frame(JA_MEM$best$MEM.select)[,2], ylim = c(-21.58556, -21.56376), xlim = c(-47.79086, -47.72289), grid=F, csize = 0.8, clegend = 1, xax = 2, yax = 1, method = "bubble")
title(main = paste("JA", colnames(data.frame(JA_MEM$best$MEM.select))[2]), line = 3, outer = F, adj = 1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-14-1.png" width="800" height="400" style="display: block; margin: auto;" />

       

## Variation Partitioning

#### Large Extent

``` r
#Variation partitioning Broad

Broad_varpart <- var_partitioning(Y = Broad_pa, 
                                       env = Broad_env_FS,
                                       clim = Broad_clim_FS,
                                       spa = data.frame(Broad_MEM$best$MEM.select), percent_r2 = F)

#Testing significande of spatially structured environment
candidates_Broad_env <- listw.candidates(Broad_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)

Broad_MEM_env <- listw.select(Broad_env_FS, candidates = candidates_Broad_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.170783 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.381126 with 7 variables (> 0.377983)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.356546 with 8 variables (> 0.352020)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.448108 with 14 variables (> 0.443978)
    ## Procedure stopped (alpha criteria): pvalue for variable 14 is 0.052995 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 16 is 0.053295 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 13 is 0.054895 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.422198 with 14 variables (> 0.418147)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.440149 with 16 variables (> 0.433487)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.421414 with 16 variables (> 0.408867)

``` r
Broad_env_spa_sig <- envspace.test(spe = Broad_pa, env = Broad_env_FS,   coord = Broad_coord,  
              MEM.spe =  data.frame(Broad_MEM$best$MEM.select), 
              listw.env = candidates_Broad_env[[Broad_MEM_env$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)

#Testing significande of spatially structured climate
Broad_MEM_clim <- listw.select(Broad_clim_FS, candidates = candidates_Broad_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.932942 with 13 variables (> 0.930662)
    ## Procedure stopped (alpha criteria): pvalue for variable 15 is 0.051095 (> 0.050000)
    ## Procedure stopped (R2more criteria): variable 23 explains only 0.000953 of the variance.
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.968621 with 18 variables (> 0.968432)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.962569 with 20 variables (> 0.961720)
    ## Procedure stopped (alpha criteria): pvalue for variable 28 is 0.054795 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.960383 with 20 variables (> 0.958395)
    ## Procedure stopped (alpha criteria): pvalue for variable 19 is 0.051595 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 29 is 0.066693 (> 0.050000)

``` r
Broad_clim_spa_sig <- envspace.test(spe = Broad_pa, env = Broad_clim_FS,   coord = Broad_coord,  
              MEM.spe =  data.frame(Broad_MEM$best$MEM.select), 
              listw.env = candidates_Broad_env[[Broad_MEM_clim$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)



Broad_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.411 27  3.454546 0.0001
    ## Env           0.166  4  5.736676 0.0001
    ## Clim          0.177  1 21.461519 0.0001
    ## Spa           0.369 22  3.527408 0.0001
    ## Pure_Env      0.012  4  1.364727 0.0122
    ## Pure_Clim     0.026  1  3.998770 0.0001
    ## Pure_Spa      0.168 22  2.165478 0.0001
    ## Env_Spa       0.054 NA        NA     NA
    ## Env_Clim      0.004 NA        NA     NA
    ## Spa_Clim      0.051 NA        NA     NA
    ## Spa_Clim_Env  0.096 NA        NA     NA
    ## Resid         0.589 NA        NA     NA

``` r
Broad_env_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.1501444 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 9.999e-05 
    ## Alternative hypothesis: greater 
    ## 
    ##     Std.Obs Expectation    Variance 
    ## 4.848439708 0.034941702 0.000564575

``` r
Broad_clim_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.1474176 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.00029997 
    ## Alternative hypothesis: greater 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 5.2302997382 0.0212520490 0.0005818731

   

#### Intermediate Extent

``` r
DRF_varpart <- var_partitioning(Y = DRF_pa, 
                                      env = DRF_env_FS,
                                      clim = DRF_clim_FS,
                                      spa = DRF_MEM_FS, percent_r2 = F)

#Testing significande of spatially structured environment
candidates_DRF_env <- listw.candidates(DRF_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)

DRF_MEM_env <- listw.select(DRF_env_FS, candidates = candidates_DRF_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.078192 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.056594 (> 0.050000)

``` r
DRF_env_spa_sig <- envspace.test(spe = DRF_pa, env = DRF_env_FS,   coord = DRF_coord,  
              MEM.spe =  data.frame(DRF_MEM$best$MEM.select), 
              listw.env = candidates_DRF_env[[DRF_MEM_env$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)


#Testing significande of spatially structured climate
DRF_MEM_clim <- listw.select(DRF_clim_FS, candidates = candidates_DRF_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.877894 with 5 variables (> 0.863802)
    ## Procedure stopped (alpha criteria): pvalue for variable 5 is 0.051295 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.954497 with 15 variables (> 0.953531)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.949753 with 15 variables (> 0.947409)
    ## Procedure stopped (alpha criteria): pvalue for variable 9 is 0.081492 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 16 is 0.058794 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 11 is 0.069593 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.899555 with 9 variables (> 0.898920)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.879873 with 14 variables (> 0.879685)

``` r
DRF_clim_spa_sig <- envspace.test(spe = DRF_pa, env = DRF_clim_FS,   coord = DRF_coord,  
              MEM.spe =  data.frame(DRF_MEM$best$MEM.select), 
              listw.env = candidates_DRF_env[[DRF_MEM_clim$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)




###########################################################################################
SSF_varpart <- var_partitioning(Y = SSF_pa, 
                                      env = SSF_env_FS,
                                      clim = SSF_clim_FS,
                                      spa = SSF_MEM_FS, percent_r2 = F)

#Testing significande of spatially structured environment
candidates_SSF_env <- listw.candidates(SSF_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)

SSF_MEM_env <- listw.select(SSF_env_FS, candidates = candidates_SSF_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.123588 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 6 is 0.090291 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.358144 with 8 variables (> 0.340049)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.053295 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 10 is 0.078892 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 7 is 0.054995 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 7 is 0.058694 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.054395 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.089191 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 3 is 0.054395 (> 0.050000)

``` r
SSF_env_spa_sig <- envspace.test(spe = SSF_pa, env = SSF_env_FS,   coord = SSF_coord,  
              MEM.spe =  data.frame(SSF_MEM$best$MEM.select), 
              listw.env = candidates_SSF_env[[SSF_MEM_env$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)


#Testing significande of spatially structured climate
SSF_MEM_clim <- listw.select(SSF_clim_FS, candidates = candidates_SSF_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (alpha criteria): pvalue for variable 6 is 0.051495 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 7 is 0.052695 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.966230 with 12 variables (> 0.965678)
    ## Procedure stopped (alpha criteria): pvalue for variable 12 is 0.071293 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 15 is 0.067293 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 14 is 0.066593 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.936295 with 8 variables (> 0.931722)
    ## Procedure stopped (alpha criteria): pvalue for variable 9 is 0.058194 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 8 is 0.054395 (> 0.050000)

``` r
SSF_clim_spa_sig <- envspace.test(spe = SSF_pa, env = SSF_clim_FS,   coord = SSF_coord,  
              MEM.spe =  data.frame(SSF_MEM$best$MEM.select), 
              listw.env = candidates_SSF_env[[SSF_MEM_clim$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)



DRF_varpart
```

    ##              Adj_R2 Df        F      p
    ## All           0.191 10 2.154592 0.0001
    ## Env           0.073  2 2.934146 0.0003
    ## Clim          0.049  1 3.509476 0.0003
    ## Spa           0.165  7 2.378271 0.0001
    ## Pure_Env      0.016  2 1.407495 0.0827
    ## Pure_Clim     0.006  1 1.291486 0.1928
    ## Pure_Spa      0.091  7 1.738194 0.0003
    ## Env_Spa       0.035 NA       NA     NA
    ## Env_Clim      0.004 NA       NA     NA
    ## Spa_Clim      0.021 NA       NA     NA
    ## Spa_Clim_Env  0.018 NA       NA     NA
    ## Resid         0.809 NA       NA     NA

``` r
DRF_env_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.05287454 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.01089891 
    ## Alternative hypothesis: greater 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 2.4182696166 0.0119710565 0.0002860955

``` r
DRF_clim_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.03861442 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.1041896 
    ## Alternative hypothesis: greater 
    ## 
    ##     Std.Obs Expectation    Variance 
    ## 1.289184542 0.023964214 0.000129139

``` r
SSF_varpart
```

    ##              Adj_R2 Df        F      p
    ## All           0.346 14 2.703080 0.0001
    ## Env           0.136  4 2.767066 0.0001
    ## Clim          0.171  2 5.653462 0.0001
    ## Spa           0.294  8 3.344447 0.0001
    ## Pure_Env      0.034  4 1.450372 0.0211
    ## Pure_Clim     0.009  2 1.224655 0.1878
    ## Pure_Spa      0.078  8 1.580929 0.0005
    ## Env_Spa       0.063 NA       NA     NA
    ## Env_Clim      0.010 NA       NA     NA
    ## Spa_Clim      0.124 NA       NA     NA
    ## Spa_Clim_Env  0.029 NA       NA     NA
    ## Resid         0.654 NA       NA     NA

``` r
SSF_env_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.09248927 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.00859914 
    ## Alternative hypothesis: greater 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 2.5157353962 0.0334944218 0.0005499183

``` r
SSF_clim_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.1528604 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.00079992 
    ## Alternative hypothesis: greater 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 4.1851154744 0.0458508808 0.0006537776

   

#### Small Extent

Small Extent allowing negative fractions

``` r
ST_varpart <- var_partitioning(ST_pa, ST_env_FS, spa = ST_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
IC_varpart <- var_partitioning(IC_pa, IC_env_FS, spa = IC_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
NI_varpart <- var_partitioning(NI_pa, NI_env_FS, spa = NI_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
#Testing significande of spatially structured environment
candidates_NI_env <- listw.candidates(NI_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)

NI_MEM_env <- listw.select(NI_env_FS, candidates = candidates_NI_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.1, p.adjust = FALSE, verbose = FALSE) 
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.175515 with 1 variables (> 0.175515)

``` r
# We are not adjusting p-values here because we assume that if there was a significant spatial and environmental strucutre in species distributions, we want to find spatial patterns in the environment

#There was no significant spatial structure in the environmental variables
NI_env_spa_sig <- envspace.test(spe = NI_pa, env = NI_env_FS,   coord = NI_coord,  
              MEM.spe =  NI_MEM_FS, 
              listw.env = candidates_NI_env$DBEM_PCNM, MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)
```

    ## No significant relation between 'spe' and 'env'; the test was not performed

``` r
MD_varpart <- var_partitioning(MD_pa, MD_env_FS, spa = MD_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
JA_varpart <- var_partitioning(JA_pa, JA_env_FS, spa = JA_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
#Testing significande of spatially structured environment
candidates_JA_env <- listw.candidates(JA_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)

JA_MEM_env <- listw.select(JA_env_FS, candidates = candidates_JA_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.1, p.adjust = FALSE, verbose = FALSE) # We are not adjusting p-values here because we assume that if there was a significant spatial and environmental strucutre in species distributions, we want to find spatial patterns in the environment
```

    ## Procedure stopped (alpha criteria): pvalue for variable 1 is 0.058094 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 1 is 0.057494 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.892776 with 2 variables (> 0.874947)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.841344 with 2 variables (> 0.808686)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.852574 with 3 variables (> 0.825464)
    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.055394 (> 0.050000)
    ## Procedure stopped (alpha criteria): pvalue for variable 2 is 0.133287 (> 0.050000)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.855931 with 2 variables (> 0.843916)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.859910 with 3 variables (> 0.837980)

``` r
#There was no significant spatial structure in the environmental variables
JA_env_spa_sig <- envspace.test(spe = JA_pa, env = JA_env_FS,   coord = JA_coord,  
              MEM.spe =  data.frame(JA_MEM$best$MEM.select), 
              listw.env = candidates_JA_env[[JA_MEM_env$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)




UBA_varpart <- var_partitioning(UBA_pa, UBA_env_FS, spa = UBA_MEM_FS, allow_negative_r2 = T, percent_r2 = F)

#Testing significande of spatially structured environment
candidates_UBA_env <- listw.candidates(UBA_coord, style = "B", nb = c("del", "gab", "rel", "pcnm"),
                                   weights = c("flin", "fup","fdown"), y_fdown = 5, y_fup = 0.5)
```

    ## Warning in nb2listw(nb.object, style = style, glist = lapply(nb.dist, f1, : zero
    ## sum general weights

    ## Warning in nb2listw(nb.object, style = style, glist =
    ## lapply(nbdists(nb.object, : zero sum general weights

    ## Warning in nb2listw(nb.object, style = style, glist = lapply(nb.dist, f1, : zero
    ## sum general weights

    ## Warning in nb2listw(nb.object, style = style, glist =
    ## lapply(nbdists(nb.object, : zero sum general weights

``` r
UBA_MEM_env <- listw.select(UBA_env_FS, candidates = candidates_UBA_env, MEM.autocor = c("positive"), method = c("FWD"),
             MEM.all = FALSE, nperm = 10000, nperm.global = 10000, alpha = 0.05, p.adjust = FALSE, verbose = FALSE)
```

    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.114080 with 1 variables (> 0.114080)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.610050 with 4 variables (> 0.570723)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.532181 with 3 variables (> 0.504917)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.531064 with 2 variables (> 0.526818)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.584362 with 3 variables (> 0.552325)
    ## Procedure stopped (adjR2thresh criteria) adjR2cum = 0.544172 with 3 variables (> 0.543362)

``` r
#There was no significant spatial structure in the environmental variables, thus we used the one selected for the species
UBA_env_spa_sig <- envspace.test(spe = UBA_pa, env = UBA_env_FS,   coord = UBA_coord,  
              MEM.spe =  data.frame(UBA_MEM$best$MEM.select), 
              listw.env = candidates_UBA_env[[UBA_MEM$best.id]], MEM.autocor = c("positive"),
              regular = FALSE, nperm = 10000, MSR.method = "singleton", alpha = 0.05)


BER_varpart <- var_partitioning(BER_pa, BER_env_FS, spa = BER_MEM_FS, allow_negative_r2 = T, percent_r2 = F)
ITA_varpart <- var_partitioning(ITA_pa, ITA_env_FS, spa = ITA_MEM_FS, allow_negative_r2 = T, percent_r2 = F)


ST_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.197  5 1.3440574 0.2655
    ## Env           0.213  3 1.6319731 0.1014
    ## Clim             NA NA        NA     NA
    ## Spa          -0.085  2 0.7256147 0.8109
    ## Pure_Env      0.282  3 1.5862081 0.2532
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa     -0.016  2 0.9605140 0.5131
    ## Env_Spa      -0.069 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         0.803 NA        NA     NA

``` r
IC_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.291  4 2.1287607 0.0037
    ## Env           0.292  2 3.2736653 0.0001
    ## Clim             NA NA        NA     NA
    ## Spa          -0.041  2 0.7833254 0.7138
    ## Pure_Env      0.332  2 3.1073625 0.0017
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa     -0.001  2 0.9906546 0.5023
    ## Env_Spa      -0.040 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         0.709 NA        NA     NA

``` r
NI_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.140  5 1.2271911 0.2881
    ## Env           0.170  4 1.3581374 0.1512
    ## Clim             NA NA        NA     NA
    ## Spa           0.159  1 2.3197422 0.0181
    ## Pure_Env     -0.019  4 0.9668644 0.5475
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa     -0.030  1 0.8944824 0.5097
    ## Env_Spa       0.189 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         0.860 NA        NA     NA

``` r
NI_env_spa_sig
```

    ## NULL

``` r
MD_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.282  5 1.5500378 0.1596
    ## Env           0.139  4 1.2819349 0.2647
    ## Clim             NA NA        NA     NA
    ## Spa          -0.005  1 0.9665789 0.3754
    ## Pure_Env      0.287  4 1.5993494 0.1710
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa      0.143  1 1.5988564 0.2754
    ## Env_Spa      -0.148 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         0.718 NA        NA     NA

``` r
JA_varpart
```

    ##              Adj_R2 Df        F      p
    ## All           0.425  3 3.220404 0.0025
    ## Env           0.210  1 3.393898 0.0012
    ## Clim             NA NA       NA     NA
    ## Spa           0.314  2 3.060817 0.0057
    ## Pure_Env      0.111  1 2.354789 0.0562
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.215  2 2.498105 0.0374
    ## Env_Spa       0.099 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.575 NA       NA     NA

``` r
JA_env_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.09888151 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.4524548 
    ## Alternative hypothesis: greater 
    ## 
    ##     Std.Obs Expectation    Variance 
    ## 0.138770640 0.092678692 0.001997944

``` r
UBA_varpart
```

    ##              Adj_R2 Df        F      p
    ## All           0.151  3 2.308338 0.0005
    ## Env           0.095  1 3.306803 0.0010
    ## Clim             NA NA       NA     NA
    ## Spa           0.117  2 2.450992 0.0024
    ## Pure_Env      0.035  1 1.821646 0.0442
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.056  2 1.699031 0.0228
    ## Env_Spa       0.060 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.849 NA       NA     NA

``` r
UBA_env_spa_sig
```

    ## Monte-Carlo test
    ## Call: as.randtest(sim = E.b, obs = R2.b, alter = alternative)
    ## 
    ## Observation: 0.06004112 
    ## 
    ## Based on 10000 replicates
    ## Simulated p-value: 0.02929707 
    ## Alternative hypothesis: greater 
    ## 
    ##      Std.Obs  Expectation     Variance 
    ## 1.7941032784 0.0069554952 0.0008755057

``` r
BER_varpart
```

    ##              Adj_R2 Df         F      p
    ## All          -0.085  4 0.7852969 0.8008
    ## Env          -0.044  3 0.8438125 0.6890
    ## Clim             NA NA        NA     NA
    ## Spa          -0.012  1 0.8660103 0.5708
    ## Pure_Env     -0.072  3 0.7776483 0.7545
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa     -0.040  1 0.7035544 0.6956
    ## Env_Spa       0.028 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         1.085 NA        NA     NA

``` r
ITA_varpart
```

    ##              Adj_R2 Df         F      p
    ## All           0.065  7 1.1390639 0.3520
    ## Env           0.048  5 1.1420928 0.3372
    ## Clim             NA NA        NA     NA
    ## Spa          -0.008  2 0.9443449 0.5007
    ## Pure_Env      0.073  5 1.1874488 0.3126
    ## Pure_Clim        NA NA        NA     NA
    ## Pure_Spa      0.017  2 1.0804479 0.3919
    ## Env_Spa      -0.025 NA        NA     NA
    ## Env_Clim         NA NA        NA     NA
    ## Spa_Clim         NA NA        NA     NA
    ## Spa_Clim_Env     NA NA        NA     NA
    ## Resid         0.935 NA        NA     NA

   

Small Extent not allowing negative fractions

``` r
Broad_varpart2 <- var_partitioning(Y = Broad_pa, 
                                       env = Broad_env_FS,
                                       clim = Broad_clim_FS,
                                       spa = data.frame(Broad_MEM$best$MEM.select), percent_r2 = F, allow_negative_r2 = F)

DRF_varpart2 <- var_partitioning(Y = DRF_pa, 
                                      env = DRF_env_FS,
                                      clim = DRF_clim_FS,
                                      spa = DRF_MEM_FS, percent_r2 = F, allow_negative_r2 = F)

SSF_varpart2 <- var_partitioning(Y = SSF_pa, 
                                      env = SSF_env_FS,
                                      clim = SSF_clim_FS,
                                      spa = SSF_MEM_FS, percent_r2 = F, allow_negative_r2 = F)



ST_varpart2 <- var_partitioning(ST_pa, ST_env_FS, spa = ST_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
IC_varpart2 <- var_partitioning(IC_pa, IC_env_FS, spa = IC_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
NI_varpart2 <- var_partitioning(NI_pa, NI_env_FS, spa = NI_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
MD_varpart2 <- var_partitioning(MD_pa, MD_env_FS, spa = MD_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
JA_varpart2 <- var_partitioning(JA_pa, JA_env_FS, spa = JA_MEM_FS, allow_negative_r2 = F, percent_r2 = F)

UBA_varpart2 <- var_partitioning(UBA_pa, UBA_env_FS, spa = UBA_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
BER_varpart2 <- var_partitioning(BER_pa, BER_env_FS, spa = BER_MEM_FS, allow_negative_r2 = F, percent_r2 = F)
ITA_varpart2 <- var_partitioning(ITA_pa, ITA_env_FS, spa = ITA_MEM_FS, allow_negative_r2 = F, percent_r2 = F)

ST_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.213  3 1.631973 0.0979
    ## Env           0.213  3 1.631973 0.0979
    ## Clim             NA NA       NA     NA
    ## Spa           0.000 NA       NA     NA
    ## Pure_Env      0.213  3 1.631973 0.0979
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.000 NA       NA     NA
    ## Env_Spa       0.000 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.787 NA       NA     NA

``` r
IC_varpart2
```

    ##              Adj_R2 Df        F     p
    ## All           0.292  2 3.273665 1e-04
    ## Env           0.292  2 3.273665 1e-04
    ## Clim             NA NA       NA    NA
    ## Spa           0.000 NA       NA    NA
    ## Pure_Env      0.292  2 3.273665 1e-04
    ## Pure_Clim        NA NA       NA    NA
    ## Pure_Spa      0.000 NA       NA    NA
    ## Env_Spa       0.000 NA       NA    NA
    ## Env_Clim         NA NA       NA    NA
    ## Spa_Clim         NA NA       NA    NA
    ## Spa_Clim_Env     NA NA       NA    NA
    ## Resid         0.708 NA       NA    NA

``` r
NI_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.170  4 1.358137 0.1485
    ## Env           0.170  4 1.358137 0.1485
    ## Clim             NA NA       NA     NA
    ## Spa           0.159  1 2.319742 0.0149
    ## Pure_Env      0.011 NA       NA     NA
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.000 NA       NA     NA
    ## Env_Spa       0.159 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.830 NA       NA     NA

``` r
MD_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.139  4 1.281935 0.2652
    ## Env           0.139  4 1.281935 0.2652
    ## Clim             NA NA       NA     NA
    ## Spa           0.000 NA       NA     NA
    ## Pure_Env      0.139  4 1.281935 0.2652
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.000 NA       NA     NA
    ## Env_Spa       0.000 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.861 NA       NA     NA

``` r
JA_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.425  3 3.220404 0.0013
    ## Env           0.210  1 3.393898 0.0022
    ## Clim             NA NA       NA     NA
    ## Spa           0.314  2 3.060817 0.0054
    ## Pure_Env      0.111  1 2.354789 0.0575
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.215  2 2.498105 0.0348
    ## Env_Spa       0.099 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.575 NA       NA     NA

``` r
UBA_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.151  3 2.308338 0.0009
    ## Env           0.095  1 3.306803 0.0013
    ## Clim             NA NA       NA     NA
    ## Spa           0.117  2 2.450992 0.0032
    ## Pure_Env      0.035  1 1.821646 0.0468
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.056  2 1.699031 0.0272
    ## Env_Spa       0.060 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.849 NA       NA     NA

``` r
BER_varpart2
```

    ##              Adj_R2 Df  F  p
    ## All               0 NA NA NA
    ## Env               0 NA NA NA
    ## Clim             NA NA NA NA
    ## Spa               0 NA NA NA
    ## Pure_Env          0 NA NA NA
    ## Pure_Clim        NA NA NA NA
    ## Pure_Spa          0 NA NA NA
    ## Env_Spa           0 NA NA NA
    ## Env_Clim         NA NA NA NA
    ## Spa_Clim         NA NA NA NA
    ## Spa_Clim_Env     NA NA NA NA
    ## Resid             1 NA NA NA

``` r
ITA_varpart2
```

    ##              Adj_R2 Df        F      p
    ## All           0.048  5 1.142093 0.3252
    ## Env           0.048  5 1.142093 0.3252
    ## Clim             NA NA       NA     NA
    ## Spa           0.000 NA       NA     NA
    ## Pure_Env      0.048  5 1.142093 0.3252
    ## Pure_Clim        NA NA       NA     NA
    ## Pure_Spa      0.000 NA       NA     NA
    ## Env_Spa       0.000 NA       NA     NA
    ## Env_Clim         NA NA       NA     NA
    ## Spa_Clim         NA NA       NA     NA
    ## Spa_Clim_Env     NA NA       NA     NA
    ## Resid         0.952 NA       NA     NA

   

Constructing a matrix with all R2 values (with negative values)

``` r
Varpart_plot_neg <- data.frame(Broad_varpart[,1],
                           SSF_varpart[,1],
                           DRF_varpart[,1],
                           ST_varpart[,1],
                           IC_varpart[,1],
                           NI_varpart[,1],
                           MD_varpart[,1],
                           JA_varpart[,1],
                           UBA_varpart[,1],
                           BER_varpart[,1],
                           ITA_varpart[,1])

colnames(Varpart_plot_neg) <- c("Broad", "SSF", "DRF", "Santa Fé do Sul", "Icém", "Nova Itapirema", "Morro do Diabo", "Jataí", "Ubatuba", "Bertioga", "Itanhaém")

for(i in 1:dim(Varpart_plot_neg)[1]){
  for(j in 1:dim(Varpart_plot_neg)[2]){
    if(is.na(Varpart_plot_neg[i,j])){Varpart_plot_neg[i,j] <- 0}
  }
}

Varpart_plot_neg
```

    ##              Broad   SSF   DRF Santa Fé do Sul   Icém Nova Itapirema
    ## All          0.411 0.346 0.191           0.197  0.291          0.140
    ## Env          0.166 0.136 0.073           0.213  0.292          0.170
    ## Clim         0.177 0.171 0.049           0.000  0.000          0.000
    ## Spa          0.369 0.294 0.165          -0.085 -0.041          0.159
    ## Pure_Env     0.012 0.034 0.016           0.282  0.332         -0.019
    ## Pure_Clim    0.026 0.009 0.006           0.000  0.000          0.000
    ## Pure_Spa     0.168 0.078 0.091          -0.016 -0.001         -0.030
    ## Env_Spa      0.054 0.063 0.035          -0.069 -0.040          0.189
    ## Env_Clim     0.004 0.010 0.004           0.000  0.000          0.000
    ## Spa_Clim     0.051 0.124 0.021           0.000  0.000          0.000
    ## Spa_Clim_Env 0.096 0.029 0.018           0.000  0.000          0.000
    ## Resid        0.589 0.654 0.809           0.803  0.709          0.860
    ##              Morro do Diabo Jataí Ubatuba Bertioga Itanhaém
    ## All                   0.282 0.425   0.151   -0.085    0.065
    ## Env                   0.139 0.210   0.095   -0.044    0.048
    ## Clim                  0.000 0.000   0.000    0.000    0.000
    ## Spa                  -0.005 0.314   0.117   -0.012   -0.008
    ## Pure_Env              0.287 0.111   0.035   -0.072    0.073
    ## Pure_Clim             0.000 0.000   0.000    0.000    0.000
    ## Pure_Spa              0.143 0.215   0.056   -0.040    0.017
    ## Env_Spa              -0.148 0.099   0.060    0.028   -0.025
    ## Env_Clim              0.000 0.000   0.000    0.000    0.000
    ## Spa_Clim              0.000 0.000   0.000    0.000    0.000
    ## Spa_Clim_Env          0.000 0.000   0.000    0.000    0.000
    ## Resid                 0.718 0.575   0.849    1.085    0.935

     
   

Constructing a matrix with all p values (with negative R2 values)

``` r
Varpart_plot_p <- data.frame(Broad_varpart[,4],
                           SSF_varpart[,4],
                           DRF_varpart[,4],
                           ST_varpart[,4],
                           IC_varpart[,4],
                           NI_varpart[,4],
                           MD_varpart[,4],
                           JA_varpart[,4],
                           UBA_varpart[,4],
                           BER_varpart[,4],
                           ITA_varpart[,4])

colnames(Varpart_plot_p) <- c("Broad", "SSF", "DRF", "Santa Fé do Sul", "Icém", "Nova Itapirema", "Morro do Diabo", "Jataí", "Ubatuba", "Bertioga", "Itanhaém")

Varpart_plot_p
```

    ##               Broad    SSF    DRF Santa Fé do Sul   Icém Nova Itapirema
    ## All          0.0001 0.0001 0.0001          0.2655 0.0037         0.2881
    ## Env          0.0001 0.0001 0.0003          0.1014 0.0001         0.1512
    ## Clim         0.0001 0.0001 0.0003              NA     NA             NA
    ## Spa          0.0001 0.0001 0.0001          0.8109 0.7138         0.0181
    ## Pure_Env     0.0122 0.0211 0.0827          0.2532 0.0017         0.5475
    ## Pure_Clim    0.0001 0.1878 0.1928              NA     NA             NA
    ## Pure_Spa     0.0001 0.0005 0.0003          0.5131 0.5023         0.5097
    ## Env_Spa          NA     NA     NA              NA     NA             NA
    ## Env_Clim         NA     NA     NA              NA     NA             NA
    ## Spa_Clim         NA     NA     NA              NA     NA             NA
    ## Spa_Clim_Env     NA     NA     NA              NA     NA             NA
    ## Resid            NA     NA     NA              NA     NA             NA
    ##              Morro do Diabo  Jataí Ubatuba Bertioga Itanhaém
    ## All                  0.1596 0.0025  0.0005   0.8008   0.3520
    ## Env                  0.2647 0.0012  0.0010   0.6890   0.3372
    ## Clim                     NA     NA      NA       NA       NA
    ## Spa                  0.3754 0.0057  0.0024   0.5708   0.5007
    ## Pure_Env             0.1710 0.0562  0.0442   0.7545   0.3126
    ## Pure_Clim                NA     NA      NA       NA       NA
    ## Pure_Spa             0.2754 0.0374  0.0228   0.6956   0.3919
    ## Env_Spa                  NA     NA      NA       NA       NA
    ## Env_Clim                 NA     NA      NA       NA       NA
    ## Spa_Clim                 NA     NA      NA       NA       NA
    ## Spa_Clim_Env             NA     NA      NA       NA       NA
    ## Resid                    NA     NA      NA       NA       NA

     
   

Constructing a matrix with all R2 values (without negative values)

``` r
Varpart_plot <- data.frame(Broad_varpart2[,1],
                           SSF_varpart2[,1],
                           DRF_varpart2[,1],
                           ST_varpart2[,1],
                           IC_varpart2[,1],
                           NI_varpart2[,1],
                           MD_varpart2[,1],
                           JA_varpart2[,1],
                           UBA_varpart2[,1],
                           BER_varpart2[,1],
                           ITA_varpart2[,1])

colnames(Varpart_plot) <- c("Broad", "SSF", "DRF", "Santa Fé do Sul", "Icém", "Nova Itapirema", "Morro do Diabo", "Jataí", "Ubatuba", "Bertioga", "Itanhaém")

for(i in 1:dim(Varpart_plot)[1]){
  for(j in 1:dim(Varpart_plot)[2]){
    if(is.na(Varpart_plot[i,j])){Varpart_plot[i,j] <- 0}
  }
}

Varpart_plot
```

    ##              Broad   SSF   DRF Santa Fé do Sul  Icém Nova Itapirema
    ## All          0.411 0.346 0.191           0.213 0.292          0.170
    ## Env          0.166 0.136 0.073           0.213 0.292          0.170
    ## Clim         0.177 0.171 0.049           0.000 0.000          0.000
    ## Spa          0.369 0.294 0.165           0.000 0.000          0.159
    ## Pure_Env     0.012 0.034 0.016           0.213 0.292          0.011
    ## Pure_Clim    0.026 0.009 0.006           0.000 0.000          0.000
    ## Pure_Spa     0.168 0.078 0.091           0.000 0.000          0.000
    ## Env_Spa      0.054 0.063 0.035           0.000 0.000          0.159
    ## Env_Clim     0.004 0.010 0.004           0.000 0.000          0.000
    ## Spa_Clim     0.051 0.124 0.021           0.000 0.000          0.000
    ## Spa_Clim_Env 0.096 0.029 0.018           0.000 0.000          0.000
    ## Resid        0.589 0.654 0.809           0.787 0.708          0.830
    ##              Morro do Diabo Jataí Ubatuba Bertioga Itanhaém
    ## All                   0.139 0.425   0.151        0    0.048
    ## Env                   0.139 0.210   0.095        0    0.048
    ## Clim                  0.000 0.000   0.000        0    0.000
    ## Spa                   0.000 0.314   0.117        0    0.000
    ## Pure_Env              0.139 0.111   0.035        0    0.048
    ## Pure_Clim             0.000 0.000   0.000        0    0.000
    ## Pure_Spa              0.000 0.215   0.056        0    0.000
    ## Env_Spa               0.000 0.099   0.060        0    0.000
    ## Env_Clim              0.000 0.000   0.000        0    0.000
    ## Spa_Clim              0.000 0.000   0.000        0    0.000
    ## Spa_Clim_Env          0.000 0.000   0.000        0    0.000
    ## Resid                 0.861 0.575   0.849        1    0.952

     
   

Plotting the Variation Partioning as barplots.

``` r
Varpart_barplot <- Varpart_plot[-c(1:4),]
Varpart_barplot_break <- Varpart_barplot; Varpart_barplot_break[8,] <- Varpart_barplot_break[8,]-0.2
Varpart_barplot_break <- Varpart_barplot_break[c(1,4,3,6,2,5,7,8),]

par(mfrow = c(1,1))
barplot(as.matrix(Varpart_barplot_break), axes = F, col = c(Pure_Env = "gold",
                                                            Env_Spa = mix_color(alpha = 0.4,"gold","cornflowerblue"),
                                                            Pure_Spa = "cornflowerblue",
                                                            Spa_Clim  = mix_color(alpha = 0.4,"brown1","cornflowerblue"),
                                                            Pure_Clim = "brown1",
                                                            Env_Clim  = mix_color(alpha = 0.4,"brown1","gold"),
                                                            Spa_Clim_Env = mix_color(alpha = 0.7,"brown1","grey"),
                                                            Resid = "grey80"), space = c(0,2,1,2,1,1,1,1,2,1,1), border = "white",
        legend.text = c("Environment","Environment-Space","Space","Climate-Space","Climate","Climate-Environment","All three","Residual"), ylim = c(0,0.8),
        args.legend = list(x = -0.5,y = 0.9, yjust = 1, xjust = 0, horiz = F, ncol = 4, box.col = "transparent",text.width = 5, border = "white"), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25)
axis(2, at = c(0,0.1,0.2,0.3,0.4,0.5,0.6, 0.7, 0.8), labels = c("0 %","10 %","20 %","30 %","40 %","50 %","60 %","70%","100%"))
axis.break(2, 0.75, style = "slash") 
axis(1,at = c(0.5,4.5,16),line = 2, labels =c("Broad","Intermediate","Fine"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
axis(1,at = c(12.5,21.5),line = 0.5, labels =c("SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.1)
text(c(0.5,   3.5,5.5   ,8.5,10.5,12.5,14.5,16.5,    19.5,21.5,23.5),rep(0.79,11), labels =c("All","SSF","DRF","Santa Fé do Sul","Icém","Nova Itapirema","Morro do Diabo","Jataí","Ubatuba","Bertioga","Itanhaém"), srt = 90, adj = 1, col = "grey45")
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-22-1.png" width="800" height="700" style="display: block; margin: auto;" />

``` r
#text(c(10.5,12.5,16.5,16.5,19.5),c(0.1,0.08,0.17,0.59,0.03), labels =c("*","*","*","*","*","*"), adj = 0.5, col = "white", cex = 2)
#text(c(0.5,0.5,0.5),c(0.01,0.12,0.299), labels =c("*","*","*"), adj = 0.5, col = "white", cex = 2)
#text(c(3.5,3.5),c(0.02,0.1), labels =c("*","*","*"), adj = 0.5, col = "white", cex = 2)
#text(c(5.5,5.5),c(0.01,0.065), labels =c("*","*","*"), adj = 0.5, col = "white", cex = 2)

#text(c(0.5,0.3,0.7,0.5),c(0.24,0.35,0.35,0.035), labels =c("*","*","*","*"), adj = 0.5, col = c("brown1","brown1","gold","gold"), cex = 2)

#text(c(3.5,3.3,3.7),c(0.22,0.3075,0.3075), labels =c("*","*","*"), adj = 0.5, col = c("brown1","brown1","gold"), cex = 2)

#text(c(5.5,5.3, 5.7),c(0.115,0.145,0.145), labels =c("*","*","*"), adj = 0.5, col = c("brown1","brown1","gold"), cex = 2)


###Mudar para simbolos

#cada simbolo representa env, clim e spa
#os shared serão os simbolos de spa com cores dos outros
```

The asterisks represent significant fractions.      

Creating some new important fractions

``` r
Env_Spatially_Structured <- Varpart_plot["Env_Spa",] + Varpart_plot["Spa_Clim_Env",]
Env_Non_Spatially_Structured <- Varpart_plot["Pure_Env",] + Varpart_plot["Env_Clim",]
Clim_Spatially_Structured <- Varpart_plot["Spa_Clim",] + Varpart_plot["Spa_Clim_Env",]
Clim_Non_Spatially_Structured <- Varpart_plot["Pure_Clim",] + Varpart_plot["Env_Clim",]
Non_Spatially_Climate_Environment <- Varpart_plot["Pure_Env",] + Varpart_plot["Pure_Clim",] + Varpart_plot["Env_Clim",]
Spatially_Climate_Environment <- Varpart_plot["Spa_Clim",] + Varpart_plot["Env_Spa",] + Varpart_plot["Spa_Clim_Env",]
Total_Climate_Environment <- Spatially_Climate_Environment + Non_Spatially_Climate_Environment
#Total_Climate_Environment["Nova Itapirema"] <- 0

Varpart_plot <- rbind(Varpart_plot,
                      Env_Spatially_Structured = Env_Spatially_Structured,
                      Env_Non_Spatially_Structured = Env_Non_Spatially_Structured,
                      Clim_Spatially_Structured = Clim_Spatially_Structured,
                      Clim_Non_Spatially_Structured = Clim_Non_Spatially_Structured,
                      Non_Spatially_Climate_Environment = Non_Spatially_Climate_Environment,
                      Spatially_Climate_Environment = Spatially_Climate_Environment,
                      Total_Climate_Environment = Total_Climate_Environment)
```

Means of fractions in all extents

``` r
se <- function(x){
  sd(x)/sqrt(length(x))
}


Means_spatial_extent <- data.frame(Broad_Extent = Varpart_plot[,1],
           Intermediate_Extent = apply(Varpart_plot[,2:3],1,mean),
           Small_Extent = apply(Varpart_plot[,4:11],1,mean))

se_spatial_extent <- data.frame(Intermediate_Extent = apply(Varpart_plot[,2:3],1,se),
           Small_Extent = apply(Varpart_plot[,4:11],1,se))

se_spatial_extent_upper <- Means_spatial_extent[,2:3] + se_spatial_extent
se_spatial_extent_lower <- Means_spatial_extent[,2:3] - se_spatial_extent

Means_spatial_extent
```

    ##                                   Broad_Extent Intermediate_Extent Small_Extent
    ## All                                      0.411              0.2685     0.179750
    ## Env                                      0.166              0.1045     0.145875
    ## Clim                                     0.177              0.1100     0.000000
    ## Spa                                      0.369              0.2295     0.073750
    ## Pure_Env                                 0.012              0.0250     0.106125
    ## Pure_Clim                                0.026              0.0075     0.000000
    ## Pure_Spa                                 0.168              0.0845     0.033875
    ## Env_Spa                                  0.054              0.0490     0.039750
    ## Env_Clim                                 0.004              0.0070     0.000000
    ## Spa_Clim                                 0.051              0.0725     0.000000
    ## Spa_Clim_Env                             0.096              0.0235     0.000000
    ## Resid                                    0.589              0.7315     0.820250
    ## Env_Spatially_Structured                 0.150              0.0725     0.039750
    ## Env_Non_Spatially_Structured             0.016              0.0320     0.106125
    ## Clim_Spatially_Structured                0.147              0.0960     0.000000
    ## Clim_Non_Spatially_Structured            0.030              0.0145     0.000000
    ## Non_Spatially_Climate_Environment        0.042              0.0395     0.106125
    ## Spatially_Climate_Environment            0.201              0.1450     0.039750
    ## Total_Climate_Environment                0.243              0.1845     0.145875

 

Means of fractions in small extent

``` r
Means_small_extent <- data.frame(SSF = apply(Varpart_plot[,4:8],1,mean),
           DRF = apply(Varpart_plot[,9:11],1,mean))

se_small_extent <- data.frame(SSF = apply(Varpart_plot[,4:8],1,se),
           DRF = apply(Varpart_plot[,9:11],1,se))

se_small_extent_upper <- Means_small_extent + se_small_extent
se_small_extent_lower <- Means_small_extent - se_small_extent

Means_small_extent
```

    ##                                      SSF        DRF
    ## All                               0.2478 0.06633333
    ## Env                               0.2048 0.04766667
    ## Clim                              0.0000 0.00000000
    ## Spa                               0.0946 0.03900000
    ## Pure_Env                          0.1532 0.02766667
    ## Pure_Clim                         0.0000 0.00000000
    ## Pure_Spa                          0.0430 0.01866667
    ## Env_Spa                           0.0516 0.02000000
    ## Env_Clim                          0.0000 0.00000000
    ## Spa_Clim                          0.0000 0.00000000
    ## Spa_Clim_Env                      0.0000 0.00000000
    ## Resid                             0.7522 0.93366667
    ## Env_Spatially_Structured          0.0516 0.02000000
    ## Env_Non_Spatially_Structured      0.1532 0.02766667
    ## Clim_Spatially_Structured         0.0000 0.00000000
    ## Clim_Non_Spatially_Structured     0.0000 0.00000000
    ## Non_Spatially_Climate_Environment 0.1532 0.02766667
    ## Spatially_Climate_Environment     0.0516 0.02000000
    ## Total_Climate_Environment         0.2048 0.04766667

 

``` r
par(mfrow = c(1,3))


barplot(as.matrix(Means_spatial_extent[c(16,15),]), axes = F, col = c("brown1",mix_color(alpha = 0.5,"brown1","grey25")), space = c(0,1,1),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Climate")
arrows(y0 = c(as.matrix(se_spatial_extent_lower[3,])),
       y1 = c(as.matrix(se_spatial_extent_upper[3,])),
       x1 = c(2.5,4.5),
       x0 = c(2.5,4.5),  code = 3, angle = 90, length = 0.1, lwd = 1)
```

    ## Warning in arrows(y0 = c(as.matrix(se_spatial_extent_lower[3, ])), y1 =
    ## c(as.matrix(se_spatial_extent_upper[3, : zero-length arrow is of indeterminate
    ## angle and so skipped

``` r
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(0.5,2.5,4.5),line = 0, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)


barplot(as.matrix(Means_spatial_extent[c(14,13),]), axes = F, col = c("gold",mix_color(alpha = 0.5,"gold","grey25")), space = c(0,1,1),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Local Environment")
arrows(y0 = c(as.matrix(se_spatial_extent_lower[2,])),
       y1 = c(as.matrix(se_spatial_extent_upper[2,])),
       x1 = c(2.5,4.5),
       x0 = c(2.5,4.5),  code = 3, angle = 90, length = 0.1, lwd = 1)
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(0.5,2.5,4.5),line = 0, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)



barplot(as.matrix(Means_spatial_extent[c(7,18),]), axes = F, col = c("cornflowerblue",mix_color(alpha = 0.5,"cornflowerblue","grey25")), space = c(0,1,1),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Space")
arrows(y0 = c(as.matrix(se_spatial_extent_lower[4,])),
       y1 = c(as.matrix(se_spatial_extent_upper[4,])),
       x1 = c(2.5,4.5),
       x0 = c(2.5,4.5),  code = 3, angle = 90, length = 0.1, lwd = 1)
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(0.5,2.5,4.5),line = 0, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-26-1.png" width="1000" height="500" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,3))

barplot(as.matrix(Varpart_plot[c("Clim_Non_Spatially_Structured", "Clim_Spatially_Structured"),c("SSF","DRF")]),
        axes = F, col = c("brown1",mix_color(alpha = 0.5,"brown1","grey25")), space = c(0,1,2,1),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Climate")
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(1.5,6.5),line = 1, labels =c("Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
axis(1,at = c(0.5,2.5,5.5,7.5),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)


barplot(as.matrix(data.frame(Varpart_plot[c("Env_Non_Spatially_Structured", "Env_Spatially_Structured"),c("SSF","DRF")],
        Means_small_extent[c("Env_Non_Spatially_Structured","Env_Spatially_Structured"),])),
        axes = F, col = c("gold",mix_color(alpha = 0.5,"gold","grey25")), space = c(0,1,2,1),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Local Environment")
arrows(y0 = c(as.matrix(se_small_extent_lower["Env",])),
       y1 = c(as.matrix(se_small_extent_upper["Env",])),
       x1 = c(5.5,7.5),
       x0 = c(5.5,7.5),  code = 3, angle = 90, length = 0.1, lwd = 1)
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(1.5,6.5),line = 1, labels =c("Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
axis(1,at = c(0.5,2.5,5.5,7.5),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)



barplot(as.matrix(data.frame(Varpart_plot[c("Pure_Spa", "Spatially_Climate_Environment"),c("SSF","DRF")],
  Means_small_extent[c("Pure_Spa","Spatially_Climate_Environment"),])),
  axes = F, col = c("cornflowerblue",mix_color(alpha = 0.5,"cornflowerblue","grey25")), space = c(0,1,2,1),
  border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Space")
arrows(y0 = c(as.matrix(se_small_extent_lower["Spa",])),
       y1 = c(as.matrix(se_small_extent_upper["Spa",])),
       x1 = c(5.5,7.5),
       x0 = c(5.5,7.5),  code = 3, angle = 90, length = 0.1, lwd = 1)
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))
axis(1,at = c(1.5,6.5),line = 1, labels =c("Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
axis(1,at = c(0.5,2.5,5.5,7.5),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-27-1.png" width="1000" height="500" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,3))

barplot(as.matrix(Varpart_plot[c("Clim_Non_Spatially_Structured", "Clim_Spatially_Structured"),c("Broad","SSF","DRF")]),
        axes = F, col = c("brown1",mix_color(alpha = 0.5,"brown1","grey25")), space = c(0.75,1.5,0.5,1.5,0.5),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Climate", xlim = c(0,9.75) )
```

    ## Warning in space + width: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

    ## Warning in w.r - delta: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

    ## Warning in w.m - delta: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

    ## Warning in space + width: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

    ## Warning in w.r - delta: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

    ## Warning in w.m - delta: comprimento do objeto maior não é múltiplo do
    ## comprimento do objeto menor

``` r
arrows(y0 = c(as.matrix(se_small_extent_lower["Clim",])),
       y1 = c(as.matrix(se_small_extent_upper["Clim",])),
       x1 = c(7.75,9.25),
       x0 = c(7.75,9.25),  code = 3, angle = 90, length = 0.1, lwd = 1)
```

    ## Warning in arrows(y0 = c(as.matrix(se_small_extent_lower["Clim", ])),
    ## y1 = c(as.matrix(se_small_extent_upper["Clim", : zero-length arrow is of
    ## indeterminate angle and so skipped

    ## Warning in arrows(y0 = c(as.matrix(se_small_extent_lower["Clim", ])),
    ## y1 = c(as.matrix(se_small_extent_upper["Clim", : zero-length arrow is of
    ## indeterminate angle and so skipped

``` r
axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))

axis(1,at = c(1.25,4.5,8.5),line = 1, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)
axis(1,at = c(3.75,5.25,7.75,9.25),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)



barplot(as.matrix(data.frame(Varpart_plot[c("Env_Non_Spatially_Structured", "Env_Spatially_Structured"),c("Broad","SSF","DRF")],
        Means_small_extent[c("Env_Non_Spatially_Structured","Env_Spatially_Structured"),])),
        axes = F, col = c("gold",mix_color(alpha = 0.5,"gold","grey25")), space = c(0.75,1.5,0.5,1.5,0.5),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "local Environment", xlim = c(0,9.75) )

arrows(y0 = c(as.matrix(se_small_extent_lower["Env",])),
       y1 = c(as.matrix(se_small_extent_upper["Env",])),
       x1 = c(7.75,9.25),
       x0 = c(7.75,9.25),  code = 3, angle = 90, length = 0.1, lwd = 1)

axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))

axis(1,at = c(1.25,4.5,8.5),line = 1, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)
axis(1,at = c(3.75,5.25,7.75,9.25),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)




barplot(as.matrix(data.frame(Varpart_plot[c("Pure_Spa", "Spatially_Climate_Environment"),c("Broad","SSF","DRF")],
        Means_small_extent[c("Pure_Spa","Spatially_Climate_Environment"),])),
        axes = F, col = c("cornflowerblue",mix_color(alpha = 0.5,"cornflowerblue","grey25")), space = c(0.75,1.5,0.5,1.5,0.5),
        border = "white", ylim = c(0,0.4), axisnames= F, ylab = "Adjusted R²", cex.lab = 1.25, main = "Space", xlim = c(0,9.75) )

arrows(y0 = c(as.matrix(se_small_extent_lower["Spa",])),
       y1 = c(as.matrix(se_small_extent_upper["Spa",])),
       x1 = c(7.75,9.25),
       x0 = c(7.75,9.25),  code = 3, angle = 90, length = 0.1, lwd = 1)

axis(2, at = c(0,0.05,0.1,0.15,0.2, 0.25,0.3,0.35,0.4), labels = c("0 %","5 %","10 %","15 %","20 %", "25 %", "30 %", "35 %", "40 %"))

axis(1,at = c(1.25,4.5,8.5),line = 1, labels =c("Broad","Intermediate","Small"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)
axis(1,at = c(3.75,5.25,7.75,9.25),line = -0.5, labels =c("SSF","DRF", "SSF","DRF"), tick = F,las = 1, hadj = 0.5, cex.axis = 1.3, gap.axis = -1)
```

<img src="Variation_Partitioning_files/figure-gfm/unnamed-chunk-28-1.png" width="1000" height="500" style="display: block; margin: auto;" />
