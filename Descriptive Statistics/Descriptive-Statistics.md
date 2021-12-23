Descriptive Statistics
================
Rodolfo Pelinson
13/04/2021

Here you can find the number and identity of species at all spatial
extents, localities and ecoregions. We also show here which species are
shared between ecoregions and localities.

We also show maximum, minimum, and mean values of each environmental
variable considered in all analyses at all spatial extents, localities
and ecoregions.      

First we need preparare all data matrices from the main dataset. These
are prepared sourcing the “Loading\_data.R” file in the Auxiliary
Scripts folder.

``` r
library(AtlanticForestMetacommunity)
source("Loading_data.R")
```

    ## Error in get(genname, envir = envir) : 
    ##   objeto 'testthat_print' não encontrado

## Number of Species

### Large Extent

``` r
ncol(Broad_pa_orig)
```

    ## [1] 57

   
 

### Intermediate Extent

#### DRF

``` r
ncol(DRF_pa_orig)
```

    ## [1] 28

 

#### SSF

``` r
ncol(SSF_pa_orig)
```

    ## [1] 34

   

### Small Extent

#### SSF - Santa Fé do Sul

``` r
ncol(ST_pa_orig)
```

    ## [1] 14

 

#### SSF - Icém

``` r
ncol(IC_pa_orig)
```

    ## [1] 22

 

#### SSF - Nova Itapirema

``` r
ncol(NI_pa_orig)
```

    ## [1] 21

 

#### SSF - Morro do Diabo

``` r
ncol(MD_pa_orig)
```

    ## [1] 17

 

#### SSF - Jataí

``` r
ncol(JA_pa_orig)
```

    ## [1] 15

 

#### DRF - Ubatuba

``` r
ncol(UBA_pa_orig)
```

    ## [1] 23

 

#### DRF - Bertioga

``` r
ncol(BER_pa_orig)
```

    ## [1] 17

 

#### DRF - Itanhaém

``` r
ncol(ITA_pa_orig)
```

    ## [1] 21

         

## Species Shared

### Intermediate Extent

``` r
sp_comm <- na.omit(match(colnames(DRF_pa_orig), colnames(SSF_pa_orig)))
length(sp_comm)
```

    ## [1] 5

``` r
colnames(SSF_pa_orig)[sp_comm]
```

    ## [1] "Bfa" "Dmi" "Llt" "Pcu" "Ror"

   

### Small Extent

#### Santa Fé do Sul - Icém

``` r
sp_comm <- na.omit(match(colnames(ST_pa_orig), colnames(IC_pa_orig)))
length(sp_comm)
```

    ## [1] 13

``` r
colnames(IC_pa_orig)[sp_comm]
```

    ##  [1] "Bal"  "Dmu"  "Dna"  "Esp1" "Lfs"  "Lpo"  "Pcu"  "Pma"  "Pna"  "Rsc" 
    ## [11] "Sfs"  "Ssi"  "Tve"

   

#### Santa Fé do Sul - Nova Itapirema

``` r
sp_comm <- na.omit(match(colnames(ST_pa_orig), colnames(NI_pa_orig)))
length(sp_comm)
```

    ## [1] 13

``` r
colnames(NI_pa_orig)[sp_comm]
```

    ##  [1] "Bal"  "Bra"  "Dna"  "Esp1" "Lfs"  "Lpo"  "Pcu"  "Pma"  "Pna"  "Rsc" 
    ## [11] "Sfs"  "Ssi"  "Tve"

   

#### Santa Fé do Sul - Morro do Diabo

``` r
sp_comm <- na.omit(match(colnames(ST_pa_orig), colnames(MD_pa_orig)))
length(sp_comm)
```

    ## [1] 11

``` r
colnames(MD_pa_orig)[sp_comm]
```

    ##  [1] "Bra"  "Dna"  "Esp1" "Lfs"  "Lpo"  "Pcu"  "Pna"  "Rsc"  "Sfs"  "Ssi" 
    ## [11] "Tve"

   

#### Santa Fé do Sul - Jataí

``` r
sp_comm <- na.omit(match(colnames(ST_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 7

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dna" "Lfs" "Lpo" "Pcu" "Sfs" "Ssi"

   

#### Icém - Nova Itapirema

``` r
sp_comm <- na.omit(match(colnames(IC_pa_orig), colnames(NI_pa_orig)))
length(sp_comm)
```

    ## [1] 17

``` r
colnames(NI_pa_orig)[sp_comm]
```

    ##  [1] "Bal"  "Dmi"  "Dna"  "Esp1" "Esp2" "Lfs"  "Llb"  "Lpo"  "Pcu"  "Pma" 
    ## [11] "Pmy"  "Pna"  "Rsc"  "Sfs"  "Sfu"  "Ssi"  "Tve"

   

#### Icém - Morro do Diabo

``` r
sp_comm <- na.omit(match(colnames(IC_pa_orig), colnames(MD_pa_orig)))
length(sp_comm)
```

    ## [1] 12

``` r
colnames(MD_pa_orig)[sp_comm]
```

    ##  [1] "Dmi"  "Dna"  "Esp1" "Lfs"  "Lpo"  "Pcu"  "Pna"  "Rsc"  "Sfs"  "Sfu" 
    ## [11] "Ssi"  "Tve"

   

#### Icém - Jataí

``` r
sp_comm <- na.omit(match(colnames(IC_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 9

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dmi" "Dna" "Lfs" "Llb" "Lpo" "Pcu" "Sfs" "Ssi"

   

#### Nova Itapirema - Morro do Diabo

``` r
sp_comm <- na.omit(match(colnames(NI_pa_orig), colnames(MD_pa_orig)))
length(sp_comm)
```

    ## [1] 14

``` r
colnames(MD_pa_orig)[sp_comm]
```

    ##  [1] "Bra"  "Dmi"  "Dna"  "Esp1" "Lfs"  "Llt"  "Lpo"  "Pcu"  "Pna"  "Rsc" 
    ## [11] "Sfs"  "Sfu"  "Ssi"  "Tve"

   

#### Nova Itapirema - Jataí

``` r
sp_comm <- na.omit(match(colnames(NI_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 9

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dmi" "Dna" "Lfs" "Llb" "Lpo" "Pcu" "Sfs" "Ssi"

   

#### Morro do Diabo - Jataí

``` r
sp_comm <- na.omit(match(colnames(MD_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 8

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bfa" "Dmi" "Dna" "Lfs" "Lpo" "Pcu" "Sfs" "Ssi"

   

#### Ubatuba - Bertioga

``` r
sp_comm <- na.omit(match(colnames(UBA_pa_orig), colnames(BER_pa_orig)))
length(sp_comm)
```

    ## [1] 15

``` r
colnames(BER_pa_orig)[sp_comm]
```

    ##  [1] "Aeu" "Bab" "Bfa" "Dbe" "Del" "Dmi" "Ila" "Llt" "Ror" "Sal" "Sar" "Sha"
    ## [13] "Spe" "Str" "Tme"

   

#### Ubatuba - Itanhaém

``` r
sp_comm <- na.omit(match(colnames(UBA_pa_orig), colnames(ITA_pa_orig)))
length(sp_comm)
```

    ## [1] 17

``` r
colnames(ITA_pa_orig)[sp_comm]
```

    ##  [1] "Bab" "Bfa" "Cca" "Dbe" "Del" "Dmc" "Dmi" "Ila" "Llt" "Ror" "Sal" "Sar"
    ## [13] "Sha" "Sli" "Spe" "Str" "Tme"

   

#### Itanhaém - Bertioga

``` r
sp_comm <- na.omit(match(colnames(ITA_pa_orig), colnames(BER_pa_orig)))
length(sp_comm)
```

    ## [1] 15

``` r
colnames(BER_pa_orig)[sp_comm]
```

    ##  [1] "Bab" "Bfa" "Bse" "Dbe" "Del" "Dmi" "Ila" "Llt" "Ror" "Sal" "Sar" "Sha"
    ## [13] "Spe" "Str" "Tme"

         

## Means, maximum and minimum of predictor variables

 

### Intermediate Extent

#### DRF

``` r
data.frame(DRF_min = apply(DRF_clim, 2, min), DRF_mean = apply(DRF_clim, 2, mean), DRF_max = apply(DRF_clim, 2, max))
```

    ##             DRF_min DRF_mean DRF_max
    ## temp_Season    2194  2308.18    2603
    ## range_temp      143   161.42     173
    ## total_prec     1911  2279.84    2747
    ## prec_season      36    41.86      48

``` r
data.frame(DRF_min = apply(DRF_env, 2, min), DRF_mean = apply(DRF_env, 2, mean), DRF_max = apply(DRF_env, 2, max))
```

    ##                DRF_min  DRF_mean DRF_max
    ## hydroperiod       0.00    0.3600     1.0
    ## area              1.15 2397.0658 66803.5
    ## depth             0.10    0.8704     4.0
    ## canopy_cover      0.00    1.4800     3.0
    ## nvt               2.00    2.9400     3.0
    ## dist_to_forest    0.00   13.1200   209.0

 

#### SSF

``` r
data.frame(SSF_min = apply(SSF_clim, 2, min), SSF_mean = apply(SSF_clim, 2, mean), SSF_max = apply(SSF_clim, 2, max))
```

    ##             SSF_min   SSF_mean SSF_max
    ## temp_Season    1890 2128.65217    2598
    ## range_temp      182  188.97826     197
    ## total_prec     1174 1263.69565    1435
    ## prec_season      41   67.06522      78

``` r
data.frame(SSF_min = apply(SSF_env, 2, min), SSF_mean = apply(SSF_env, 2, mean), SSF_max = apply(SSF_env, 2, max))
```

    ##                SSF_min    SSF_mean SSF_max
    ## hydroperiod       0.00   0.4782609     1.0
    ## area              4.00 695.7586957 10000.0
    ## depth             0.05   0.6728261     2.4
    ## canopy_cover      0.00   0.3478261     3.0
    ## nvt               1.00   2.0000000     3.0
    ## dist_to_forest    0.00 177.5652174   834.0

   

### Small Extent

#### Santa Fé do Sul

``` r
data.frame(ST_min = apply(ST_env, 2, min),
  ST_mean = apply(ST_env, 2, mean),
  ST_max = apply(ST_env, 2, max))
```

    ##                ST_min   ST_mean ST_max
    ## hydroperiod       0.0   0.37500      1
    ## area             50.0 270.31250    960
    ## depth             0.1   0.86625      2
    ## nvt               1.0   1.62500      2
    ## dist_to_forest   65.0 140.37500    230

 

#### Icém

``` r
data.frame(IC_min = apply(IC_env, 2, min),
  IC_mean = apply(IC_env, 2, mean),
  IC_max = apply(IC_env, 2, max))
```

    ##                IC_min     IC_mean IC_max
    ## hydroperiod       0.0   0.5833333    1.0
    ## area             60.0 398.4166667  800.0
    ## depth             0.1   0.4225000    0.7
    ## nvt               1.0   1.9166667    3.0
    ## dist_to_forest   10.0 223.9166667  551.0

 

#### Nova Itapirema

``` r
data.frame(NI_min = apply(NI_env, 2, min),
  NI_mean = apply(NI_env, 2, mean),
  NI_max = apply(NI_env, 2, max))
```

    ##                NI_min NI_mean NI_max
    ## hydroperiod      0.00   0.625    1.0
    ## area            32.00 338.550  800.0
    ## depth            0.27   0.665    1.5
    ## nvt              1.00   1.875    2.0
    ## dist_to_forest  12.00 422.500  834.0

 

#### Morro do Diabo

``` r
data.frame(MD_min = apply(MD_env, 2, min),
  MD_mean = apply(MD_env, 2, mean),
  MD_max = apply(MD_env, 2, max))
```

    ##                MD_min   MD_mean  MD_max
    ## hydroperiod       0.0    0.3750     1.0
    ## area            100.0 2300.2500 10000.0
    ## depth             0.3    0.8975     2.1
    ## canopy_cover      0.0    0.7500     3.0
    ## nvt               1.0    2.2500     3.0
    ## dist_to_forest    0.0   17.8750    36.0

 

#### Jataí

``` r
data.frame(JA_min = apply(JA_env, 2, min),
  JA_mean = apply(JA_env, 2, mean),
  JA_max = apply(JA_env, 2, max))
```

    ##                JA_min JA_mean JA_max
    ## hydroperiod      0.00   0.400    1.0
    ## area             4.00 395.100 1833.0
    ## depth            0.05   0.645    2.4
    ## canopy_cover     0.00   1.000    3.0
    ## nvt              1.00   2.300    3.0
    ## dist_to_forest   0.00  83.500  214.0

 

#### Ubatuba

``` r
data.frame(UBA_min = apply(UBA_env, 2, min),
  UBA_mean = apply(UBA_env, 2, mean),
  UBA_max = apply(UBA_env, 2, max))
```

    ##                UBA_min     UBA_mean UBA_max
    ## hydroperiod       0.00    0.4347826     1.0
    ## area              1.15 3530.8978261 66803.5
    ## depth             0.10    0.6217391     2.0
    ## canopy_cover      0.00    1.9565217     3.0
    ## dist_to_forest    0.00    5.6521739    90.0

 

#### Bertioga

``` r
data.frame(BER_min = apply(BER_env, 2, min),
  BER_mean = apply(BER_env, 2, mean),
  BER_max = apply(BER_env, 2, max))
```

    ##                BER_min    BER_mean  BER_max
    ## hydroperiod       0.00    0.250000     1.00
    ## area             70.65 1935.707500 12363.75
    ## depth             0.30    1.066667     2.00
    ## canopy_cover      0.00    1.083333     3.00
    ## dist_to_forest    0.00   43.833333   209.00

 

#### Itanhaém

``` r
data.frame(ITA_min = apply(ITA_env, 2, min),
  ITA_mean = apply(ITA_env, 2, mean),
  ITA_max = apply(ITA_env, 2, max))
```

    ##              ITA_min     ITA_mean ITA_max
    ## hydroperiod     0.00    0.3333333    1.00
    ## area           53.44 1027.6100000 4592.25
    ## depth           0.25    1.0946667    4.00
    ## canopy_cover    0.00    1.0666667    3.00
    ## nvt             2.00    2.8000000    3.00

         

## Coefficient of Variation of predictor variables

### Intermediate Extent

``` r
Int_Clim_cv <- data.frame(DRF_cv = apply(DRF_clim, 2, coef_var), SSF_cv = apply(SSF_clim, 2, coef_var))
Int_Env_cv <- data.frame(DRF_cv = apply(DRF_env, 2, coef_var), SSF_cv = apply(SSF_env, 2, coef_var))

Int_Clim_cv
```

    ##                 DRF_cv     SSF_cv
    ## temp_Season 0.05215230 0.10576190
    ## range_temp  0.05661376 0.02396501
    ## total_prec  0.09938354 0.06859710
    ## prec_season 0.06739787 0.18438093

``` r
Int_Env_cv
```

    ##                    DRF_cv    SSF_cv
    ## hydroperiod    1.34687006 1.0560073
    ## area           3.98121316 2.2115741
    ## depth          0.80900349 0.8667245
    ## canopy_cover   0.84413943 2.7252676
    ## nvt            0.08159794 0.3651484
    ## dist_to_forest 3.22261540 1.3173039

``` r
apply(Int_Clim_cv, 2, mean)
```

    ##     DRF_cv     SSF_cv 
    ## 0.06888687 0.09567624

``` r
apply(Int_Env_cv, 2, mean)
```

    ##   DRF_cv   SSF_cv 
    ## 1.714240 1.423671

   

### Small Extent

``` r
UBA_cv <- c(apply(UBA_env, 2, coef_var), nvt = NA)
BER_cv <- c(apply(BER_env, 2, coef_var), nvt = NA)
ITA_cv <- c(apply(ITA_env, 2, coef_var), dist_to_forest = NA)
ST_cv <- c(apply(ST_env, 2, coef_var), canopy_cover = NA)
IC_cv <- c(apply(IC_env, 2, coef_var), canopy_cover = NA)
NI_cv <- c(apply(NI_env, 2, coef_var), canopy_cover = NA)
MD_cv <- c(apply(MD_env, 2, coef_var))
JA_cv <- c(apply(JA_env, 2, coef_var))

UBA_cv <- UBA_cv[order(names(UBA_cv))]
BER_cv <- BER_cv[order(names(BER_cv))]
ITA_cv <- ITA_cv[order(names(ITA_cv))]
ST_cv <- ST_cv[order(names(ST_cv))]
IC_cv <- IC_cv[order(names(IC_cv))]
NI_cv <- NI_cv[order(names(NI_cv))]
MD_cv <- MD_cv[order(names(MD_cv))]
JA_cv <- JA_cv[order(names(JA_cv))]

Sm_Env_cv <- data.frame(UBA_cv,BER_cv,ITA_cv,ST_cv,IC_cv,NI_cv,MD_cv,JA_cv)

Sm_Env_cv
```

    ##                   UBA_cv    BER_cv    ITA_cv     ST_cv     IC_cv     NI_cv
    ## area           3.9205225 1.9703002 1.3393552 1.2354847 0.6800164 0.9218279
    ## canopy_cover   0.6255682 1.1447191 1.0310471        NA        NA        NA
    ## depth          0.8307725 0.4633208 0.8837021 0.8303653 0.4495484 0.5676852
    ## dist_to_forest 3.5715239 1.7326073        NA 0.4482230 0.8530787 0.9661925
    ## hydroperiod    1.1658005 1.8090681 1.4638501 1.3801311 0.8827348 0.8280787
    ## nvt                   NA        NA 0.1478712 0.3184918 0.4697409 0.1885618
    ##                    MD_cv     JA_cv
    ## area           1.4337193 1.4884257
    ## canopy_cover   1.8516402 1.4142136
    ## depth          0.7715253 1.2083799
    ## dist_to_forest 0.8586156 0.9410876
    ## hydroperiod    1.3801311 1.2909944
    ## nvt            0.3142697 0.3579446

``` r
apply(na.omit(Sm_Env_cv), 2, mean)
```

    ##    UBA_cv    BER_cv    ITA_cv     ST_cv     IC_cv     NI_cv     MD_cv     JA_cv 
    ## 1.9723652 1.4142297 1.2289691 1.1486604 0.6707665 0.7725306 1.1951252 1.3292667
