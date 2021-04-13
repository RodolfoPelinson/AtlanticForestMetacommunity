Descriptive Statistics
================
Rodolfo Pelinson
13/04/2021

These are just means and variability of the environmental predictors
used both in the EMS and Variation partitioning analyses.

First we have to load the separated data matrices we need. The data
matrices are prepared sourcing the “Loading\_data.R” file in the
Auxiliary Scripts folder.

``` r
library(AtlanticForestMetacommunity)
source("Loading_data.R")
```

### Number of Species

#### Broad Extent

``` r
ncol(Broad_pa_orig)
```

    ## [1] 57

   
 

#### Intermediate Extent

DRF

``` r
ncol(DRF_pa_orig)
```

    ## [1] 28

   

SSF

``` r
ncol(SSF_pa_orig)
```

    ## [1] 34

   

Species in commom SSF-DRF

``` r
sp_comm <- na.omit(match(colnames(DRF_pa_orig), colnames(SSF_pa_orig)))
length(sp_comm)
```

    ## [1] 5

``` r
colnames(SSF_pa_orig)[sp_comm]
```

    ## [1] "Bfa" "Dmi" "Llt" "Pcu" "Ror"

       

#### Small Extent

DRF - Ubatuba

``` r
ncol(UBA_pa_orig)
```

    ## [1] 23

   

DRF - Bertioga

``` r
ncol(BER_pa_orig)
```

    ## [1] 17

   

DRF - Itanhaém

``` r
ncol(ITA_pa_orig)
```

    ## [1] 21

   

Species in commom UBA-BER

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

   

Species in commom UBA-ITA

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

   

Species in commom ITA-BER

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

       

SSF - Santa Fé do Sul

``` r
ncol(ST_pa_orig)
```

    ## [1] 14

   

SSF - Icém

``` r
ncol(IC_pa_orig)
```

    ## [1] 22

   

SSF - Nova Itapirema

``` r
ncol(NI_pa_orig)
```

    ## [1] 21

   

SSF - Morro do Diabo

``` r
ncol(MD_pa_orig)
```

    ## [1] 17

   

SSF - Jataí

``` r
ncol(JA_pa_orig)
```

    ## [1] 15

   

Species in commom ST-IC

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

   

Species in commom ST-NI

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

   

Species in commom ST-MD

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

   

Species in commom ST-JA

``` r
sp_comm <- na.omit(match(colnames(ST_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 7

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dna" "Lfs" "Lpo" "Pcu" "Sfs" "Ssi"

   

Species in commom IC-NI

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

   

Species in commom IC-MD

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

   

Species in commom IC-JA

``` r
sp_comm <- na.omit(match(colnames(IC_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 9

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dmi" "Dna" "Lfs" "Llb" "Lpo" "Pcu" "Sfs" "Ssi"

   

Species in commom NI-MD

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

   

Species in commom NI-MD

``` r
sp_comm <- na.omit(match(colnames(NI_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 9

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bal" "Dmi" "Dna" "Lfs" "Llb" "Lpo" "Pcu" "Sfs" "Ssi"

   

Species in commom JA-MD

``` r
sp_comm <- na.omit(match(colnames(MD_pa_orig), colnames(JA_pa_orig)))
length(sp_comm)
```

    ## [1] 8

``` r
colnames(JA_pa_orig)[sp_comm]
```

    ## [1] "Bfa" "Dmi" "Dna" "Lfs" "Lpo" "Pcu" "Sfs" "Ssi"

               

### Means of continuous variables

Intermediate Extent

``` r
data.frame(DRF_mean = apply(DRF_clim, 2, mean), SSF_mean = apply(SSF_clim, 2, mean))
```

    ##                DRF_mean   SSF_mean
    ## temp_Season     2307.70 2128.65217
    ## range_temp       161.32  188.97826
    ## total_prec      2283.88 1263.69565
    ## prec_season       41.80   67.06522
    ## prec_wet_quart   878.72  615.91304
    ## prec_dry_quart   276.76   80.65217

``` r
data.frame(DRF_mean = apply(DRF_env[,3:5], 2, mean), SSF_mean = apply(SSF_env[,3:5], 2, mean))
```

    ##        DRF_mean    SSF_mean
    ## area  2397.0658 695.7586957
    ## depth    0.8704   0.6728261
    ## nvt      2.9400   2.0000000

   

Small Extent

``` r
data.frame(UBA_mean = apply(UBA_env[,3:5], 2, mean),
           BER_mean = apply(BER_env[,3:5], 2, mean),
           ITA_mean = apply(ITA_env[,3:5], 2, mean),
           ST_mean = apply(ST_env[,3:5], 2, mean),
           IC_mean = apply(IC_env[,3:5], 2, mean),
           NI_mean = apply(NI_env[,3:5], 2, mean),
           MD_mean = apply(MD_env[,3:5], 2, mean),
           JA_mean = apply(JA_env[,3:5], 2, mean))
```

    ##           UBA_mean    BER_mean    ITA_mean   ST_mean    IC_mean NI_mean
    ## area  3530.8978261 1935.707500 1027.610000 270.31250 398.416667 338.550
    ## depth    0.6217391    1.066667    1.094667   0.86625   0.422500   0.665
    ## nvt      3.0000000    3.000000    2.800000   1.62500   1.916667   1.875
    ##         MD_mean JA_mean
    ## area  2300.2500 395.100
    ## depth    0.8975   0.645
    ## nvt      2.2500   2.300

       

### Coefficient of Variation of continuous variables

Intermediate Extent

``` r
Int_Clim_cv <- data.frame(DRF_cv = apply(DRF_clim, 2, coef_var), SSF_cv = apply(SSF_clim, 2, coef_var))
Int_Env_cv <- data.frame(DRF_cv = apply(DRF_env[,3:5], 2, coef_var), SSF_cv = apply(SSF_env[,3:5], 2, coef_var))

Int_Clim_cv
```

    ##                    DRF_cv     SSF_cv
    ## temp_Season    0.05233907 0.10576190
    ## range_temp     0.05608330 0.02396501
    ## total_prec     0.09952700 0.06859710
    ## prec_season    0.06714587 0.18438093
    ## prec_wet_quart 0.06420671 0.14778671
    ## prec_dry_quart 0.14612379 0.38828156

``` r
Int_Env_cv
```

    ##           DRF_cv    SSF_cv
    ## area  3.98121316 2.2115741
    ## depth 0.80900349 0.8667245
    ## nvt   0.08159794 0.3651484

``` r
apply(Int_Clim_cv, 2, mean)
```

    ##     DRF_cv     SSF_cv 
    ## 0.08090429 0.15312887

``` r
apply(Int_Env_cv, 2, mean)
```

    ##   DRF_cv   SSF_cv 
    ## 1.623938 1.147816

   

Small Extent

``` r
Sm_Env_cv <- data.frame(UBA_cv = apply(UBA_env[,3:5], 2, coef_var),
           BER_cv = apply(BER_env[,3:5], 2, coef_var),
           ITA_cv = apply(ITA_env[,3:5], 2, coef_var),
           ST_cv = apply(ST_env[,3:5], 2, coef_var),
           IC_cv = apply(IC_env[,3:5], 2, coef_var),
           NI_cv = apply(NI_env[,3:5], 2, coef_var),
           MD_cv = apply(MD_env[,3:5], 2, coef_var),
           JA_cv = apply(JA_env[,3:5], 2, coef_var))

Sm_Env_cv
```

    ##          UBA_cv    BER_cv    ITA_cv     ST_cv     IC_cv     NI_cv     MD_cv
    ## area  3.9205225 1.9703002 1.3393552 1.2354847 0.6800164 0.9218279 1.4337193
    ## depth 0.8307725 0.4633208 0.8837021 0.8303653 0.4495484 0.5676852 0.7715253
    ## nvt   0.0000000 0.0000000 0.1478712 0.3184918 0.4697409 0.1885618 0.3142697
    ##           JA_cv
    ## area  1.4884257
    ## depth 1.2083799
    ## nvt   0.3579446

``` r
apply(Sm_Env_cv, 2, mean)
```

    ##    UBA_cv    BER_cv    ITA_cv     ST_cv     IC_cv     NI_cv     MD_cv     JA_cv 
    ## 1.5837650 0.8112070 0.7903095 0.7947806 0.5331019 0.5593583 0.8398381 1.0182501
