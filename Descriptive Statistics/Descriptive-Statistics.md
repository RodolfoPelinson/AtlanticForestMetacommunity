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
data.frame(DRF_cv = apply(DRF_clim, 2, cv), SSF_cv = apply(SSF_clim, 2, cv))
```

    ##                   DRF_cv    SSF_cv
    ## temp_Season     5.233907 10.576190
    ## range_temp      5.608330  2.396501
    ## total_prec      9.952700  6.859710
    ## prec_season     6.714587 18.438093
    ## prec_wet_quart  6.420671 14.778671
    ## prec_dry_quart 14.612379 38.828156

``` r
data.frame(DRF_cv = apply(DRF_env[,3:5], 2, cv), SSF_cv = apply(SSF_env[,3:5], 2, cv))
```

    ##           DRF_cv    SSF_cv
    ## area  398.121316 221.15741
    ## depth  80.900349  86.67245
    ## nvt     8.159794  36.51484

   

Small Extent

``` r
data.frame(UBA_cv = apply(UBA_env[,3:5], 2, cv),
           BER_cv = apply(BER_env[,3:5], 2, cv),
           ITA_cv = apply(ITA_env[,3:5], 2, cv),
           ST_cv = apply(ST_env[,3:5], 2, cv),
           IC_cv = apply(IC_env[,3:5], 2, cv),
           NI_cv = apply(NI_env[,3:5], 2, cv),
           MD_cv = apply(MD_env[,3:5], 2, cv),
           JA_cv = apply(JA_env[,3:5], 2, cv))
```

    ##          UBA_cv    BER_cv    ITA_cv     ST_cv    IC_cv    NI_cv     MD_cv
    ## area  392.05225 197.03002 133.93552 123.54847 68.00164 92.18279 143.37193
    ## depth  83.07725  46.33208  88.37021  83.03653 44.95484 56.76852  77.15253
    ## nvt     0.00000   0.00000  14.78712  31.84918 46.97409 18.85618  31.42697
    ##           JA_cv
    ## area  148.84257
    ## depth 120.83799
    ## nvt    35.79446
