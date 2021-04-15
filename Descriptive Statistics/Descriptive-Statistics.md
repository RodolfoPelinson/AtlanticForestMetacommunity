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

               

### Means, maximum and minimum of continuous variables

 

#### Intermediate Extent

  DRF

``` r
data.frame(DRF_min = apply(DRF_clim, 2, min), DRF_mean = apply(DRF_clim, 2, mean), DRF_max = apply(DRF_clim, 2, max))
```

    ##                DRF_min DRF_mean DRF_max
    ## temp_Season       2194  2307.70    2603
    ## range_temp         143   161.32     173
    ## total_prec        1911  2283.88    2747
    ## prec_season         36    41.80      48
    ## prec_wet_quart     762   878.72    1001
    ## prec_dry_quart     219   276.76     369

``` r
data.frame(DRF_min = apply(DRF_env, 2, min), DRF_mean = apply(DRF_env, 2, mean), DRF_max = apply(DRF_env, 2, max))
```

    ##              DRF_min  DRF_mean DRF_max
    ## hydroperiod     0.00    0.3600     1.0
    ## canopy_cover    0.00    0.4400     1.0
    ## area            1.15 2397.0658 66803.5
    ## depth           0.10    0.8704     4.0
    ## nvt             2.00    2.9400     3.0

 

SSF

``` r
data.frame(SSF_min = apply(SSF_clim, 2, min), SSF_mean = apply(SSF_clim, 2, mean), SSF_max = apply(SSF_clim, 2, max))
```

    ##                SSF_min   SSF_mean SSF_max
    ## temp_Season       1890 2128.65217    2598
    ## range_temp         182  188.97826     197
    ## total_prec        1174 1263.69565    1435
    ## prec_season         41   67.06522      78
    ## prec_wet_quart     444  615.91304     717
    ## prec_dry_quart      50   80.65217     148

``` r
data.frame(SSF_min = apply(SSF_env, 2, min), SSF_mean = apply(SSF_env, 2, mean), SSF_max = apply(SSF_env, 2, max))
```

    ##              SSF_min    SSF_mean SSF_max
    ## hydroperiod     0.00   0.4782609     1.0
    ## canopy_cover    0.00   0.1086957     1.0
    ## area            4.00 695.7586957 10000.0
    ## depth           0.05   0.6728261     2.4
    ## nvt             1.00   2.0000000     3.0

   

#### Small Extent

Ubatuba

``` r
data.frame(UBA_min = apply(UBA_env, 2, min),
  UBA_mean = apply(UBA_env, 2, mean),
  UBA_max = apply(UBA_env, 2, max))
```

    ##              UBA_min     UBA_mean UBA_max
    ## hydroperiod     0.00    0.4347826     1.0
    ## canopy_cover    0.00    0.6956522     1.0
    ## area            1.15 3530.8978261 66803.5
    ## depth           0.10    0.6217391     2.0
    ## nvt             3.00    3.0000000     3.0

 

Bertioga

``` r
data.frame(BER_min = apply(BER_env, 2, min),
  BER_mean = apply(BER_env, 2, mean),
  BER_max = apply(BER_env, 2, max))
```

    ##              BER_min    BER_mean  BER_max
    ## hydroperiod     0.00    0.250000     1.00
    ## canopy_cover    0.00    0.250000     1.00
    ## area           70.65 1935.707500 12363.75
    ## depth           0.30    1.066667     2.00
    ## nvt             3.00    3.000000     3.00

 

Itanhaém

``` r
data.frame(ITA_min = apply(ITA_env, 2, min),
  ITA_mean = apply(ITA_env, 2, mean),
  ITA_max = apply(ITA_env, 2, max))
```

    ##              ITA_min     ITA_mean ITA_max
    ## hydroperiod     0.00    0.3333333    1.00
    ## canopy_cover    0.00    0.2000000    1.00
    ## area           53.44 1027.6100000 4592.25
    ## depth           0.25    1.0946667    4.00
    ## nvt             2.00    2.8000000    3.00

 

Santa Fé do Sul

``` r
data.frame(ST_min = apply(ST_env, 2, min),
  ST_mean = apply(ST_env, 2, mean),
  ST_max = apply(ST_env, 2, max))
```

    ##              ST_min   ST_mean ST_max
    ## hydroperiod     0.0   0.37500      1
    ## canopy_cover    0.0   0.00000      0
    ## area           50.0 270.31250    960
    ## depth           0.1   0.86625      2
    ## nvt             1.0   1.62500      2

 

Icém

``` r
data.frame(IC_min = apply(IC_env, 2, min),
  IC_mean = apply(IC_env, 2, mean),
  IC_max = apply(IC_env, 2, max))
```

    ##              IC_min     IC_mean IC_max
    ## hydroperiod     0.0   0.5833333    1.0
    ## canopy_cover    0.0   0.0000000    0.0
    ## area           60.0 398.4166667  800.0
    ## depth           0.1   0.4225000    0.7
    ## nvt             1.0   1.9166667    3.0

 

Nova Itapirema

``` r
data.frame(NI_min = apply(NI_env, 2, min),
  NI_mean = apply(NI_env, 2, mean),
  NI_max = apply(NI_env, 2, max))
```

    ##              NI_min NI_mean NI_max
    ## hydroperiod    0.00   0.625    1.0
    ## canopy_cover   0.00   0.000    0.0
    ## area          32.00 338.550  800.0
    ## depth          0.27   0.665    1.5
    ## nvt            1.00   1.875    2.0

 

``` r
data.frame(MD_min = apply(MD_env, 2, min),
  MD_mean = apply(MD_env, 2, mean),
  MD_max = apply(MD_env, 2, max))
```

    ##              MD_min   MD_mean  MD_max
    ## hydroperiod     0.0    0.3750     1.0
    ## canopy_cover    0.0    0.2500     1.0
    ## area          100.0 2300.2500 10000.0
    ## depth           0.3    0.8975     2.1
    ## nvt             1.0    2.2500     3.0

 

``` r
data.frame(JA_min = apply(JA_env, 2, min),
  JA_mean = apply(JA_env, 2, mean),
  JA_max = apply(JA_env, 2, max))
```

    ##              JA_min JA_mean JA_max
    ## hydroperiod    0.00   0.400    1.0
    ## canopy_cover   0.00   0.300    1.0
    ## area           4.00 395.100 1833.0
    ## depth          0.05   0.645    2.4
    ## nvt            1.00   2.300    3.0

       

### Coefficient of Variation of continuous variables

Intermediate Extent

``` r
Int_Clim_cv <- data.frame(DRF_cv = apply(DRF_clim, 2, coef_var), SSF_cv = apply(SSF_clim, 2, coef_var))
Int_Env_cv <- data.frame(DRF_cv = apply(DRF_env, 2, coef_var), SSF_cv = apply(SSF_env, 2, coef_var))

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

    ##                  DRF_cv    SSF_cv
    ## hydroperiod  1.34687006 1.0560073
    ## canopy_cover 1.13960576 2.8952068
    ## area         3.98121316 2.2115741
    ## depth        0.80900349 0.8667245
    ## nvt          0.08159794 0.3651484

``` r
apply(Int_Clim_cv, 2, mean)
```

    ##     DRF_cv     SSF_cv 
    ## 0.08090429 0.15312887

``` r
apply(Int_Env_cv, 2, mean)
```

    ##   DRF_cv   SSF_cv 
    ## 1.471658 1.478932

   

Small Extent

``` r
Sm_Env_cv <- data.frame(UBA_cv = apply(UBA_env, 2, coef_var),
           BER_cv = apply(BER_env, 2, coef_var),
           ITA_cv = apply(ITA_env, 2, coef_var),
           ST_cv = apply(ST_env, 2, coef_var),
           IC_cv = apply(IC_env, 2, coef_var),
           NI_cv = apply(NI_env, 2, coef_var),
           MD_cv = apply(MD_env, 2, coef_var),
           JA_cv = apply(JA_env, 2, coef_var))

Sm_Env_cv
```

    ##                 UBA_cv    BER_cv    ITA_cv     ST_cv     IC_cv     NI_cv
    ## hydroperiod  1.1658005 1.8090681 1.4638501 1.3801311 0.8827348 0.8280787
    ## canopy_cover 0.6763035 1.8090681 2.0701967       NaN       NaN       NaN
    ## area         3.9205225 1.9703002 1.3393552 1.2354847 0.6800164 0.9218279
    ## depth        0.8307725 0.4633208 0.8837021 0.8303653 0.4495484 0.5676852
    ## nvt          0.0000000 0.0000000 0.1478712 0.3184918 0.4697409 0.1885618
    ##                  MD_cv     JA_cv
    ## hydroperiod  1.3801311 1.2909944
    ## canopy_cover 1.8516402 1.6101530
    ## area         1.4337193 1.4884257
    ## depth        0.7715253 1.2083799
    ## nvt          0.3142697 0.3579446

``` r
apply(Sm_Env_cv, 2, mean)
```

    ##   UBA_cv   BER_cv   ITA_cv    ST_cv    IC_cv    NI_cv    MD_cv    JA_cv 
    ## 1.318680 1.210351 1.180995      NaN      NaN      NaN 1.150257 1.191180
