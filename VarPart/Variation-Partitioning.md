Variation Partitioning
================
Rodolfo Pelinson
07/04/2021

First we have to load the separated data matrices we need. It all comes
from the Data.csv file in the Data folder. The data matrices are
prepared sourcing the “Loading\_data.R” file in the Auxiliary Scripts
folder.

``` r
source("Loading_data.R")
```

The used packages to run this analysis are:

`vegan` version 2.5-6 `ade4` version 1.7-15 `adespatial` version 0.3-8

To remove multicolinear variables from our Environmental Matrix I built
a function that removes variables with the variance inflation factor
(VIF) higher than 3. The function `VIF_selection()` simply uses the
`vif.cca()` from the `vegan` package.

``` r
Broad_env_VIF <- VIF_selection(Broad_pa,  Broad_env_st[,-6])

SSF_env_VIF <- VIF_selection(SSF_pa,  SSF_env_st)
ST_env_VIF <- VIF_selection(ST_pa,  ST_env_st)
IC_env_VIF <- VIF_selection(IC_pa,  IC_env_st)
NI_env_VIF <- VIF_selection(NI_pa,  NI_env_st)
MD_env_VIF <- VIF_selection(MD_pa,  MD_env_st)
JA_env_VIF <- VIF_selection(JA_pa,  JA_env_st)

DRF_env_VIF <- VIF_selection(DRF_pa,  DRF_env_st)
UBA_env_VIF <- VIF_selection(UBA_pa,  UBA_env_st)
BER_env_VIF <- VIF_selection(BER_pa,  BER_env_st)
ITA_env_VIF <- VIF_selection(ITA_pa,  ITA_env_st)
```

Same thing for climate variables

``` r
Broad_clim_VIF <- VIF_selection(Broad_pa,  Broad_clim_st[,-6])
SSF_clim_VIF <- VIF_selection(SSF_pa,  SSF_clim_st)
DRF_clim_VIF <- VIF_selection(DRF_pa,  DRF_clim_st)
```

Constructing matrices of spatial filters

``` r
Broad_dbMEM <- dbmem(dist(Broad_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 2.740265 
    ## Time to compute dbMEMs = 0.030000  sec

``` r
SSF_dbMEM <- dbmem(dist(SSF_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 2.608894 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
ST_dbMEM <- dbmem(dist(ST_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.006243823 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
IC_dbMEM <- dbmem(dist(IC_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.03260922 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
NI_dbMEM <- dbmem(dist(NI_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.01866376 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
MD_dbMEM <- dbmem(dist(MD_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.1459119 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
JA_dbMEM <- dbmem(dist(JA_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.03033848 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
DRF_dbMEM <- dbmem(dist(DRF_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 1.072626 
    ## Time to compute dbMEMs = 0.020000  sec

``` r
UBA_dbMEM <- dbmem(dist(UBA_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.08923115 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
BER_dbMEM <- dbmem(dist(BER_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.09885289 
    ## Time to compute dbMEMs = 0.000000  sec

``` r
ITA_dbMEM <- dbmem(dist(ITA_coord), silent = F, thresh = NULL)
```

    ## Truncation level = 0.254218 
    ## Time to compute dbMEMs = 0.000000  sec
