## LOADING DATA
############################

library("vegan")
data = read.csv("C:/Users/rodol/OneDrive/Trabalho/Papers/Analysis/AtlanticForestMetacommunity/Data csv/Data.csv",header=TRUE,row.names = 1, encoding = "UTF-8")
#data(Tadpoles_Atlantic_Forest)
#data <- Tadpoles_Atlantic_Forest
##################################################################################################################
###########################################    SPECIES DATA    ###################################################
##################################################################################################################

Broad=data[2:nrow(data),1:57]

for(i in 1:ncol(Broad)){
  Broad[,i]<- as.numeric(Broad[,i])
}

SSFoc<-data[which(data$ecoregion == "SSF"),1:57]
for(i in 1:ncol(SSFoc)){
  SSFoc[,i]<- as.numeric(SSFoc[,i])
}

DRFoc=data[which(data$ecoregion == "DRF"),1:57]
for(i in 1:ncol(DRFoc)){
  DRFoc[,i]<- as.numeric(DRFoc[,i])
}


SToc=(SSFoc[1:8,])
ICoc=(SSFoc[9:20,])
NIoc=(SSFoc[21:28,])
MDoc=(SSFoc[29:36,])
JAoc=(SSFoc[37:46,])


DRFoc
UBAoc=(DRFoc[c(10,11,12,13,14,15,16,17,18,19,38,39,40,41,42,43,44,45,46,47,48,49,50),])
BERoc=(DRFoc[26:37,])
ITAoc=(DRFoc[c(1,2,3,4,5,6,7,8,9,20,21,22,23,24,25),])


#Transforming data to presence absence
Broad_pa <- decostand(Broad, method = "pa")

SSF_pa <- decostand(SSFoc, method = "pa")
IC_pa <- decostand(ICoc, method = "pa")
NI_pa <- decostand(NIoc, method = "pa")
ST_pa <- decostand(SToc, method = "pa")
MD_pa <- decostand(MDoc, method = "pa")
JA_pa <- decostand(JAoc, method = "pa")

DRF_pa <- decostand(DRFoc, method = "pa")
UBA_pa <- decostand(UBAoc, method = "pa")
BER_pa <- decostand(BERoc, method = "pa")
ITA_pa <- decostand(ITAoc, method = "pa")


Broad_pa_orig <- Broad_pa[,which (colSums(Broad_pa)>0)]

SSF_pa_orig <- SSF_pa[,which (colSums(SSF_pa)>0)]
ST_pa_orig <- ST_pa[,which (colSums(ST_pa)>0)]
IC_pa_orig <- IC_pa[,which (colSums(IC_pa)>0)]
NI_pa_orig <- NI_pa[,which (colSums(NI_pa)>0)]
MD_pa_orig <- MD_pa[,which (colSums(MD_pa)>0)]
JA_pa_orig <- JA_pa[,which (colSums(JA_pa)>0)]

DRF_pa_orig <- DRF_pa[,which (colSums(DRF_pa)>0)]
UBA_pa_orig <- UBA_pa[,which (colSums(UBA_pa)>0)]
BER_pa_orig <- BER_pa[,which (colSums(BER_pa)>0)]
ITA_pa_orig <- ITA_pa[,which (colSums(ITA_pa)>0)]



#Removing singletons
Broad_pa=Broad_pa_orig[,which (colSums(Broad_pa_orig)>1)]

SSF_pa=SSF_pa_orig[,which (colSums(SSF_pa_orig)>1)]
ST_pa=ST_pa_orig[,which (colSums(ST_pa_orig)>1)]
IC_pa=IC_pa_orig[,which (colSums(IC_pa_orig)>1)]
NI_pa=NI_pa_orig[,which (colSums(NI_pa_orig)>1)]
MD_pa=MD_pa_orig[,which (colSums(MD_pa_orig)>1)]
JA_pa=JA_pa_orig[,which (colSums(JA_pa_orig)>1)]

DRF_pa=DRF_pa_orig[,which (colSums(DRF_pa_orig)>1)]
UBA_pa=UBA_pa_orig[,which (colSums(UBA_pa_orig)>1)]
BER_pa=BER_pa_orig[,which (colSums(BER_pa_orig)>1)]
ITA_pa=ITA_pa_orig[,which (colSums(ITA_pa_orig)>1)]



##################################################################################################################
########################################    ENVIRONMENTAL DATA    ################################################
##################################################################################################################
Broad_locality <- data$locality

Broadenv<-data.frame(hydroperiod = data$hydroperiod[2:nrow(data)],
                    canopy_cover = data$canopy_cover[2:nrow(data)],
                    area = data$area[2:nrow(data)],
                    depth = data$depth[2:nrow(data)],
                    #nvt = data$nvt[2:nrow(data)],
                    arbust_veg = data$arbust_veg[2:nrow(data)],
                    arbor_veg = data$arbor_veg[2:nrow(data)],
                    dist_to_forest = data$dist_to_forest[2:nrow(data)],
                    ecoregion = as.factor(data$ecoregion[2:nrow(data)]))

rownames(Broadenv)<- rownames(data)[2:nrow(data)]
for(i in 1:ncol(Broadenv)){
  Broadenv[,i]<- as.numeric(Broadenv[,i])
}

SSF_locality <- data$locality[2:47]

#SSFenv<-data.frame(hydroperiod = data$hydroperiod,
#                   canopy_cover = data$canopy_cover,
#                   area = data$area,
#                   depth = data$depth,
#                   nvt = data$nvt)[2:47,]

SSFenv<-Broadenv[1:46,-dim(Broadenv)[2]]

rownames(SSFenv)<- rownames(data)[2:47]

for(i in 1:ncol(SSFenv)){
  SSFenv[,i]<- as.numeric(SSFenv[,i])
}

DRF_locality <- data$locality[48:nrow(data)]

DRFenv<-Broadenv[47:nrow(Broadenv),-dim(Broadenv)[2]]

#DRFenv<-data.frame(hydroperiod = data$hydroperiod,
#                   canopy_cover = data$canopy_cover,
#                   area = data$area,
#                   depth = data$depth,
#                   nvt = data$nvt)[48:nrow(data),]
rownames(DRFenv)<- rownames(data)[48:nrow(data)]
for(i in 1:ncol(DRFenv)){
  DRFenv[,i]<- as.numeric(DRFenv[,i])
}


Broad_env <- Broadenv

SSF_env <- SSFenv
IC_env <- SSFenv[9:20,]
NI_env <- SSFenv[21:28,]
ST_env <- SSFenv[1:8,]
MD_env <- SSFenv[29:36,]
JA_env <- SSFenv[37:46,]


DRF_env <- DRFenv
UBA_env=(DRFenv[c(10,11,12,13,14,15,16,17,18,19,38,39,40,41,42,43,44,45,46,47,48,49,50),])
BER_env=(DRFenv[26:37,])
ITA_env=(DRFenv[c(1,2,3,4,5,6,7,8,9,20,21,22,23,24,25),])


Broad_env <- Broad_env[,which(apply(na.omit(Broad_env),2,sd) != 0)]

SSF_env <- SSF_env[,which(apply(na.omit(SSF_env),2,sd) != 0)]
IC_env <- IC_env[,which(apply(na.omit(IC_env),2,sd) != 0)]
NI_env <- NI_env[,which(apply(na.omit(NI_env),2,sd) != 0)]
ST_env <- ST_env[,which(apply(na.omit(ST_env),2,sd) != 0)]
MD_env <- MD_env[,which(apply(na.omit(MD_env),2,sd) != 0)]
JA_env <- JA_env[,which(apply(na.omit(JA_env),2,sd) != 0)]


DRF_env <- DRF_env[,which(apply(na.omit(DRF_env),2,sd) != 0)]
UBA_env <- UBA_env[,which(apply(na.omit(UBA_env),2,sd) != 0)]
BER_env <- BER_env[,which(apply(na.omit(BER_env),2,sd) != 0)]
ITA_env <- ITA_env[,which(apply(na.omit(ITA_env),2,sd) != 0)]




Broad_env_st<- decostand(Broad_env, method = "standardize",na.rm=TRUE)

SSF_env_st<- decostand(SSF_env, method = "standardize",na.rm=TRUE)
IC_env_st <- decostand(IC_env, method = "standardize",na.rm=TRUE)
NI_env_st <- decostand(NI_env, method = "standardize",na.rm=TRUE)
ST_env_st <- decostand(ST_env, method = "standardize",na.rm=TRUE)
MD_env_st <- decostand(MD_env, method = "standardize",na.rm=TRUE)
JA_env_st <- decostand(JA_env, method = "standardize",na.rm=TRUE)

DRF_env_st <- decostand(DRF_env, method = "standardize",na.rm=TRUE)
UBA_env_st <- decostand(UBA_env, method = "standardize",na.rm=TRUE)
BER_env_st <- decostand(BER_env, method = "standardize",na.rm=TRUE)
ITA_env_st <- decostand(ITA_env, method = "standardize",na.rm=TRUE)







##################################################################################################################
#########################################    SPATIAL DATA    #####################################################
##################################################################################################################

Broad_coord <- data.frame(lat = data$lat[2:nrow(data)], long = data$long[2:nrow(data)])

rownames(Broad_coord) <- rownames(data[2:nrow(data),])

Broad_coord[,1] <- jitter(Broad_coord[,1], amount = 0.001)
Broad_coord[,2] <- jitter(Broad_coord[,2], amount = 0.001)

SSF_coord <- data.frame(lat = data$lat, long = data$long)[2:47,]
rownames(SSF_coord) <- rownames(data[2:47,])

IC_coord=SSF_coord[9:20,]
ST_coord=SSF_coord[1:8,]
NI_coord=SSF_coord[21:28,]
MD_coord=SSF_coord[29:36,]
JA_coord=SSF_coord[37:46,]


DRF_coord=data.frame(lat = data$lat, long = data$long)[48:nrow(data),]
rownames(DRF_coord) <- rownames(data[48:nrow(data),])

UBA_coord=(DRF_coord[c(10,11,12,13,14,15,16,17,18,19,38,39,40,41,42,43,44,45,46,47,48,49,50),])
BER_coord=(DRF_coord[26:37,])
ITA_coord=(DRF_coord[c(1,2,3,4,5,6,7,8,9,20,21,22,23,24,25),])




##################################################################################################################
###########################################    CLIMATE DATA    ###################################################
##################################################################################################################

library(dismo)
library(raster)
library(rgdal)

r <- getData("worldclim",var="bio",res=5)

#bio4 = Temperature seasonality (standard deviation *100)
#bio6 = Min temperature of coldest month
#bio7 = Temperature annual range (MAX-MIN)
#bio12 = Total (annual) precipitation
#bio15 = Precipitation seasonality (coefficient of variation)

r <- r[[c(4,7,12,15,16,17)]]

Broad_coord_plot <- SpatialPoints(Broad_coord[,c(2,1)], proj4string = r@crs)

names(r) <- c("temp_Season","range_temp", "total_prec", "prec_season", "prec_wet_quart", "prec_dry_quart")

Broad_clim <- extract(r,Broad_coord_plot)
rownames(Broad_clim) <- rownames(Broad_pa)

Broad_clim<- Broad_clim[,-c(5,6)]

Broad_clim_st<- decostand(Broad_clim[,-c(5,6)], method = "standardize")
Broad_clim_rank <- apply(Broad_clim[,-c(5,6)], 2, rank)


SSF_clim <- Broad_clim[which(Broad_env$ecoregion == 2),]
DRF_clim <- Broad_clim[which(Broad_env$ecoregion == 1),]

SSF_clim_st <- Broad_clim_st[which(Broad_env$ecoregion == 2),]
DRF_clim_st <- Broad_clim_st[which(Broad_env$ecoregion == 1),]

SSF_clim_rank <- Broad_clim_rank[which(Broad_env$ecoregion == 2),]
DRF_clim_rank <- Broad_clim_rank[which(Broad_env$ecoregion == 1),]



