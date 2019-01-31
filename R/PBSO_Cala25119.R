## Code credits########################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection
## fitting PO PA and Integrated models the data
## 19/08/2016
##
## GBIF and wolrdclim data reading capability by Diego J. Lizcano. 
## 12/07/2018
##
## Editing structures, comments and adaptation for Coyote by Gabriel Andrade-Ponce
##

######### Data required for the function ##########################################################

## All the data required for this function is stored in a file data.rda
## s.intensity - raster with background covatiates that effect occupancy
## s.detection - raster with background covariates that effect detection
## pb.intensity - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
## pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
##
## y.so - matrix of detection/non detection of the PA surveys (Camera trap data)
## so.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
## so.detection - matrix with covariates that effect detection in the locations of PA survey sites in each survey

###Canis latrans-V2. Gabriel Andrade Ponce###

########Required libraries and functions ######
library(raster)
library(fields)
library(mvtnorm)
library(matrixStats)
library(rgbif) # get GBIF data
library(maptools) 
source("R/functions.r") # PoANDso.fit functions of Koshkina et al. (2017) supplementary material


####### Ready to start #############

#################################### PO data #######################
# 1. Preparing PO the data ========================================================
#1.1 Get GBIF Data per species========================================================

#sp_in_GBIF <- occ_search(scientificName = "Canis latrans", return='data') to search in GBIF data base
sp_in_GBIF <- read.csv("Cala_gbif.csv") # read prescence points file

sp_points <- SpatialPoints(cbind(sp_in_GBIF$longitude, 
                                 sp_in_GBIF$latitude))

# put georeference
proj4string(sp_points) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
#centroid <- c(mean(coordinates(sp_points)[,1]),mean(coordinates(sp_points)[,2])) 

#1.2 Get altitude and worldclim Data for the species======================================================== 

# La M del coyote es muy grande por lo que los datos se bajaron manualmente
altitude <- raster("altitude.asc")#tomada como M

#bio_vars <- getData('worldclim', var='bio' , res=0.5, lon=centroid[1], lat=centroid[2]) Me toco descargarlas manualmente, por que la M es muy grande
bio_vars <- stack("bio_vars.tif") #ya estan cortadas y ajustadas a la M


#1.3 cut by area of interest ========================================================   
## in this case mask and extent by M           
# Como los raster ya pasaron el proceso de cut y mask no realiza de nuevo
# read M for Canis latrans
# La altitud funciona como M "M <- raster('Mcala.asc')" # equivale al M
# Prec <- readShapePoly('shp/Unidad_analisis.shp')
#proj4string(M) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
### cut Bioclim## La capa de Bioclim ya esta ajustada
#bio_vars_crop <- crop(bio_vars, extent(altitude)) #crop by a raster. M
#extent(bio_vars_crop) <- alignExtent(bio_vars_crop,altitude) #Aligne extent
#bio_vars_cut <- mask(bio_vars_crop, altitude)# mask


### cut altutude# la capa de altitud ya esta ajustada
#alti_crop <- crop(alt_mosaic, extent(Prec))
#alti_crop_r <- resample(alti_crop, bio_vars_cut, method='bilinear')
#extent(alti_crop_r) <- alignExtent(alti_crop_r, Prec)
#alti_cut <- mask(alti_crop_r, Prec)
#names(alti_cut) <- "alti" # put name


plot(bio_vars[[2]])
plot(sp_points, add=T)

#1.3.1 read spatial covariates and standarize using scale ========================================================   
s.intensity <- stack(scale(altitude, center=TRUE, scale=TRUE),
                     scale(bio_vars[[12]], center=TRUE, scale=TRUE), # Annual Precipitation
                     scale(bio_vars[[1]], center=TRUE, scale=TRUE),  # Annual Mean Temperature
                     scale(bio_vars[[3]], center=TRUE, scale=TRUE))#, # Isothermality (BIO2/BIO7) (* 100)
# scale(bio_vars_cut[[4]], center=TRUE, scale=TRUE), 
# scale(bio_vars_cut[[5]], center=TRUE, scale=TRUE),
# scale(bio_vars_cut[[6]], center=TRUE, scale=TRUE), 
# scale(bio_vars_cut[[7]], center=TRUE, scale=TRUE),
# scale(bio_vars_cut[[12]], center=TRUE, scale=TRUE)) # Annual Precipitation
#writeRaster(s.intensity,"s.intensity.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

#s.intensity <- stack("s.intensity.tif")

s.detection <-  stack(scale(altitude, center=TRUE, scale=TRUE),
                    scale(bio_vars[[12]], center=TRUE, scale=TRUE), # Annual Precipitation
                    scale(bio_vars[[1]], center=TRUE, scale=TRUE),  # Annual Mean temperature
                    scale(bio_vars[[3]], center=TRUE, scale=TRUE))#, # Isothermality (BIO2/BIO7) (* 100)
# scale(bio_vars_cut[[13]], center=TRUE, scale=TRUE), # BIO13 = Precipitation of Wettest Month
# scale(bio_vars_cut[[14]], center=TRUE, scale=TRUE), # BIO14 = Precipitation of Driest Month
# scale(bio_vars_cut[[15]], center=TRUE, scale=TRUE)) # BIO15 = Precipitation Seasonality (Coefficient of Variation)
#writeRaster(s.detection,"s.detection.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
#s.detection <- stack("s.detection.tif")

# plot... just to check 
plot(s.intensity[[1]]) 
plot(sp_points, add=T)

#1.3.2 Get the covariates from scaled rasters ========================================================

pb <- sp_points # coming from GBIF
pb.loc <- pb

pb.intensity=extract(s.intensity, pb.loc, method='bilinear') 
pb.detection=extract(s.detection, pb.loc, method='bilinear') 
#write.csv(pb.detection, "pb.detection.csv")
#write.csv(pb.detection, "pb.intensity.csv")

#pb.intensity <- (read.csv("pb.intensity.csv")[,-1])
#pb.detection <- (read.csv("pb.detection.csv")[,-1])

### deja solo puntos dentro de la ventana
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.intensity)
pb.detection=pb.detection[is.complete.pb,]
pb.intensity=pb.intensity[is.complete.pb,]

# 1.4 Preparing the data ========================================================

#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.intensity)), values(s.intensity))
colnames(X.back)=c("",names(s.intensity))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))
colnames(W.back)=c("",names(s.detection))
# remove all NA values
tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km 
area.back = rep((xres(s.intensity)/1000)*(yres(s.intensity)/1000), nrow(X.back))# each cell
s.area=area.back*nrow(X.back) #study area

# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.intensity))), pb.intensity)
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)

#### PO done

#################################### PA data #######################
# 2. Preparing PA the data ========================================================
#2.1 Read data, Occupancy data measured in the field ========================================================

PA <- (read.csv("data/Cala.hist3noNA.csv",sep=",")[,-1]) #sp capture hist

#Read occupancy covariates
#bait <- read.csv("data/bait_covsep8119.csv", sep = ";")
#rain <- read.csv("data/rain_cov8119.csv", sep = ";",  na.strings = "NA")
#site_cov <- read.csv("data/site_cov.sep20112018.csv", sep=";")

CovOcc <-  as.matrix(read.csv("cov.oc.csv")[,-1])
CovOccPB <- as.matrix(read.csv("cov.det.csv")[,-1]) 


#2.2 Prepare PA data ========================================================
#These 3 files must have the same nùmber of rows
so.occupancy <- CovOcc #Occupancy covariates
so.detection <- CovOccPB # Detection covariates
y.so <- PA # sp capture matrix


print("removing NA")
is.complete=complete.cases(so.occupancy)&complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*39.4


 
#2.3 Transform data to Occ and det Matrix ========================================================

J.so=ncol(y.so)
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)




W.so = array(dim=c(nrow(as.matrix(so.detection)), J.so, 5))# El nùmero de matrices depende del nùmero de covariables
W.so[,,1] = 1
W.so[,,2] = matrix(rep(so.detection$s.detection.1,J.so),ncol = J.so)# use rep to site covariates
W.so[,,3] = matrix(rep(so.detection$s.detection.2,J.so),ncol = J.so)
W.so[,,4] = matrix(rep(so.detection$s.detection.3,J.so),ncol = J.so)
W.so[,,5] = matrix(rep(so.detection$s.detection.4,J.so),ncol = J.so)
# You just can use W.so[,,2] = so.detection, when your detection covariate change in time

#################################### Fit functions #######################
# 3. Analising the data ========================================================

# 3.1 Analyzing Presence-Only (PO) data ========================================================
pb.fit=pb.ipp(X.po, W.po,X.back, W.back)

# 3.2 Analyzing presence-absence (PA) data ========================================================
so.fit=so.model(X.so,W.so,y.so)

# 3.3 Analyzing presence-only data AND presence-absence data ========================================================
poANDso.fit=pbso.integrated(X.po, W.po,X.back, W.back,X.so,W.so,y.so)
