###################################################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection
## fitting PO PA and Integrated models the data
## 19/08/2016
###################################################################################################
## GBIF and wolrdclim data reading capability by Diego J. Lizcano. 
## 12/07/2018
###################################################################################################
######### Data required for the function ##########################################################
###################################################################################################
## All the data required for this function is stored in a file data.rda
## s.occupancy - raster with background covatiates that effect occupancy
## s.detection - raster with background covariates that effect detection
## pb.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
## pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
##
## y.so - matrix of detection/non detection of the PA surveys (Camera trap data)
## so.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
## so.detection - matrix with covariates that effect detection in the locations of PA survey sites in each survey
######################################################################################################
###Canis latrans-V1. Gabriel Andrade Ponce###

library(raster)
library(fields)
library(mvtnorm)
library(matrixStats)
library(rgbif) # get GBIF data
library(maptools) 

# install.packages(c("raster", "fields", "mvtnorm", "matrixStats", "rgbif","maptools" ))
source("R/functions.r")

####################################
### Get GBIF Data per species
####################################
#sp_in_GBIF <- occ_search(scientificName = "Canis latrans", return='data')
sp_in_GBIF <- read.csv("Cala_gbif.csv") # read prescence points file

sp_points <- SpatialPoints(cbind(sp_in_GBIF$longitude, 
                                 sp_in_GBIF$latitude))

# put georeference
proj4string(sp_points) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
centroid <- c(mean(coordinates(sp_points)[,1]),mean(coordinates(sp_points)[,2])) 
####################################
### Get altitude and worldclim Data for the species
####################################
altitude <- raster("altitude.asc")#tomada como M

#bio_vars <- getData('worldclim', var='bio' , res=0.5, lon=centroid[1], lat=centroid[2]) Me toco descargarlas manualmente, por que la M es muy grande
bio_vars <- stack("bio_vars.tif") #ya estan cortadas y ajustadas a la M

########################################################################
## cut by area of interest   
## in this case mask and extent by M          
########################################################################

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

# Data Preparation -----------------------------------------------------------------------------
#### read spatial covariates from BIOCLIM and standarize using scale
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
                      scale(bio_vars[[1]], center=TRUE, scale=TRUE),  # Annual Precipitation
                      scale(bio_vars[[3]], center=TRUE, scale=TRUE))#, # Isothermality (BIO2/BIO7) (* 100)
# scale(bio_vars_cut[[13]], center=TRUE, scale=TRUE), # BIO13 = Precipitation of Wettest Month
# scale(bio_vars_cut[[14]], center=TRUE, scale=TRUE), # BIO14 = Precipitation of Driest Month
# scale(bio_vars_cut[[15]], center=TRUE, scale=TRUE)) # BIO15 = Precipitation Seasonality (Coefficient of Variation)
#writeRaster(s.detection,"s.detection.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
#s.detection <- stack("s.detection.tif")

# plot... just to check 
plot(s.intensity[[1]]) 
plot(sp_points, add=T)

# read pb and PA survey data. This is the only presence oportunistic data!
pb <- sp_points # coming from GBIF... no depurado


pb.loc <- pb
#### Get the covariates from rasters (stacked alti and Bioclim)
pb.intensity=extract(s.intensity, pb.loc, method='bilinear') # deja 4 puntos fuera de la ventana por fuera
pb.detection=extract(s.detection, pb.loc, method='bilinear') # deja 4 puntos fuera de la ventana por fuera
pb.intensity <- (read.csv("pb.intensity.csv")[,-1])
pb.detection <- (read.csv("pb.detection.csv")[,-1])

### deja solo puntos dentro de la ventana
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.intensity)
pb.detection=pb.detection[is.complete.pb,]
pb.intensity=pb.intensity[is.complete.pb,]

#write.csv(pb.detection, "pb.detection.csv")
#write.csv(pb.detection, "pb.intensity.csv")

################ ################ ################ ################ 
################ READ DATA, Occupancy data measured in the field
################ ################ ################ ################ 

#### estos 3 archivos deben tener el mismo numero de columnas
# PA <- PA <- as.matrix(read.csv("data/Cala.hist.csv"))
# CovOcc <- read.table('CovOcc.txt', header=T, sep='\t') # En PC_lenovo_ADP: CovOcc <- read.table('CovOcc.txt')
# CovOccPB <- read.table('CovOccPB.txt', header=T, sep='\t') # En PC_lenovo_ADP: CovOccPB <- read.table('CovOccPB.txt')
PA <- as.matrix(read.csv("data/Cala.hist.3length2.csv",sep=",")[,-1])

#bait <- read.csv("data/bait_covsep8119.csv", sep = ";")
#rain <- read.csv("data/rain_cov8119.csv", sep = ";",  na.strings = "NA")
#site_cov <- read.csv("data/site_cov.sep20112018.csv", sep=";")

CovOcc <-  (read.csv("data/cov.oc22119.csv")[,-1])# En PC_lenovo_ADP: CovOcc <- read.table('CovOcc.txt')
CovOccPB <- (read.csv("data/cov.det22119.csv")[,-1]) # En PC_lenovo_ADP: CovOccPB <- read.table('CovOccPB.txt')
# CovOccPB <- scale(CovOccPB) #

#### estos 3 archivos deben tener el mismo numero de columnas
# so.occupancy <- scale(CovOcc) #read.csv("data/so_occupancy.csv")
so.occupancy <- (CovOcc[,-c(2:4)])

#### estos 2 archivos deben tener el mismo numero de filas
# so.detection <- cbind(scale(CovOccPB[,1:4]),CovOccPB[,5]) #CovOccPB[,5] puede ser categorica, en este caso no se scala
so.detection <- CovOccPB #CovOccPB[,5] puede ser categorica, en este caso no se scala

y.so <- PA #read.csv("data/yb-so.csv")


print("removing NA")
is.complete=complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*3.20
  

print('removing raster files from memory')
#removing rasters to free the memory
remove(bio_vars) 


print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.intensity),values(s.intensity)))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))

# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

print("specifying area ")
#area in squared km ??????????????????????? -----------------------------------
area.back = rep((xres(s.intensity)/1000)*(yres(s.intensity)/1000), nrow(X.back))# each cell

s.area=area.back*nrow(X.back) #study area

# 1. Preparing the data ========================================================
#Preparing Presence-Only data ------------------------------------------------

#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.intensity)), values(s.intensity))
colnames(X.back)=c("",names(s.intensity))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))
colnames(W.back)=c("",names(s.detection))
# remove all NA values
tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km -----------------------------------
area.back = rep((xres(s.intensity)/1000)*(yres(s.intensity)/1000), nrow(X.back))# each cell
s.area=area.back*nrow(X.back) #study area

# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.intensity))), pb.intensity)
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)


#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)

W.so = array(dim=c(nrow(as.matrix(so.detection)), J.so, 4))
W.so[,,1] = 1
W.so[,,2] = matrix(rep(so.detection$s.detection.1,J.so),ncol = J.so)# if it changes
W.so[,,3] = matrix(rep(so.detection$s.detection.2,J.so),ncol = J.so)
W.so[,,4] = matrix(rep(so.detection$s.detection.3,J.so),ncol = J.so)
W.so[,,5] = matrix(rep(so.detection$s.detection.4,J.so),ncol = J.so)
# 2. Analising the data ========================================================

#Analyzing Presence-Only data
pb.fit=pb.ipp(X.po, W.po,X.back, W.back)

# Analyzing presence-absence data
so.fit=so.model(X.so,W.so,y.so)

# Analyzing presence-only data AND presence-absence data
poANDso.fit=pbso.integrated(X.po, W.po,X.back, W.back,X.so,W.so,y.so)
