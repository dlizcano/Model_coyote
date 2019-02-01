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

source("R/functions.r") # read functions

####################################
### Get GBIF Data per species
####################################
sp_in_GBIF <- occ_search(scientificName = "Canis latrans", return='data')
sp_in_GBIF <- read.csv("Cala_gbif.csv") # read prescence points file

# Plot to check GBIF data... maybe some point tweeking required !
# gbifmap(sp_in_GBIF, region = "Colombia") # plot in a world map

# make spatial point
sp_points <- SpatialPoints(cbind(sp_in_GBIF$longitude, 
                           sp_in_GBIF$latitude))

# put georeference
proj4string(sp_points) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
centroid <- c(mean(coordinates(sp_points)[,1]),mean(coordinates(sp_points)[,2])) 
####################################
### Get altitude and worldclim Data for the species
####################################
altitude1 <- getData('SRTM', lon=centroid[1], lat=centroid[2])
altitude2 <- getData('SRTM', lon=centroid[1], lat=centroid[2]-2) # to get the lower part
alt_mosaic <- mosaic(altitude1, altitude2, fun=mean) # unir las dos altitudes
# tmin <-     getData('worldclim', var='tmin', res=0.5, lon=centroid[1], lat=centroid[2])
# tmax <-     getData('worldclim', var='tmax', res=0.5, lon=centroid[1], lat=centroid[2])
# precipit <- getData('worldclim', var='prec', res=0.5, lon=centroid[1], lat=centroid[2])
bio_vars <- getData('worldclim', var='bio' , res=0.5, lon=centroid[1], lat=centroid[2])

########################################################################
## cut by area of interest   
## in this case mask and extent by Prec (mapa probided by Angelica)              
########################################################################

# read Angelica´s raster
Prec <- raster('Prec.asc') # equivale al M
# Prec <- readShapePoly('shp/Unidad_analisis.shp')
proj4string(Prec) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
### cut Bioclim
bio_vars_crop <- crop(bio_vars, extent(Prec)) #crop by a raster. Precipitation
extent(bio_vars_crop) <- alignExtent(bio_vars_crop, Prec) #Aligne extent
bio_vars_cut <- mask(bio_vars_crop, Prec)# mask

### cut altutude
alti_crop <- crop(alt_mosaic, extent(Prec))
alti_crop_r <- resample(alti_crop, bio_vars_cut, method='bilinear')
extent(alti_crop_r) <- alignExtent(alti_crop_r, Prec)
alti_cut <- mask(alti_crop_r, Prec)
names(alti_cut) <- "alti" # put name


plot(bio_vars_cut[[2]])
plot(sp_points, add=T)

# Data Preparation -----------------------------------------------------------------------------
#### read spatial covariates from BIOCLIM and standarize using scale
s.occupancy <- stack(scale(alti_cut, center=TRUE, scale=TRUE),
                     scale(bio_vars_cut[[12]], center=TRUE, scale=TRUE), # Annual Precipitation
                     scale(bio_vars_cut[[1]], center=TRUE, scale=TRUE),  # Annual Mean Temperature
                     scale(bio_vars_cut[[3]], center=TRUE, scale=TRUE))#, # Isothermality (BIO2/BIO7) (* 100)
                     # scale(bio_vars_cut[[4]], center=TRUE, scale=TRUE), 
                     # scale(bio_vars_cut[[5]], center=TRUE, scale=TRUE),
                     # scale(bio_vars_cut[[6]], center=TRUE, scale=TRUE), 
                     # scale(bio_vars_cut[[7]], center=TRUE, scale=TRUE),
                     # scale(bio_vars_cut[[12]], center=TRUE, scale=TRUE)) # Annual Precipitation

s.detection <-  stack(scale(alti_cut, center=TRUE, scale=TRUE),
                      scale(bio_vars_cut[[12]], center=TRUE, scale=TRUE), # Annual Precipitation
                      scale(bio_vars_cut[[1]], center=TRUE, scale=TRUE),  # Annual Precipitation
                      scale(bio_vars_cut[[3]], center=TRUE, scale=TRUE))#, # Isothermality (BIO2/BIO7) (* 100)
                      # scale(bio_vars_cut[[13]], center=TRUE, scale=TRUE), # BIO13 = Precipitation of Wettest Month
                      # scale(bio_vars_cut[[14]], center=TRUE, scale=TRUE), # BIO14 = Precipitation of Driest Month
                      # scale(bio_vars_cut[[15]], center=TRUE, scale=TRUE)) # BIO15 = Precipitation Seasonality (Coefficient of Variation)
  
  
# plot... just to check 
plot(s.occupancy[[1]]) 
plot(sp_points, add=T)

# read pb and PA survey data. This is the only presence oportunistic data!
pb <- sp_points # coming from GBIF... no depurado


pb.loc <- pb
#### Get the covariates from rasters (stacked alti and Bioclim)
pb.occupancy=extract(s.occupancy, pb.loc, method='bilinear') # deja 4 puntos fuera de la ventana por fuera
pb.detection=extract(s.detection, pb.loc, method='bilinear') # deja 4 puntos fuera de la ventana por fuera

### deja solo puntos dentro de la ventana
is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,]
pb.occupancy=pb.occupancy[is.complete.pb,]

################ ################ ################ ################ 
################ READ DATA, Occupancy data measured in the field
################ ################ ################ ################ 

#### estos 3 archivos deben tener el mismo numero de columnas
# PA <- as.matrix(read.table('PA.txt', header=T, sep='\t')) # En PC_lenovo_ADP: PA <- as.matrix(read.table('PA.txt'))
# CovOcc <- read.table('CovOcc.txt', header=T, sep='\t') # En PC_lenovo_ADP: CovOcc <- read.table('CovOcc.txt')
# CovOccPB <- read.table('CovOccPB.txt', header=T, sep='\t') # En PC_lenovo_ADP: CovOccPB <- read.table('CovOccPB.txt')

humedad <- read.csv("data/Localidad_visita_humedad.txt", sep="" )
temperatura <- read.csv("data/localidad_visita_temp.txt", sep="" )
otras_cov <- read.csv("data/matriz_covariables_F.csv", sep=",")


PA <- as.matrix(read.csv("data/Ituango_final_conteos.txt", sep="")[,-1]) # En PC_lenovo_ADP: PA <- as.matrix(read.table('PA.txt'))
CovOcc <-  scale(otras_cov[-1])# En PC_lenovo_ADP: CovOcc <- read.table('CovOcc.txt')
CovOccPB <- scale(temperatura[-1]) # En PC_lenovo_ADP: CovOccPB <- read.table('CovOccPB.txt')


#### estos 3 archivos deben tener el mismo numero de columnas
# so.occupancy <- scale(CovOcc) #read.csv("data/so_occupancy.csv")
so.occupancy <- CovOcc

#### estos 2 archivos deben tener el mismo numero de filas
# so.detection <- cbind(scale(CovOccPB[,1:4]),CovOccPB[,5]) #CovOccPB[,5] puede ser categorica, en este caso no se scala
so.detection <- CovOccPB #CovOccPB[,5] puede ser categorica, en este caso no se scala

y.so <- PA #read.csv("data/yb-so.csv")


print("removing NA")
is.complete=complete.cases(so.occupancy)&complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*0.04

print('removing raster files from memory')
#removing rasters to free the memory
remove(bio_vars) 


print("allocating background")
#turning rasters into tables - background, adding column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))

# remove all NA values

tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

print("specifying area ")
#area in squared km ??????????????????????? -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell


s.area=area.back*nrow(X.back) #study area


# 1. Preparing the data ========================================================
#Preparing Presence-Only data ------------------------------------------------

#turning rasters into tables - background, adding a column of ones
X.back = cbind(rep(1, ncell(s.occupancy)), values(s.occupancy))
colnames(X.back)=c("",names(s.occupancy))
W.back = cbind(rep(1, ncell(s.detection)), values(s.detection))
colnames(W.back)=c("",names(s.detection))
# remove all NA values
tmp=X.back[complete.cases(X.back)&complete.cases(W.back),]
W.back=W.back[complete.cases(X.back)&complete.cases(W.back),]
X.back=tmp

#area in squared km -----------------------------------
area.back = rep((xres(s.occupancy)/1000)*(yres(s.occupancy)/1000), nrow(X.back))# each cell
s.area=area.back*nrow(X.back) #study area

# adding column of ones - po locations
X.po=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy)
W.po=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)



#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
W.so = array(dim=c(nrow(as.matrix(so.detection)), J.so, 2))
W.so[,,1] = 1
W.so[,,2] = so.detection# if it changes

# 2. Analising the data ========================================================

#Analyzing Presence-Only data
pb.fit=pb.ipp(X.po, W.po, X.back, W.back)

# Analyzing presence-absence data
so.fit=so.model(X.so, W.so, y.so)

# Analyzing presence-only data AND presence-absence data
poANDso.fit=pbso.integrated(X.po, W.po, X.back, W.back, X.so, W.so, y.so)



