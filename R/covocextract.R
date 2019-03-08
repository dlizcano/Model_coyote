library(raster)
library(fields)
library(mvtnorm)
library(matrixStats)
library(rgbif) # get GBIF data
library(maptools)

cam <- read.table("data/camarasWGS84.txt", sep=",", header = TRUE)
coor <- cam[,c(3,4)]


campoints <- SpatialPoints(cbind(coor$lon, 
                                 coor$lat))

proj4string(campoints) <- CRS('+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0')
centroid <- c(mean(coordinates(campoints)[,1]),mean(coordinates(campoints)[,2])) 

Mexico <- getData(name="GADM", country="MEX", level=1)
plot(Mexico)
plot(campoints, add=T)

bio_varsc <- getData('worldclim', var='bio' , res=0.5, lon=centroid[1], lat=centroid[2]) 

plot(bio_varsc[[3]])
plot(campoints, add=T)

s.intensity <- stack("s.intensity.tif")
s.detection <- stack("s.detection.tif")

cov.oc=extract(s.intensity, campoints, method='bilinear') # deja 4 puntos fuera de la ventana por fuera
cov.det=extract(s.detection, campoints, method='bilinear')



write.csv(cov.oc, "cov.oc.csv")
write.csv(cov.det, "cov.det.csv")

