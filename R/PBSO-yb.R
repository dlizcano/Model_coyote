###################################################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection
## fitting PO PA and Integrated models the data
## 19/08/2016
###################################################################################################
######### Data required for the function ##########################################################
###################################################################################################
## All the data required for this function is stored in a file data.rda
## s.occupancy - raster with background covatiates that effect occupancy
## s.detection - raster with background covariates that effect detection
## pb.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
## pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
##
## y.so - matrix of detection/non detection of the PA surveys
## so.occupancy - matrix with covariates that effect occupancy in the locations of PA survey sites
## so.detection - matrix with covariates that effect detection in the locations of PA survey sites in each survey
###################################################################################################

require(raster)
require(fields)
require(mvtnorm)
require(matrixStats)

source("R/functions.r")



# Data Preparation -----------------------------------------------------------------------------
#### read spatial covariates
x.files=c( "raster/Dem75mInteger_100_no"
					 # ,"MinTemp_Jul_100_no"
					 # ,"MaxTemp_Jan_100_no"
					 # ,"log_vertical_distance_major_streams_sept2012_100_no"
					 # ,"wetness_index_saga_sept2012_100_no"
					 # ,'Evaporation_Jul_100_no.grd'
					 ,'raster/Evaporation_Jan_100_no.grd'
					 # ,'RainDays_Jul_100_no.gri'
					 # ,'Rainfall_Jul_100_no.grd'
					 # ,'RainDays_Jan_100_no.gri'
					 # ,'Rainfall_Jan_100_no.grd'
					 # ,'visible_sky_sept2012_100_no'
)

w.files=c(#"distance_to_roads_100_no",
          "raster/terrain_ruggedness_index_sept2012_100_no" )

for (i in c(x.files,w.files)) {
	do.call('=',list(i, raster(i)))
}
s.occupancy= raster(get(x.files[1]))

for (i in x.files) {
	temp=get(i)
	names(temp) = i
	s.occupancy = addLayer(s.occupancy, temp)
}

s.detection= raster(get(w.files[1]))
for (i in w.files) {
	temp=get(i)
	names(temp) = i
	s.detection = addLayer(s.detection, temp)
}

# #plotting occupancy and detection covariates distributions in the study area
# print('plotting occupancy and detection rasters')
# ppi = 300
# png('occupancy-covariates.png', width=9*ppi, height=round(length(x.files)/2+2)*ppi, res=ppi)
# plot(s.occupancy)
# dev.off()
# 
# png('pb-detection-covariates.png', width=9*ppi, height=round(length(w.files)/2+2)*ppi, res=ppi)
# plot(s.detection)
# dev.off()

# read pb and PA survey data
pb=read.csv("data/yb-pb-location.csv")


pb.loc=SpatialPoints(pb)
pb.occupancy=extract(s.occupancy,pb.loc)
pb.detection=extract(s.detection,pb.loc)

is.complete.pb=complete.cases(pb.detection)&complete.cases(pb.occupancy)
pb.detection=pb.detection[is.complete.pb,]
pb.occupancy=pb.occupancy[is.complete.pb,]

so.occupancy=read.csv("data/so_occupancy.csv")
so.detection=read.csv("data/so_detection.csv")
y.so=read.csv("data/yb-so.csv")

print("removing NA")
is.complete=complete.cases(so.occupancy)&complete.cases(so.detection)&complete.cases(y.so)
so.occupancy=so.occupancy[is.complete,]# only use complete cases
so.detection=so.detection[is.complete,]# only use complete cases
y.so=y.so[is.complete,]# only use complete cases
area.so =pi*0.04

print('removing raster files')
#removing rasters to free the memory
for (i in c(x.files,w.files)) {
	do.call(remove,list(i))
}
remove(is.complete.pb,po,yb.so,po.loc)


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

# adding column of ones - po locations
# Error change po by pb
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy)
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)



#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
### Error where is PA? <-  changeg to pb
X.so=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy)
W.so = array(dim=c(nrow(as.matrix(pb.detection)), J.so, 3))
W.so[,,1] = 1
W.so[,,2] = pb.detection# if it changes 
W.so[,,3] = pb.detection# if it changes





#Checking whether occupancy and detection rasters have the same resolution -----
if(sum(res(s.occupancy)!=res(s.detection)))
	stop("Occupancy and detection raster layers have different resolution")

if(ncell(s.occupancy)!=ncell(s.detection))
	stop("Occupancy and detection have different number of cells")


# # Plotting covariates that drive occupancy and detection in PO
# ppi = 300
# png('occupancy-covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.occupancy)
# dev.off()

# png('PO-detection -covariates.png', width=9*ppi, height=3*ppi, res=ppi)
plot(s.detection)
# dev.off()

# 1. Preparing the data ========================================================
#Preparing Presence-Only data ------------------------------------------------

#turning rasters into tables - background, adding column of ones
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
X.pb=cbind(rep(1, nrow(as.matrix(pb.occupancy))), pb.occupancy)
W.pb=cbind(rep(1, nrow(as.matrix(pb.detection))), pb.detection)



#Preparing Presence-Absence data ------------------------------------------------

#add a column of ones to the PA covariat
#y.so # matrix of presences and absences (when the species was and wasn't present)
J.so=ncol(y.so)
X.so=cbind(rep(1, nrow(as.matrix(so.occupancy))), so.occupancy)
W.so = array(dim=c(nrow(as.matrix(so.detection)), J.so, 3))
W.so[,,1] = 1
W.so[,,2] = so.detection[,1:2]# if it changes
# W.so[,,3] = so.detection[,3:4]# if it changes
# 2. Analising the data ========================================================

#Analyzing Presence-Only data
pb.fit=pb.ipp(X.pb, W.pb,X.back, W.back)

# Analyzing presence-absence data
so.fit=so.model(X.so,W.so,y.so)

# Analyzing presence-only data AND presence-absence data
poANDso.fit=pbso.integrated(X.pb, W.pb,X.back, W.back,X.so,W.so,y.so)


