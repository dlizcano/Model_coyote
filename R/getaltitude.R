library(raster)

altitude1 <- getData('alt', country="USA")
USA <- mosaic(altitude1[[1]],altitude1[[2]], fun=min)
Canada <- getData('alt', country="CAN")
Mexico <- getData('alt', country="MEX")
Guatemala <- getData('alt', country="GTM")
Belice <- getData('alt', country="BLZ")
Salvador <- getData('alt', country="SLV")
Honduras <- getData('alt', country="HND")
Nicaragua <- getData('alt', country="NIC")
CostaR <- getData('alt', country="CRI")
Panama <- getData('alt', country="PAN")

altitude <- mosaic(USA,Canada,Mexico,Guatemala,Belice,Salvador,Honduras,Nicaragua,CostaR, Panama, fun=min)

writeRaster(altitude,"altitude.asc")
