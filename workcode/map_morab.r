############################################################################

options(java.parameters = "-Xmx4g" )

library(XLConnect)
library(raster)
library(sp)
library(rgeos)
library(leaflet)
library(htmlwidgets)
library(RColorBrewer)

############################################################################
# Define a working dir

setwd("")


############################################################################
# Define a function for sp clipping

crop_raster_pols=function(sppolydf,r) {
                                      require(raster)
                                      require(sp)
 
                                     if (class(r)!="RasterLayer") {stop("Is not a R raster object")}
                                     if (proj4string(r)!=proj4string(sppolydf)) {stop("Check projection of spatial objects")}
                                     r.sub <- crop(r, extent(sppolydf))
                                     r.sub <- mask(r.sub, sppolydf)
                                     return(r.sub)
}

############################################################################


#############################################################################
firenze_data = read.csv("firenze_cs_LST.csv")
names(firenze_data)[2] = "built-up"

firenze_LST_poly = readRDS("firenze_LST_poly.rds")


firenze_AATmap_LST = readRDS("firenze_AATmap_LST.rds")

firenze_LST_center=gCentroid(firenze_LST_poly,byid=TRUE)

firenze_index_LST=as.data.frame(extract(stack(firenze_AATmap_LST),firenze_LST_center))

names(firenze_index_LST)<-names(firenze_AATmap_LST)

writeWorksheetToFile("firenze_slope_LST.xls",firenze_index_LST,sheet="Parametri Trend")

firenze_LST_center=SpatialPointsDataFrame(firenze_LST_center,cbind(firenze_data,firenze_index_LST))

pixels_LST_firenze <- SpatialPixelsDataFrame(firenze_LST_center, tolerance = 0.01, firenze_LST_center@data)


firenze_LST_center_df=as.data.frame(firenze_LST_center)

firenze_LST_poly_bound=gUnionCascaded(firenze_LST_poly, id = NULL)
proj4string(firenze_LST_poly_bound)=proj4string(firenze_LST_poly)
firenze_LST_poly_bound_p<- spTransform(firenze_LST_poly_bound, CRS("+init=epsg:3857"))


raster_LST_firenze_p <- projectRaster(stack(pixels_LST_firenze), crs="+init=epsg:3857")

raster_LST_firenze=raster_LST_firenze_p$Annual
raster_pvalue_firenze=raster_LST_firenze_p$PvalSEG1
raster_trend_firenze=raster_LST_firenze_p$SlopeSEG1

raster_LST_firenze_n=crop_raster_pols(firenze_LST_poly_bound_p,raster_LST_firenze)
raster_pvalue_firenze_n=crop_raster_pols(firenze_LST_poly_bound_p,raster_pvalue_firenze)
raster_trend_firenze_n=crop_raster_pols(firenze_LST_poly_bound_p,raster_trend_firenze)


pvalues_color=gray(seq(0,1,0.2))
pvalues_color[1]="#ff0000"
pal_LST <- colorNumeric(rev(heat.colors(14)),seq(14,21),na.color = "transparent")
pal_pvalue <- colorNumeric(pvalues_color,c(seq(0.05,1,0.1),1),na.color = "transparent")
pal_trend <- colorNumeric(rev(colorRampPalette(brewer.pal(9,"RdBu"))(100)),rev(seq(-0.02,0.02,0.01)),na.color = "transparent")

lat_mean_firenze=mean(firenze_LST_center_df$y)
lon_mean_firenze=mean(firenze_LST_center_df$x)

m_LST_firenze=leaflet() %>% setView(lon_mean_firenze,lat_mean_firenze,zoom=12) %>% addTiles() %>% addRasterImage(raster_LST_firenze_n, colors = pal_LST, opacity = 0.7,project=F) %>%
                           addCircleMarkers(data=firenze_LST_center_df, lng = ~x, lat = ~y,radius = ~0.5*sqrt(built.up), stroke=FALSE, fillOpacity=0.5,color = c('black')) %>%
                           addLegend(pal = pal_LST, values = seq(14,21),title = "Annual Mean LST (°C)")

saveWidget(m_LST_firenze, file="firenze_map_csuolo.html")


m_pvalue_firenze=leaflet() %>% setView(lon_mean_firenze,lat_mean_firenze,zoom=12) %>% addTiles() %>% addRasterImage(raster_pvalue_firenze_n, colors = pal_pvalue, opacity = 0.7,project=F) %>%
  addLegend(pal = pal_pvalue, values = c(seq(0.05,1,0.1),1),labFormat =labelFormat(prefix=" ",digits=5,transform= function(x) x-0.05),title = "p-value LST linear trend")

saveWidget(m_pvalue_firenze, file="firenze_pvalue_LST.html")

m_trend_firenze=leaflet() %>% setView(lon_mean_firenze,lat_mean_firenze,zoom=12) %>% addTiles() %>% addRasterImage(raster_trend_firenze_n, colors = pal_trend, opacity = 0.7,project=F) %>%
  addLegend(pal = pal_trend, values = rev(seq(-0.02,0.02,0.01)),labFormat =labelFormat(prefix="  ",digits=5),title = "Slope LST linear trend")

saveWidget(m_trend_firenze, file="firenze_trend_LST.html")



############################################################################

bologna_data = read.csv("bologna_cs_LST.csv")
names(bologna_data)[2] = "built-up"

bologna_LST_poly = readRDS("bologna_LST_poly.rds")


bologna_AATmap_LST = readRDS("bologna_AATmap_LST.rds")

bologna_LST_center=gCentroid(bologna_LST_poly,byid=TRUE)

bologna_index_LST=as.data.frame(extract(stack(bologna_AATmap_LST),bologna_LST_center))

names(bologna_index_LST)<-names(bologna_AATmap_LST)

writeWorksheetToFile("bologna_slope_LST.xls",bologna_index_LST,sheet="Parametri Trend")

bologna_LST_center=SpatialPointsDataFrame(bologna_LST_center,cbind(bologna_data,bologna_index_LST))

pixels_LST_bologna <- SpatialPixelsDataFrame(bologna_LST_center, tolerance = 0.01, bologna_LST_center@data)


bologna_LST_center_df=as.data.frame(bologna_LST_center)
bologna_LST_poly_bound=gUnionCascaded(bologna_LST_poly, id = NULL)
proj4string(bologna_LST_poly_bound)=proj4string(bologna_LST_poly)
bologna_LST_poly_bound_p<- spTransform(bologna_LST_poly_bound, CRS("+init=epsg:3857"))


raster_LST_bologna_p <- projectRaster(stack(pixels_LST_bologna), crs="+init=epsg:3857")

raster_LST_bologna=raster_LST_bologna_p$Annual
raster_pvalue_bologna=raster_LST_bologna_p$PvalSEG1
raster_trend_bologna=raster_LST_bologna_p$SlopeSEG1

raster_LST_bologna_n=crop_raster_pols(bologna_LST_poly_bound_p,raster_LST_bologna)
raster_pvalue_bologna_n=crop_raster_pols(bologna_LST_poly_bound_p,raster_pvalue_bologna)
raster_trend_bologna_n=crop_raster_pols(bologna_LST_poly_bound_p,raster_trend_bologna)

lat_mean_bologna=mean(bologna_LST_center_df$y)
lon_mean_bologna=mean(bologna_LST_center_df$x)

m_LST_bologna=leaflet() %>% setView(lon_mean_bologna,lat_mean_bologna,zoom=12)  %>% addTiles() %>% addRasterImage(raster_LST_bologna_n, colors = pal_LST, opacity = 0.7,project=F) %>%
  addCircleMarkers(data=bologna_LST_center_df, lng = ~x, lat = ~y,radius = ~0.5*sqrt(built.up), stroke=FALSE, fillOpacity=0.5,color = c('black')) %>%
  addLegend(pal = pal_LST, values = seq(14,21),title = "Annual Mean LST (°C)")
saveWidget(m_LST_bologna, file="bologna_map_csuolo.html")

m_pvalue_bologna=leaflet() %>% setView(lon_mean_bologna,lat_mean_bologna,zoom=12) %>% addTiles() %>% addRasterImage(raster_pvalue_bologna_n, colors = pal_pvalue, opacity = 0.7,project=F) %>%
  addLegend(pal = pal_pvalue, values = c(seq(0.05,1,0.1),1),labFormat =labelFormat(prefix=" ",digits=5,transform= function(x) x-0.05),title = "p-value LST linear trend")
saveWidget(m_pvalue_bologna, file="bologna_pvalue_LST.html")

m_trend_bologna=leaflet() %>% setView(lon_mean_bologna,lat_mean_bologna,zoom=12) %>% addTiles() %>% addRasterImage(raster_trend_bologna_n, colors = pal_trend, opacity = 0.7,project=F) %>%
 addLegend(pal = pal_trend, values = rev(seq(-0.02,0.02,0.01)),labFormat =labelFormat(prefix=" ",digits=5),title = "Slope LST linear trend")
saveWidget(m_trend_bologna, file="bologna_trend_LST.html")


############################################################################


milano_data = read.csv("milano_cs_LST.csv")
names(milano_data)[2] = "built-up"

milano_LST_poly = readRDS("milano_LST_poly.rds")


milano_AATmap_LST = readRDS("milano_AATmap_LST.rds")

milano_LST_center=gCentroid(milano_LST_poly,byid=TRUE)

milano_index_LST=as.data.frame(extract(stack(milano_AATmap_LST),milano_LST_center))

names(milano_index_LST)<-names(milano_AATmap_LST)

writeWorksheetToFile("milano_slope_LST.xls",milano_index_LST,sheet="Parametri Trend")

milano_LST_center=SpatialPointsDataFrame(milano_LST_center,cbind(milano_data,milano_index_LST))

pixels_LST_milano <- SpatialPixelsDataFrame(milano_LST_center, tolerance = 0.01, milano_LST_center@data)


milano_LST_center_df=as.data.frame(milano_LST_center)
milano_LST_poly_bound=gUnionCascaded(milano_LST_poly, id = NULL)
proj4string(milano_LST_poly_bound)=proj4string(milano_LST_poly)
milano_LST_poly_bound_p<- spTransform(milano_LST_poly_bound, CRS("+init=epsg:3857"))


raster_LST_milano_p <- projectRaster(stack(pixels_LST_milano), crs="+init=epsg:3857")

raster_LST_milano=raster_LST_milano_p$Annual
raster_pvalue_milano=raster_LST_milano_p$PvalSEG1
raster_trend_milano=raster_LST_milano_p$SlopeSEG1

raster_LST_milano_n=crop_raster_pols(milano_LST_poly_bound_p,raster_LST_milano)
raster_pvalue_milano_n=crop_raster_pols(milano_LST_poly_bound_p,raster_pvalue_milano)
raster_trend_milano_n=crop_raster_pols(milano_LST_poly_bound_p,raster_trend_milano)

lat_mean_milano=mean(milano_LST_center_df$y)
lon_mean_milano=mean(milano_LST_center_df$x)

m_LST_milano=leaflet() %>% setView(lon_mean_milano,lat_mean_milano,zoom=12)  %>% addTiles() %>% addRasterImage(raster_LST_milano_n, colors = pal_LST, opacity = 0.7,project=F) %>%
  addCircleMarkers(data=milano_LST_center_df, lng = ~x, lat = ~y,radius = ~0.5*sqrt(built.up), stroke=FALSE, fillOpacity=0.5,color = c('black')) %>%
  addLegend(pal = pal_LST, values = seq(14,21),title = "Annual Mean LST (°C)")
 saveWidget(m_LST_milano, file="milano_map_csuolo.html")

m_pvalue_milano=leaflet() %>% setView(lon_mean_milano,lat_mean_milano,zoom=12) %>% addTiles() %>% addRasterImage(raster_pvalue_milano_n, colors = pal_pvalue, opacity = 0.7,project=F) %>%
   addLegend(pal = pal_pvalue, values = c(seq(0.05,1,0.1),1),labFormat =labelFormat(prefix=" ",digits=5,transform= function(x) x-0.05),title = "p-value LST linear trend")

saveWidget(m_pvalue_milano, file="milano_pvalue_LST.html")

m_trend_milano=leaflet()%>% setView(lon_mean_milano,lat_mean_milano,zoom=12) %>% addTiles() %>% addRasterImage(raster_trend_milano_n, colors = pal_trend, opacity = 0.7,project=F) %>%
   addLegend(pal = pal_trend, values = rev(seq(-0.02,0.02,0.01)),labFormat =labelFormat(prefix="  ",digits=5),title = "Slope LST linear trend")

saveWidget(m_trend_milano, file="milano_trend_LST.html")


############################################################################

roma_data=read.csv("roma_cs_LST.csv")

names(roma_data)[2] = "built-up"

roma_LST_poly = readRDS("roma_LST_poly.rds")


roma_AATmap_LST = readRDS("roma_AATmap_LST.rds")

roma_LST_center=gCentroid(roma_LST_poly,byid=TRUE)

roma_index_LST=as.data.frame(extract(stack(roma_AATmap_LST),roma_LST_center))

names(roma_index_LST)<-names(roma_AATmap_LST)

writeWorksheetToFile("roma_slope_LST.xls",roma_index_LST,sheet="Parametri Trend")

roma_LST_center=SpatialPointsDataFrame(roma_LST_center,cbind(roma_data,roma_index_LST))

pixels_LST_roma <- SpatialPixelsDataFrame(roma_LST_center, tolerance = 0.01, roma_LST_center@data)


roma_LST_center_df=as.data.frame(roma_LST_center)
roma_LST_poly_bound=gUnionCascaded(roma_LST_poly, id = NULL)
proj4string(roma_LST_poly_bound)=proj4string(roma_LST_poly)
roma_LST_poly_bound_p<- spTransform(roma_LST_poly_bound, CRS("+init=epsg:3857"))


raster_LST_roma_p <- projectRaster(stack(pixels_LST_roma), crs="+init=epsg:3857")

raster_LST_roma=raster_LST_roma_p$Annual
raster_pvalue_roma=raster_LST_roma_p$PvalSEG1
raster_trend_roma=raster_LST_roma_p$SlopeSEG1

raster_LST_roma_n=crop_raster_pols(roma_LST_poly_bound_p,raster_LST_roma)
raster_pvalue_roma_n=crop_raster_pols(roma_LST_poly_bound_p,raster_pvalue_roma)
raster_trend_roma_n=crop_raster_pols(roma_LST_poly_bound_p,raster_trend_roma)

lat_mean_roma=mean(roma_LST_center_df$y)
lon_mean_roma=mean(roma_LST_center_df$x)

m_LST_roma=leaflet() %>% setView(lon_mean_roma,lat_mean_roma,zoom=12)  %>% addTiles() %>% addRasterImage(raster_LST_roma_n, colors = pal_LST, opacity = 0.7,project=F) %>%
  addCircleMarkers(data=roma_LST_center_df, lng = ~x, lat = ~y,radius = ~0.8*sqrt(built.up), stroke=FALSE, fillOpacity=0.5,color = c('black')) %>%
  addLegend(pal = pal_LST, values = seq(14,21),title = "Annual Mean LST (°C)")
 setView(m_LST_roma,40.0,11.0,zoom =15)
saveWidget(m_LST_roma, file="roma_map_csuolo.html")


m_pvalue_roma=leaflet() %>% setView(lon_mean_roma,lat_mean_roma,zoom=12) %>% addTiles() %>% addRasterImage(raster_pvalue_roma_n, colors = pal_pvalue, opacity = 0.7,project=F) %>%
  addLegend(pal = pal_pvalue, values =c(seq(0.05,1,0.1),1),labFormat =labelFormat(prefix=" ",digits=5,transform= function(x) x-0.05),title = "p-value LST linear trend")

saveWidget(m_pvalue_roma, file="roma_pvalue_LST.html")

m_trend_roma=leaflet() %>% setView(lon_mean_roma,lat_mean_roma,zoom=12)  %>% addTiles() %>% addRasterImage(raster_trend_roma_n, colors = pal_trend, opacity = 0.7,project=F) %>%
  addLegend(pal = pal_trend, values = rev(seq(-0.02,0.02,0.01)),labFormat =labelFormat(prefix=" < ",digits=5),title = "Slope LST linear trend")

saveWidget(m_trend_roma, file="roma_trend_LST.html")


############################################################################
# Build a legend in leaf let of built up

lat=seq(43.7503, 43.82650,0.00762)[1:9]
lon=rep(as.data.frame(firenze_LST_center_df)$x[54],9)
legend_cs=data.frame(lat=lat,lon=lon)
coordinates(legend_cs)=~ lon+lat
legend_cs$cs=seq(10,90,10)
leg=leaflet() %>% addCircleMarkers(data=legend_cs,radius = ~sqrt(cs), stroke=FALSE, fillOpacity=0.5,color = c('black'))
saveWidget(leg, file="legend_csuolo.html")



#
