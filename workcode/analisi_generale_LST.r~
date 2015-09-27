############################################################

library(raster)
library(rts)
library(OpenStreetMap)
library(rasterVis)
library(ggplot2)
library(rgdal)
library(ggmap)
library(grid)
library(latticeExtra)
library(maptools)
library(plyr)
library(greenbrown)
library(plotKML)


############################################################
# Setup directory

setwd("")

############################################################
# Load supplementary code

source("auxillary_functions.r")
data(SAGA_pal)
plotKML.env(convert=Sys.which("convert"),gdal_translate=Sys.which("gdal_translate"),gdalwarp=Sys.which("gdalwarp"),python=Sys.which("python"))
############################################################
# define proj4 modis

proj_modis="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

###############################################################################################
# milano break and trend analisys NDVI MODIS 

milano_istat=readRDS("data/milano_istat.rds")

milano_AATmap_ndvi=readRDS("data/milano_AATmap_ndvi.rds")

milano_slope_ndvi_b=milano_AATmap_ndvi[["SlopeSEG1"]]  
milano_break_ndvi_s_b<- TrendClassification(milano_AATmap_ndvi[["PvalSEG1"]], min.length=10, max.pval=0.05)

milano_slope_ndvi=rectify_matlab_array_mat(milano_slope_ndvi_b,inCRS=proj_modis)
milano_break_ndvi_s=rectify_matlab_array_mat(milano_break_ndvi_s_b,inCRS=proj_modis)

milano_slope_ndvi=crop_raster(milano_istat,milano_slope_ndvi)
milano_break_ndvi=crop_raster(milano_istat,milano_break_ndvi_s)


plotKML::kml(milano_slope_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 


plotKML::kml(milano_break_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 



##########################################################################################
# milano break and trend analisys LST MODIS Day

milano_AATmap_LST=readRDS("data/milano_AATmap_LST.rds")

milano_slope_LST=milano_AATmap_LST[["SlopeSEG1"]]  
milano_trend_LST_s<- TrendClassification(milano_AATmap_LST[["PvalSEG1"]], min.length=10, max.pval=0.1)

milano_slope_LST=crop_raster(milano_istat,milano_slope_LST)
milano_trend_LST_s<-crop_raster(milano_istat,milano_trend_LST_s)

plotKML::kml(milano_slope_LST, colour_scale = SAGA_pal[[1]], colour = SlopeSEG1,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 
plotKML::kml(milano_trend_LST_s, colour_scale = SAGA_pal[[1]], colour = TrendSEG1,width=25,height=200, z.lim = c(0,1),kmz=T) 


#######################################################################################################
# milano LST climate 

milano_LST_monthly_ts=readRDS("data/milano_LST_monthly.rds")

milano_LST_DJF<-subset(milano_LST_monthly_ts,grep("-01-|-02-|-12-",index(milano_LST_monthly_ts)))
milano_LST_DJF_mean=calc(milano_LST_DJF@raster,function(x) mean(x,na.rm=T))
milano_LST_DJF_sd=calc(milano_LST_DJF@raster,function(x) sd(x,na.rm=T))
saveRDS(milano_LST_DJF_mean,"data/milano_LST_DJF_mean.rds")
saveRDS(milano_LST_DJF_sd,"data/milano_LST_DJF_sd.rds")
plotKML::kml(milano_LST_DJF_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(milano_LST_DJF_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

             
milano_LST_MAM<-subset(milano_LST_monthly_ts,grep("-04-|-05-|-06-",index(milano_LST_monthly_ts)))
milano_LST_MAM_mean=calc(milano_LST_MAM@raster,function(x) mean(x,na.rm=T))
milano_LST_MAM_sd=calc(milano_LST_MAM@raster,function(x) sd(x,na.rm=T))
saveRDS(milano_LST_MAM_mean,"data/milano_LST_MAM_mean.rds")
saveRDS(milano_LST_MAM_sd,"data/milano_LST_MAM_sd.rds")
plotKML::kml(milano_LST_MAM_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(milano_LST_MAM_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


milano_LST_JLA<-subset(milano_LST_monthly_ts,grep("-06-|-07-|-08-",index(milano_LST_monthly_ts)))
milano_LST_JLA_mean=calc(milano_LST_JLA@raster,function(x) mean(x,na.rm=T))
milano_LST_JLA_sd=calc(milano_LST_JLA@raster,function(x) sd(x,na.rm=T))
saveRDS(milano_LST_JLA_mean,"data/milano_LST_JLA_mean.rds")
saveRDS(milano_LST_JLA_sd,"data/milano_LST_JLA_sd.rds")
plotKML::kml(milano_LST_JLA_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(milano_LST_JLA_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


milano_LST_SON<-subset(milano_LST_monthly_ts,grep("-09-|-10-|-11-",index(milano_LST_monthly_ts)))
milano_LST_SON_mean=calc(milano_LST_SON@raster,function(x) mean(x,na.rm=T))
milano_LST_SON_sd=calc(milano_LST_SON@raster,function(x) sd(x,na.rm=T))
saveRDS(milano_LST_SON_mean,"data/milano_LST_SON_mean.rds")
saveRDS(milano_LST_SON_sd,"data/milano_LST_SON_sd.rds")
plotKML::kml(milano_LST_SON_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(milano_LST_SON_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


###############################################################################################
# palermo break and trend analisys NDVI MODIS 

palermo_istat=readRDS("data/palermo_istat.rds")

palermo_AATmap_ndvi=readRDS("data/palermo_AATmap_ndvi.rds")

palermo_slope_ndvi_b=palermo_AATmap_ndvi[["SlopeSEG1"]]  
palermo_break_ndvi_s_b<- TrendClassification(palermo_AATmap_ndvi[["PvalSEG1"]], min.length=10, max.pval=0.05)

palermo_slope_ndvi=rectify_matlab_array_mat(palermo_slope_ndvi_b,inCRS=proj_modis)
palermo_break_ndvi_s=rectify_matlab_array_mat(palermo_break_ndvi_s_b,inCRS=proj_modis)

palermo_slope_ndvi=crop_raster(palermo_istat,palermo_slope_ndvi)
palermo_break_ndvi=crop_raster(palermo_istat,palermo_break_ndvi_s)


plotKML::kml(palermo_slope_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 


plotKML::kml(palermo_break_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 



##########################################################################################
# palermo break and trend analisys LST MODIS Day

palermo_AATmap_LST=readRDS("data/palermo_AATmap_LST.rds")

palermo_slope_LST=palermo_AATmap_LST[["SlopeSEG1"]]  
palermo_trend_LST_s<- TrendClassification(palermo_AATmap_LST[["PvalSEG1"]], min.length=10, max.pval=0.1)

palermo_slope_LST=crop_raster(palermo_istat,palermo_slope_LST)
palermo_trend_LST_s<-crop_raster(palermo_istat,palermo_trend_LST_s)

plotKML::kml(palermo_slope_LST, colour_scale = SAGA_pal[[1]], colour = SlopeSEG1,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 
plotKML::kml(palermo_trend_LST_s, colour_scale = SAGA_pal[[1]], colour = TrendSEG1,width=25,height=200, z.lim = c(0,1),kmz=T) 


#######################################################################################################
# palermo LST climate 

palermo_LST_monthly_ts=readRDS("data/palermo_LST_monthly.rds")

palermo_LST_DJF<-subset(palermo_LST_monthly_ts,grep("-01-|-02-|-12-",index(palermo_LST_monthly_ts)))
palermo_LST_DJF_mean=calc(palermo_LST_DJF@raster,function(x) mean(x,na.rm=T))
palermo_LST_DJF_sd=calc(palermo_LST_DJF@raster,function(x) sd(x,na.rm=T))
saveRDS(palermo_LST_DJF_mean,"data/palermo_LST_DJF_mean.rds")
saveRDS(palermo_LST_DJF_sd,"data/palermo_LST_DJF_sd.rds")
plotKML::kml(palermo_LST_DJF_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(palermo_LST_DJF_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

             
palermo_LST_MAM<-subset(palermo_LST_monthly_ts,grep("-04-|-05-|-06-",index(palermo_LST_monthly_ts)))
palermo_LST_MAM_mean=calc(palermo_LST_MAM@raster,function(x) mean(x,na.rm=T))
palermo_LST_MAM_sd=calc(palermo_LST_MAM@raster,function(x) sd(x,na.rm=T))
saveRDS(palermo_LST_MAM_mean,"data/palermo_LST_MAM_mean.rds")
saveRDS(palermo_LST_MAM_sd,"data/palermo_LST_MAM_sd.rds")
plotKML::kml(palermo_LST_MAM_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(palermo_LST_MAM_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


palermo_LST_JLA<-subset(palermo_LST_monthly_ts,grep("-06-|-07-|-08-",index(palermo_LST_monthly_ts)))
palermo_LST_JLA_mean=calc(palermo_LST_JLA@raster,function(x) mean(x,na.rm=T))
palermo_LST_JLA_sd=calc(palermo_LST_JLA@raster,function(x) sd(x,na.rm=T))
saveRDS(palermo_LST_JLA_mean,"data/palermo_LST_JLA_mean.rds")
saveRDS(palermo_LST_JLA_sd,"data/palermo_LST_JLA_sd.rds")
plotKML::kml(palermo_LST_JLA_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(palermo_LST_JLA_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


palermo_LST_SON<-subset(palermo_LST_monthly_ts,grep("-09-|-10-|-11-",index(palermo_LST_monthly_ts)))
palermo_LST_SON_mean=calc(palermo_LST_SON@raster,function(x) mean(x,na.rm=T))
palermo_LST_SON_sd=calc(palermo_LST_SON@raster,function(x) sd(x,na.rm=T))
saveRDS(palermo_LST_SON_mean,"data/palermo_LST_SON_mean.rds")
saveRDS(palermo_LST_SON_sd,"data/palermo_LST_SON_sd.rds")
plotKML::kml(palermo_LST_SON_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(palermo_LST_SON_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


##################################################################################

###############################################################################################
# roma break and trend analisys NDVI MODIS 

roma_istat=readRDS("data/roma_istat.rds")

roma_AATmap_ndvi=readRDS("data/roma_AATmap_ndvi.rds")

roma_slope_ndvi_b=roma_AATmap_ndvi[["SlopeSEG1"]]  
roma_break_ndvi_s_b<- TrendClassification(roma_AATmap_ndvi[["PvalSEG1"]], min.length=10, max.pval=0.05)

roma_slope_ndvi=rectify_matlab_array_mat(roma_slope_ndvi_b,inCRS=proj_modis)
roma_break_ndvi_s=rectify_matlab_array_mat(roma_break_ndvi_s_b,inCRS=proj_modis)

roma_slope_ndvi=crop_raster(roma_istat,roma_slope_ndvi)
roma_break_ndvi=crop_raster(roma_istat,roma_break_ndvi_s)


plotKML::kml(roma_slope_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 


plotKML::kml(roma_break_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 



##########################################################################################
# roma break and trend analisys LST MODIS Day

roma_AATmap_LST=readRDS("data/roma_AATmap_LST.rds")

roma_slope_LST=roma_AATmap_LST[["SlopeSEG1"]]  
roma_trend_LST_s<- TrendClassification(roma_AATmap_LST[["PvalSEG1"]], min.length=10, max.pval=0.1)

roma_slope_LST=crop_raster(roma_istat,roma_slope_LST)
roma_trend_LST_s<-crop_raster(roma_istat,roma_trend_LST_s)

plotKML::kml(roma_slope_LST, colour_scale = SAGA_pal[[1]], colour = SlopeSEG1,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 
plotKML::kml(roma_trend_LST_s, colour_scale = SAGA_pal[[1]], colour = TrendSEG1,width=25,height=200, z.lim = c(0,1),kmz=T) 


#######################################################################################################
# roma LST climate 

roma_LST_monthly_ts=readRDS("data/roma_LST_monthly.rds")

roma_LST_DJF<-subset(roma_LST_monthly_ts,grep("-01-|-02-|-12-",index(roma_LST_monthly_ts)))
roma_LST_DJF_mean=calc(roma_LST_DJF@raster,function(x) mean(x,na.rm=T))
roma_LST_DJF_sd=calc(roma_LST_DJF@raster,function(x) sd(x,na.rm=T))
saveRDS(roma_LST_DJF_mean,"data/roma_LST_DJF_mean.rds")
saveRDS(roma_LST_DJF_sd,"data/roma_LST_DJF_sd.rds")
plotKML::kml(roma_LST_DJF_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(roma_LST_DJF_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

             
roma_LST_MAM<-subset(roma_LST_monthly_ts,grep("-04-|-05-|-06-",index(roma_LST_monthly_ts)))
roma_LST_MAM_mean=calc(roma_LST_MAM@raster,function(x) mean(x,na.rm=T))
roma_LST_MAM_sd=calc(roma_LST_MAM@raster,function(x) sd(x,na.rm=T))
saveRDS(roma_LST_MAM_mean,"data/roma_LST_MAM_mean.rds")
saveRDS(roma_LST_MAM_sd,"data/roma_LST_MAM_sd.rds")
plotKML::kml(roma_LST_MAM_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(roma_LST_MAM_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


roma_LST_JLA<-subset(roma_LST_monthly_ts,grep("-06-|-07-|-08-",index(roma_LST_monthly_ts)))
roma_LST_JLA_mean=calc(roma_LST_JLA@raster,function(x) mean(x,na.rm=T))
roma_LST_JLA_sd=calc(roma_LST_JLA@raster,function(x) sd(x,na.rm=T))
saveRDS(roma_LST_JLA_mean,"data/roma_LST_JLA_mean.rds")
saveRDS(roma_LST_JLA_sd,"data/roma_LST_JLA_sd.rds")
plotKML::kml(roma_LST_JLA_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(roma_LST_JLA_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


roma_LST_SON<-subset(roma_LST_monthly_ts,grep("-09-|-10-|-11-",index(roma_LST_monthly_ts)))
roma_LST_SON_mean=calc(roma_LST_SON@raster,function(x) mean(x,na.rm=T))
roma_LST_SON_sd=calc(roma_LST_SON@raster,function(x) sd(x,na.rm=T))
saveRDS(roma_LST_SON_mean,"data/roma_LST_SON_mean.rds")
saveRDS(roma_LST_SON_sd,"data/roma_LST_SON_sd.rds")
plotKML::kml(roma_LST_SON_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(roma_LST_SON_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


###############################################################################################
# bologna break and trend analisys NDVI MODIS 

bologna_istat=readRDS("data/bologna_istat.rds")

bologna_AATmap_ndvi=readRDS("data/bologna_AATmap_ndvi.rds")

bologna_slope_ndvi_b=bologna_AATmap_ndvi[["SlopeSEG1"]]  
bologna_break_ndvi_s_b<- TrendClassification(bologna_AATmap_ndvi[["PvalSEG1"]], min.length=10, max.pval=0.05)

bologna_slope_ndvi=rectify_matlab_array_mat(bologna_slope_ndvi_b,inCRS=proj_modis)
bologna_break_ndvi_s=rectify_matlab_array_mat(bologna_break_ndvi_s_b,inCRS=proj_modis)

bologna_slope_ndvi=crop_raster(bologna_istat,bologna_slope_ndvi)
bologna_break_ndvi=crop_raster(bologna_istat,bologna_break_ndvi_s)


plotKML::kml(bologna_slope_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 


plotKML::kml(bologna_break_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 



##########################################################################################
# bologna break and trend analisys LST MODIS Day

bologna_AATmap_LST=readRDS("data/bologna_AATmap_LST.rds")

bologna_slope_LST=bologna_AATmap_LST[["SlopeSEG1"]]  
bologna_trend_LST_s<- TrendClassification(bologna_AATmap_LST[["PvalSEG1"]], min.length=10, max.pval=0.1)

bologna_slope_LST=crop_raster(bologna_istat,bologna_slope_LST)
bologna_trend_LST_s<-crop_raster(bologna_istat,bologna_trend_LST_s)

plotKML::kml(bologna_slope_LST, colour_scale = SAGA_pal[[1]], colour = SlopeSEG1,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 
plotKML::kml(bologna_trend_LST_s, colour_scale = SAGA_pal[[1]], colour = TrendSEG1,width=25,height=200, z.lim = c(0,1),kmz=T) 


#######################################################################################################
# bologna LST climate 

bologna_LST_monthly_ts=readRDS("data/bologna_LST_monthly.rds")

bologna_LST_DJF<-subset(bologna_LST_monthly_ts,grep("-01-|-02-|-12-",index(bologna_LST_monthly_ts)))
bologna_LST_DJF_mean=calc(bologna_LST_DJF@raster,function(x) mean(x,na.rm=T))
bologna_LST_DJF_sd=calc(bologna_LST_DJF@raster,function(x) sd(x,na.rm=T))
saveRDS(bologna_LST_DJF_mean,"data/bologna_LST_DJF_mean.rds")
saveRDS(bologna_LST_DJF_sd,"data/bologna_LST_DJF_sd.rds")
plotKML::kml(bologna_LST_DJF_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(bologna_LST_DJF_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

             
bologna_LST_MAM<-subset(bologna_LST_monthly_ts,grep("-04-|-05-|-06-",index(bologna_LST_monthly_ts)))
bologna_LST_MAM_mean=calc(bologna_LST_MAM@raster,function(x) mean(x,na.rm=T))
bologna_LST_MAM_sd=calc(bologna_LST_MAM@raster,function(x) sd(x,na.rm=T))
saveRDS(bologna_LST_MAM_mean,"data/bologna_LST_MAM_mean.rds")
saveRDS(bologna_LST_MAM_sd,"data/bologna_LST_MAM_sd.rds")
plotKML::kml(bologna_LST_MAM_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(bologna_LST_MAM_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


bologna_LST_JLA<-subset(bologna_LST_monthly_ts,grep("-06-|-07-|-08-",index(bologna_LST_monthly_ts)))
bologna_LST_JLA_mean=calc(bologna_LST_JLA@raster,function(x) mean(x,na.rm=T))
bologna_LST_JLA_sd=calc(bologna_LST_JLA@raster,function(x) sd(x,na.rm=T))
saveRDS(bologna_LST_JLA_mean,"data/bologna_LST_JLA_mean.rds")
saveRDS(bologna_LST_JLA_sd,"data/bologna_LST_JLA_sd.rds")
plotKML::kml(bologna_LST_JLA_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(bologna_LST_JLA_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


bologna_LST_SON<-subset(bologna_LST_monthly_ts,grep("-09-|-10-|-11-",index(bologna_LST_monthly_ts)))
bologna_LST_SON_mean=calc(bologna_LST_SON@raster,function(x) mean(x,na.rm=T))
bologna_LST_SON_sd=calc(bologna_LST_SON@raster,function(x) sd(x,na.rm=T))
saveRDS(bologna_LST_SON_mean,"data/bologna_LST_SON_mean.rds")
saveRDS(bologna_LST_SON_sd,"data/bologna_LST_SON_sd.rds")
plotKML::kml(bologna_LST_SON_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(bologna_LST_SON_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

###############################################################################################
# firenze break and trend analisys NDVI MODIS 

firenze_istat=readRDS("data/firenze_istat.rds")

firenze_AATmap_ndvi=readRDS("data/firenze_AATmap_ndvi.rds")

firenze_slope_ndvi_b=firenze_AATmap_ndvi[["SlopeSEG1"]]  
firenze_break_ndvi_s_b<- TrendClassification(firenze_AATmap_ndvi[["PvalSEG1"]], min.length=10, max.pval=0.05)

firenze_slope_ndvi=rectify_matlab_array_mat(firenze_slope_ndvi_b,inCRS=proj_modis)
firenze_break_ndvi_s=rectify_matlab_array_mat(firenze_break_ndvi_s_b,inCRS=proj_modis)

firenze_slope_ndvi=crop_raster(firenze_istat,firenze_slope_ndvi)
firenze_break_ndvi=crop_raster(firenze_istat,firenze_break_ndvi_s)


plotKML::kml(firenze_slope_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 


plotKML::kml(firenze_break_ndvi_FI, colour_scale = SAGA_pal[[1]], colour = layer,png.width=1000,png.height=833,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 



##########################################################################################
# firenze break and trend analisys LST MODIS Day

firenze_AATmap_LST=readRDS("data/firenze_AATmap_LST.rds")

firenze_slope_LST=firenze_AATmap_LST[["SlopeSEG1"]]  
firenze_trend_LST_s<- TrendClassification(firenze_AATmap_LST[["PvalSEG1"]], min.length=10, max.pval=0.1)

firenze_slope_LST=crop_raster(firenze_istat,firenze_slope_LST)
firenze_trend_LST_s<-crop_raster(firenze_istat,firenze_trend_LST_s)

plotKML::kml(firenze_slope_LST, colour_scale = SAGA_pal[[1]], colour = SlopeSEG1,width=25,height=200, z.lim = c(-0.05,0.05), factor.labels=c("-0.05","0","0.05"),kmz=T) 
plotKML::kml(firenze_trend_LST_s, colour_scale = SAGA_pal[[1]], colour = TrendSEG1,width=25,height=200, z.lim = c(0,1),kmz=T) 


#######################################################################################################
# firenze LST climate 

firenze_LST_monthly_ts=readRDS("data/firenze_LST_monthly.rds")

firenze_LST_DJF<-subset(firenze_LST_monthly_ts,grep("-01-|-02-|-12-",index(firenze_LST_monthly_ts)))
firenze_LST_DJF_mean=calc(firenze_LST_DJF@raster,function(x) mean(x,na.rm=T))
firenze_LST_DJF_sd=calc(firenze_LST_DJF@raster,function(x) sd(x,na.rm=T))
saveRDS(firenze_LST_DJF_mean,"data/firenze_LST_DJF_mean.rds")
saveRDS(firenze_LST_DJF_sd,"data/firenze_LST_DJF_sd.rds")
plotKML::kml(firenze_LST_DJF_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(firenze_LST_DJF_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 

             
firenze_LST_MAM<-subset(firenze_LST_monthly_ts,grep("-04-|-05-|-06-",index(firenze_LST_monthly_ts)))
firenze_LST_MAM_mean=calc(firenze_LST_MAM@raster,function(x) mean(x,na.rm=T))
firenze_LST_MAM_sd=calc(firenze_LST_MAM@raster,function(x) sd(x,na.rm=T))
saveRDS(firenze_LST_MAM_mean,"data/firenze_LST_MAM_mean.rds")
saveRDS(firenze_LST_MAM_sd,"data/firenze_LST_MAM_sd.rds")
plotKML::kml(firenze_LST_MAM_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(firenze_LST_MAM_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


firenze_LST_JLA<-subset(firenze_LST_monthly_ts,grep("-06-|-07-|-08-",index(firenze_LST_monthly_ts)))
firenze_LST_JLA_mean=calc(firenze_LST_JLA@raster,function(x) mean(x,na.rm=T))
firenze_LST_JLA_sd=calc(firenze_LST_JLA@raster,function(x) sd(x,na.rm=T))
saveRDS(firenze_LST_JLA_mean,"data/firenze_LST_JLA_mean.rds")
saveRDS(firenze_LST_JLA_sd,"data/firenze_LST_JLA_sd.rds")
plotKML::kml(firenze_LST_JLA_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(firenze_LST_JLA_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


firenze_LST_SON<-subset(firenze_LST_monthly_ts,grep("-09-|-10-|-11-",index(firenze_LST_monthly_ts)))
firenze_LST_SON_mean=calc(firenze_LST_SON@raster,function(x) mean(x,na.rm=T))
firenze_LST_SON_sd=calc(firenze_LST_SON@raster,function(x) sd(x,na.rm=T))
saveRDS(firenze_LST_SON_mean,"data/firenze_LST_SON_mean.rds")
saveRDS(firenze_LST_SON_sd,"data/firenze_LST_SON_sd.rds")
plotKML::kml(firenze_LST_SON_mean, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(5,40), factor.labels=c("5","0","40"),kmz=T) 
plotKML::kml(firenze_LST_SON_sd, colour_scale = SAGA_pal[[1]], colour = layer,width=25,height=200, z.lim = c(0,5),kmz=T) 


##################################################################################


