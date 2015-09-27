############################################################

library(raster)
library(rts)
library(rasterVis)
library(ggplot2)
library(rgdal)
library(grid)
library(latticeExtra)
library(maptools)
library(plyr)
library(plotKML)

# require rasterstats python

############################################################
# Setup directory

setwd("/home/alf/Scrivania/lav_morab_LST")

############################################################
# Load supplementary code

source("auxillary_functions.r")




firenze_istat=readRDS("../r_data/firenze_istat.rds")
bologna_istat=readRDS("../r_data/bologna_istat.rds")
roma_istat=readRDS("../r_data/roma_istat.rds")
milano_istat=readRDS("../r_data/milano_istat.rds")
palermo_istat=readRDS("../r_data/palermo_istat.rds")

writeOGR(milano_istat, ".", "milano_istat", driver="ESRI Shapefile")
system("gdalwarp  -cutline milano_istat.shp data/Milano_ispra_csuolo.tif data/milano_ispra_csuolo_clip.tif")
milano_csuolo=milano_csuolo==1
writeRaster(milano_csuolo,"data/milano_ispra_csuolo_clip.tif",overwrite=TRUE)
milano_csuolo=raster("../r_data/milano_ispra_csuolo_clip.tif")

saveRDS(milano_csuolo,"milano_csuolo.rds")

writeOGR(palermo_istat, ".", "palermo_istat", driver="ESRI Shapefile")
system("gdalwarp  -cutline palermo_istat.shp data/Palermo_ispra_csuolo.tif data/palermo_ispra_csuolo_clip.tif")
palermo_csuolo=raster("../r_data/palermo_ispra_csuolo_clip.tif")
palermo_csuolo=palermo_csuolo==1
saveRDS(palermo_csuolo,"palermo_csuolo.rds")


writeOGR(bologna_istat, ".", "bologna_istat", driver="ESRI Shapefile")
system("gdalwarp  -cutline bologna_istat.shp data/Bologna_ispra_csuolo.tif data/bologna_ispra_csuolo_clip.tif")
bologna_csuolo=raster("../r_data/bologna_ispra_csuolo_clip.tif")
bologna_csuolo=bologna_csuolo==1
saveRDS(bologna_csuolo,"bologna_csuolo.rds")


writeOGR(roma_istat, ".", "roma_istat", driver="ESRI Shapefile")
system("gdalwarp  -cutline roma_istat.shp data/Roma_ispra_csuolo.tif data/Roma_ispra_csuolo_clip.tif")
roma_csuolo=raster("../r_data/Roma_ispra_csuolo_clip.tif")
roma_csuolo=roma_csuolo==1
saveRDS(roma_csuolo,"roma_csuolo.rds")


writeOGR(firenze_istat, ".", "firenze_istat", driver="ESRI Shapefile")
system("gdalwarp  -cutline firenze_istat.shp data/Firenze_ispra_csuolo.tif data/firenze_ispra_csuolo_clip.tif")
firenze_csuolo=raster("../r_data/firenze_ispra_csuolo_clip.tif")
firenze_csuolo=firenze_csuolo==1
saveRDS(firenze_csuolo,"firenze_csuolo.rds")
##################################################################################

##################################################################################

firenze_csuolo=readRDS("firenze_csuolo.rds")
milano_csuolo=readRDS("milano_csuolo.rds")
palermo_csuolo=readRDS("palermo_csuolo.rds")
bologna_csuolo=readRDS("bologna_csuolo.rds")
roma_csuolo=readRDS("roma_csuolo.rds")


##################################################################################


firenze_AATmap_LST=readRDS("../r_data/firenze_AATmap_LST.rds")
firenze_slope_LST=firenze_AATmap_LST[["SlopeSEG1"]]
firenze_slope_LST=crop_raster(firenze_istat,firenze_slope_LST)
firenze_LST_DJF_mean=readRDS("../r_data/firenze_LST_DJF_mean.rds")
firenze_LST_SON_mean=readRDS("../r_data/firenze_LST_SON_mean.rds")
firenze_LST_MAM_mean=readRDS("../r_data/firenze_LST_MAM_mean.rds")
firenze_LST_JLA_mean=readRDS("../r_data/firenze_LST_JLA_mean.rds")
firenze_LST_DJF_sd=readRDS("../r_data/firenze_LST_DJF_sd.rds")
firenze_LST_SON_sd=readRDS("../r_data/firenze_LST_SON_sd.rds")
firenze_LST_MAM_sd=readRDS("../r_data/firenze_LST_MAM_sd.rds")
firenze_LST_JLA_sd=readRDS("../r_data/firenze_LST_JLA_sd.rds")

firenze_LST_poly=rasterToPolygons(firenze_slope_LST)

saveRDS(firenze_LST_poly,"firenze_LST_poly.rds")
writeOGR(firenze_LST_poly, ".", "firenze_LST_poly", driver="ESRI Shapefile")


writeRaster(firenze_slope_LST,"data/firenze_slope_LST.tif")
writeRaster(firenze_LST_DJF_mean,"data/firenze_LST_DJF_mean.tif")
writeRaster(firenze_LST_SON_mean,"data/firenze_LST_SON_mean.tif")
writeRaster(firenze_LST_JLA_mean,"data/firenze_LST_JLA_mean.tif")
writeRaster(firenze_LST_MAM_mean,"data/firenze_LST_MAM_mean.tif")
writeRaster(firenze_LST_DJF_sd,"data/firenze_LST_DJF_sd.tif")
writeRaster(firenze_LST_SON_sd,"data/firenze_LST_SON_sd.tif")
writeRaster(firenze_LST_JLA_sd,"data/firenze_LST_JLA_sd.tif")
writeRaster(firenze_LST_MAM_sd,"data/firenze_LST_MAM_sd.tif")

##########################################################################################################
milano_AATmap_LST=readRDS("../r_data/milano_AATmap_LST.rds")
milano_slope_LST=milano_AATmap_LST[["SlopeSEG1"]]
milano_slope_LST=crop_raster(milano_istat,milano_slope_LST)
milano_LST_DJF_mean=readRDS("../r_data/milano_LST_DJF_mean.rds")
milano_LST_SON_mean=readRDS("../r_data/milano_LST_SON_mean.rds")
milano_LST_MAM_mean=readRDS("../r_data/milano_LST_MAM_mean.rds")
milano_LST_JLA_mean=readRDS("../r_data/milano_LST_JLA_mean.rds")
milano_LST_DJF_sd=readRDS("../r_data/milano_LST_DJF_sd.rds")
milano_LST_SON_sd=readRDS("../r_data/milano_LST_SON_sd.rds")
milano_LST_MAM_sd=readRDS("../r_data/milano_LST_MAM_sd.rds")
milano_LST_JLA_sd=readRDS("../r_data/milano_LST_JLA_sd.rds")

milano_LST_poly=rasterToPolygons(milano_slope_LST)

saveRDS(milano_LST_poly,"milano_LST_poly.rds")
writeOGR(milano_LST_poly, ".", "milano_LST_poly", driver="ESRI Shapefile")


writeRaster(milano_slope_LST,"data/milano_slope_LST.tif")
writeRaster(milano_LST_DJF_mean,"data/milano_LST_DJF_mean.tif")
writeRaster(milano_LST_SON_mean,"data/milano_LST_SON_mean.tif")
writeRaster(milano_LST_JLA_mean,"data/milano_LST_JLA_mean.tif")
writeRaster(milano_LST_MAM_mean,"data/milano_LST_MAM_mean.tif")
writeRaster(milano_LST_DJF_sd,"data/milano_LST_DJF_sd.tif")
writeRaster(milano_LST_SON_sd,"data/milano_LST_SON_sd.tif")
writeRaster(milano_LST_JLA_sd,"data/milano_LST_JLA_sd.tif")
writeRaster(milano_LST_MAM_sd,"data/milano_LST_MAM_sd.tif")


################################################################################################

roma_AATmap_LST=readRDS("../r_data/roma_AATmap_LST.rds")
roma_slope_LST=roma_AATmap_LST[["SlopeSEG1"]]
roma_slope_LST=crop_raster(roma_istat,roma_slope_LST)
roma_LST_DJF_mean=readRDS("../r_data/roma_LST_DJF_mean.rds")
roma_LST_SON_mean=readRDS("../r_data/roma_LST_SON_mean.rds")
roma_LST_MAM_mean=readRDS("../r_data/roma_LST_MAM_mean.rds")
roma_LST_JLA_mean=readRDS("../r_data/roma_LST_JLA_mean.rds")
roma_LST_DJF_sd=readRDS("../r_data/roma_LST_DJF_sd.rds")
roma_LST_SON_sd=readRDS("../r_data/roma_LST_SON_sd.rds")
roma_LST_MAM_sd=readRDS("../r_data/roma_LST_MAM_sd.rds")
roma_LST_JLA_sd=readRDS("../r_data/roma_LST_JLA_sd.rds")

roma_LST_poly=rasterToPolygons(roma_slope_LST)

saveRDS(roma_LST_poly,"roma_LST_poly.rds")
writeOGR(roma_LST_poly, ".", "roma_LST_poly", driver="ESRI Shapefile")


writeRaster(roma_slope_LST,"data/roma_slope_LST.tif")
writeRaster(roma_LST_DJF_mean,"data/roma_LST_DJF_mean.tif")
writeRaster(roma_LST_SON_mean,"data/roma_LST_SON_mean.tif")
writeRaster(roma_LST_JLA_mean,"data/roma_LST_JLA_mean.tif")
writeRaster(roma_LST_MAM_mean,"data/roma_LST_MAM_mean.tif")
writeRaster(roma_LST_DJF_sd,"data/roma_LST_DJF_sd.tif")
writeRaster(roma_LST_SON_sd,"data/roma_LST_SON_sd.tif")
writeRaster(roma_LST_JLA_sd,"data/roma_LST_JLA_sd.tif")
writeRaster(roma_LST_MAM_sd,"data/roma_LST_MAM_sd.tif")

################################################################################################
bologna_AATmap_LST=readRDS("../r_data/bologna_AATmap_LST.rds")
bologna_slope_LST=bologna_AATmap_LST[["SlopeSEG1"]]
bologna_slope_LST=crop_raster(bologna_istat,bologna_slope_LST)
bologna_LST_DJF_mean=readRDS("../r_data/bologna_LST_DJF_mean.rds")
bologna_LST_SON_mean=readRDS("../r_data/bologna_LST_SON_mean.rds")
bologna_LST_MAM_mean=readRDS("../r_data/bologna_LST_MAM_mean.rds")
bologna_LST_JLA_mean=readRDS("../r_data/bologna_LST_JLA_mean.rds")
bologna_LST_DJF_sd=readRDS("../r_data/bologna_LST_DJF_sd.rds")
bologna_LST_SON_sd=readRDS("../r_data/bologna_LST_SON_sd.rds")
bologna_LST_MAM_sd=readRDS("../r_data/bologna_LST_MAM_sd.rds")
bologna_LST_JLA_sd=readRDS("../r_data/bologna_LST_JLA_sd.rds")

bologna_LST_poly=rasterToPolygons(bologna_slope_LST)

saveRDS(bologna_LST_poly,"bologna_LST_poly.rds")
writeOGR(bologna_LST_poly, ".", "bologna_LST_poly", driver="ESRI Shapefile")


writeRaster(bologna_slope_LST,"data/bologna_slope_LST.tif")
writeRaster(bologna_LST_DJF_mean,"data/bologna_LST_DJF_mean.tif")
writeRaster(bologna_LST_SON_mean,"data/bologna_LST_SON_mean.tif")
writeRaster(bologna_LST_JLA_mean,"data/bologna_LST_JLA_mean.tif")
writeRaster(bologna_LST_MAM_mean,"data/bologna_LST_MAM_mean.tif")
writeRaster(bologna_LST_DJF_sd,"data/bologna_LST_DJF_sd.tif")
writeRaster(bologna_LST_SON_sd,"data/bologna_LST_SON_sd.tif")
writeRaster(bologna_LST_JLA_sd,"data/bologna_LST_JLA_sd.tif")
writeRaster(bologna_LST_MAM_sd,"data/bologna_LST_MAM_sd.tif")

################################################################################################
palermo_AATmap_LST=readRDS("../r_data/palermo_AATmap_LST.rds")
palermo_slope_LST=palermo_AATmap_LST[["SlopeSEG1"]]
palermo_slope_LST=crop_raster(palermo_istat,palermo_slope_LST)
palermo_LST_DJF_mean=readRDS("../r_data/palermo_LST_DJF_mean.rds")
palermo_LST_SON_mean=readRDS("../r_data/palermo_LST_SON_mean.rds")
palermo_LST_MAM_mean=readRDS("../r_data/palermo_LST_MAM_mean.rds")
palermo_LST_JLA_mean=readRDS("../r_data/palermo_LST_JLA_mean.rds")
palermo_LST_DJF_sd=readRDS("../r_data/palermo_LST_DJF_sd.rds")
palermo_LST_SON_sd=readRDS("../r_data/palermo_LST_SON_sd.rds")
palermo_LST_MAM_sd=readRDS("../r_data/palermo_LST_MAM_sd.rds")
palermo_LST_JLA_sd=readRDS("../r_data/palermo_LST_JLA_sd.rds")

palermo_LST_poly=rasterToPolygons(palermo_slope_LST)

saveRDS(palermo_LST_poly,"palermo_LST_poly.rds")
writeOGR(palermo_LST_poly, ".", "palermo_LST_poly", driver="ESRI Shapefile")


writeRaster(palermo_slope_LST,"data/palermo_slope_LST.tif")
writeRaster(palermo_LST_DJF_mean,"data/palermo_LST_DJF_mean.tif")
writeRaster(palermo_LST_SON_mean,"data/palermo_LST_SON_mean.tif")
writeRaster(palermo_LST_JLA_mean,"data/palermo_LST_JLA_mean.tif")
writeRaster(palermo_LST_MAM_mean,"data/palermo_LST_MAM_mean.tif")
writeRaster(palermo_LST_DJF_sd,"data/palermo_LST_DJF_sd.tif")
writeRaster(palermo_LST_SON_sd,"data/palermo_LST_SON_sd.tif")
writeRaster(palermo_LST_JLA_sd,"data/palermo_LST_JLA_sd.tif")
writeRaster(palermo_LST_MAM_sd,"data/palermo_LST_MAM_sd.tif")

################################################################################################

system("python rasterstats firenze_LST_poly.shp data/firenze_slope_LST.tif >firenze_slope_LST_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_DJF_mean.tif >firenze_LST_DJF_mean_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_SON_mean.tif >firenze_LST_SON_mean_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_MAM_mean.tif >firenze_LST_MAM_mean_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_JLA_mean.tif >firenze_LST_JLA_mean_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_DJF_sd.tif >firenze_LST_DJF_sd_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_SON_sd.tif >firenze_LST_SON_sd_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_MAM_sd.tif >firenze_LST_MAM_sd_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_JLA_sd.tif >firenze_LST_JLA_sd_stats.csv")


system("python rasterstats roma_LST_poly.shp data/roma_slope_LST.tif >roma_slope_LST_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_DJF_mean.tif >roma_LST_DJF_mean_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_SON_mean.tif >roma_LST_SON_mean_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_MAM_mean.tif >roma_LST_MAM_mean_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_JLA_mean.tif >roma_LST_JLA_mean_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_DJF_sd.tif >roma_LST_DJF_sd_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_SON_sd.tif >roma_LST_SON_sd_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_MAM_sd.tif >roma_LST_MAM_sd_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_JLA_sd.tif >roma_LST_JLA_sd_stats.csv")

system("python rasterstats milano_LST_poly.shp data/milano_slope_LST.tif >milano_slope_LST_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_DJF_mean.tif >milano_LST_DJF_mean_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_SON_mean.tif >milano_LST_SON_mean_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_MAM_mean.tif >milano_LST_MAM_mean_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_JLA_mean.tif >milano_LST_JLA_mean_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_DJF_sd.tif >milano_LST_DJF_sd_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_SON_sd.tif >milano_LST_SON_sd_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_MAM_sd.tif >milano_LST_MAM_sd_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_JLA_sd.tif >milano_LST_JLA_sd_stats.csv")


system("python rasterstats bologna_LST_poly.shp data/bologna_slope_LST.tif >bologna_slope_LST_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_DJF_mean.tif >bologna_LST_DJF_mean_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_SON_mean.tif >bologna_LST_SON_mean_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_MAM_mean.tif >bologna_LST_MAM_mean_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_JLA_mean.tif >bologna_LST_JLA_mean_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_DJF_sd.tif >bologna_LST_DJF_sd_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_SON_sd.tif >bologna_LST_SON_sd_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_MAM_sd.tif >bologna_LST_MAM_sd_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_JLA_sd.tif >bologna_LST_JLA_sd_stats.csv")


system("python rasterstats palermo_LST_poly.shp data/palermo_slope_LST.tif >palermo_slope_LST_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_DJF_mean.tif >palermo_LST_DJF_mean_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_SON_mean.tif >palermo_LST_SON_mean_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_MAM_mean.tif >palermo_LST_MAM_mean_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_JLA_mean.tif >palermo_LST_JLA_mean_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_DJF_sd.tif >palermo_LST_DJF_sd_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_SON_sd.tif >palermo_LST_SON_sd_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_MAM_sd.tif >palermo_LST_MAM_sd_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_JLA_sd.tif >palermo_LST_JLA_sd_stats.csv")

#######################################################################################################################
system("python rasterstats data/palermo_LST_poly_32n.shp data/palermo_ispra_csuolo_clip_32n.tif >palermo_csuolo_stats.csv")
system("python rasterstats data/roma_LST_poly_32n.shp data/roma_ispra_csuolo_clip_32n.tif >roma_csuolo_stats.csv")
system("python rasterstats data/milano_LST_poly_32n.shp data/milano_ispra_csuolo_clip_32n.tif >milano_csuolo_stats.csv")
system("python rasterstats data/bologna_LST_poly_32n.shp data/bologna_ispra_csuolo_clip_32n.tif >bologna_csuolo_stats.csv")
system("python rasterstats data/firenze_LST_poly_32n.shp data/firenze_ispra_csuolo_clip_32n.tif >firenze_csuolo_stats.csv")
#######################################################################################################################



# firenze_ex_suolo <- extract(firenze_csuolo,firenze_LST_poly, fun = sum, na.rm = TRUE)
# milano_ex_suolo <- extract(milano_csuolo,milano_LST_poly, fun = sum, na.rm = TRUE)
# roma_ex_suolo <- extract(roma_csuolo,roma_LST_poly, fun = sum, na.rm = TRUE)
# bologna_ex_suolo <- extract(bologna_csuolo,bologna_LST_poly, fun = sum, na.rm = TRUE)
# palermo_ex_suolo <- extract(palermo_csuolo,palermo_LST_poly, fun = sum, na.rm = TRUE)
# 
# saveRDS(firenze_ex_suolo,"firenze_ex_suolo_s.rds")
# saveRDS(milano_ex_suolo,"milano_ex_suolo_s.rds")
# saveRDS(roma_ex_suolo,"roma_ex_suolo_s.rds")
# saveRDS(roma_ex_suolo,"bologna_ex_suolo_s.rds")
# saveRDS(roma_ex_suolo,"palermo_ex_suolo_s.rds")


################################################################################
# Names
# [1] "roma_csuolo_stats"           "roma_LST_DJF_mean_stats"    
# [3] "roma_LST_DJF_sd_stats"       "roma_LST_JLA_SON_mean_stats"
# [5] "roma_LST_JLA_SON_sd_stats"   "roma_LST_MAM_mean_stats"    
# [7] "roma_LST_MAM_sd_stats"       "roma_LST_SON_mean_stats"    
# [9] "roma_LST_SON_sd_stats"       "roma_slope_LST_stats"       

# [1] "X__fid__" "count"    "majority" "max"      "mean"     "median"   "min"     
# [8] "minority" "range"    "std"      "sum"      "unique"  




###########################################################################
# References 
# https://github.com/dwtkns/gdal-cheat-sheet
