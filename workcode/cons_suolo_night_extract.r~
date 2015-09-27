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

setwd("")

############################################################
# Load supplementary code

source("auxillary_functions.r")

##################################################################################

firenze_istat=readRDS("data/firenze_istat.rds")
bologna_istat=readRDS("data/bologna_istat.rds")
roma_istat=readRDS("data/roma_istat.rds")
milano_istat=readRDS("data/milano_istat.rds")
palermo_istat=readRDS("data/palermo_istat.rds")


firenze_AATmap_LST=readRDS("data/firenze_AATmap_LST_night.rds")
firenze_slope_LST=firenze_AATmap_LST[["SlopeSEG1"]]
firenze_slope_LST=crop_raster(firenze_istat,firenze_slope_LST)
firenze_LST_DJF_mean=readRDS("data/firenze_LST_DJF_mean_night.rds")
firenze_LST_SON_mean=readRDS("data/firenze_LST_SON_mean_night.rds")
firenze_LST_MAM_mean=readRDS("data/firenze_LST_MAM_mean_night.rds")
firenze_LST_JLA_mean=readRDS("data/firenze_LST_JLA_mean_night.rds")
firenze_LST_DJF_sd=readRDS("data/firenze_LST_DJF_sd_night.rds")
firenze_LST_SON_sd=readRDS("data/firenze_LST_SON_sd_night.rds")
firenze_LST_MAM_sd=readRDS("data/firenze_LST_MAM_sd_night.rds")
firenze_LST_JLA_sd=readRDS("data/firenze_LST_JLA_sd_night.rds")


writeRaster(firenze_slope_LST,"data/firenze_slope_LST_night.tif")
writeRaster(firenze_LST_DJF_mean,"data/firenze_LST_DJF_mean_night.tif")
writeRaster(firenze_LST_SON_mean,"data/firenze_LST_SON_mean_night.tif")
writeRaster(firenze_LST_JLA_mean,"data/firenze_LST_JLA_mean_night.tif")
writeRaster(firenze_LST_MAM_mean,"data/firenze_LST_MAM_mean_night.tif")
writeRaster(firenze_LST_DJF_sd,"data/firenze_LST_DJF_sd_night.tif")
writeRaster(firenze_LST_SON_sd,"data/firenze_LST_SON_sd_night.tif")
writeRaster(firenze_LST_JLA_sd,"data/firenze_LST_JLA_sd_night.tif")
writeRaster(firenze_LST_MAM_sd,"data/firenze_LST_MAM_sd_night.tif")

##########################################################################################################
milano_AATmap_LST=readRDS("data/milano_AATmap_LST_night.rds")
milano_slope_LST=milano_AATmap_LST[["SlopeSEG1"]]
milano_slope_LST=crop_raster(milano_istat,milano_slope_LST)
milano_LST_DJF_mean=readRDS("data/milano_LST_DJF_mean_night.rds")
milano_LST_SON_mean=readRDS("data/milano_LST_SON_mean_night.rds")
milano_LST_MAM_mean=readRDS("data/milano_LST_MAM_mean_night.rds")
milano_LST_JLA_mean=readRDS("data/milano_LST_JLA_mean_night.rds")
milano_LST_DJF_sd=readRDS("data/milano_LST_DJF_sd_night.rds")
milano_LST_SON_sd=readRDS("data/milano_LST_SON_sd_night.rds")
milano_LST_MAM_sd=readRDS("data/milano_LST_MAM_sd_night.rds")
milano_LST_JLA_sd=readRDS("data/milano_LST_JLA_sd_night.rds")


writeRaster(milano_slope_LST,"data/milano_slope_LST_night.tif")
writeRaster(milano_LST_DJF_mean,"data/milano_LST_DJF_mean_night.tif")
writeRaster(milano_LST_SON_mean,"data/milano_LST_SON_mean_night.tif")
writeRaster(milano_LST_JLA_mean,"data/milano_LST_JLA_mean_night.tif")
writeRaster(milano_LST_MAM_mean,"data/milano_LST_MAM_mean_night.tif")
writeRaster(milano_LST_DJF_sd,"data/milano_LST_DJF_sd_night.tif")
writeRaster(milano_LST_SON_sd,"data/milano_LST_SON_sd_night.tif")
writeRaster(milano_LST_JLA_sd,"data/milano_LST_JLA_sd_night.tif")
writeRaster(milano_LST_MAM_sd,"data/milano_LST_MAM_sd_night.tif")


################################################################################################

roma_AATmap_LST=readRDS("data/roma_AATmap_LST_night.rds")
roma_slope_LST=roma_AATmap_LST[["SlopeSEG1"]]
roma_slope_LST=crop_raster(roma_istat,roma_slope_LST)
roma_LST_DJF_mean=readRDS("data/roma_LST_DJF_mean_night.rds")
roma_LST_SON_mean=readRDS("data/roma_LST_SON_mean_night.rds")
roma_LST_MAM_mean=readRDS("data/roma_LST_MAM_mean_night.rds")
roma_LST_JLA_mean=readRDS("data/roma_LST_JLA_mean_night.rds")
roma_LST_DJF_sd=readRDS("data/roma_LST_DJF_sd_night.rds")
roma_LST_SON_sd=readRDS("data/roma_LST_SON_sd_night.rds")
roma_LST_MAM_sd=readRDS("data/roma_LST_MAM_sd_night.rds")
roma_LST_JLA_sd=readRDS("data/roma_LST_JLA_sd_night.rds")

writeRaster(roma_slope_LST,"data/roma_slope_LST_night.tif")
writeRaster(roma_LST_DJF_mean,"data/roma_LST_DJF_mean_night.tif")
writeRaster(roma_LST_SON_mean,"data/roma_LST_SON_mean_night.tif")
writeRaster(roma_LST_JLA_mean,"data/roma_LST_JLA_mean_night.tif")
writeRaster(roma_LST_MAM_mean,"data/roma_LST_MAM_mean_night.tif")
writeRaster(roma_LST_DJF_sd,"data/roma_LST_DJF_sd_night.tif")
writeRaster(roma_LST_SON_sd,"data/roma_LST_SON_sd_night.tif")
writeRaster(roma_LST_JLA_sd,"data/roma_LST_JLA_sd_night.tif")
writeRaster(roma_LST_MAM_sd,"data/roma_LST_MAM_sd_night.tif")

################################################################################################
bologna_AATmap_LST=readRDS("data/bologna_AATmap_LST_night.rds")
bologna_slope_LST=bologna_AATmap_LST[["SlopeSEG1"]]
bologna_slope_LST=crop_raster(bologna_istat,bologna_slope_LST)
bologna_LST_DJF_mean=readRDS("data/bologna_LST_DJF_mean_night.rds")
bologna_LST_SON_mean=readRDS("data/bologna_LST_SON_mean_night.rds")
bologna_LST_MAM_mean=readRDS("data/bologna_LST_MAM_mean_night.rds")
bologna_LST_JLA_mean=readRDS("data/bologna_LST_JLA_mean_night.rds")
bologna_LST_DJF_sd=readRDS("data/bologna_LST_DJF_sd_night.rds")
bologna_LST_SON_sd=readRDS("data/bologna_LST_SON_sd_night.rds")
bologna_LST_MAM_sd=readRDS("data/bologna_LST_MAM_sd_night.rds")
bologna_LST_JLA_sd=readRDS("data/bologna_LST_JLA_sd_night.rds")



writeRaster(bologna_slope_LST,"data/bologna_slope_LST_night.tif")
writeRaster(bologna_LST_DJF_mean,"data/bologna_LST_DJF_mean_night.tif")
writeRaster(bologna_LST_SON_mean,"data/bologna_LST_SON_mean_night.tif")
writeRaster(bologna_LST_JLA_mean,"data/bologna_LST_JLA_mean_night.tif")
writeRaster(bologna_LST_MAM_mean,"data/bologna_LST_MAM_mean_night.tif")
writeRaster(bologna_LST_DJF_sd,"data/bologna_LST_DJF_sd_night.tif")
writeRaster(bologna_LST_SON_sd,"data/bologna_LST_SON_sd_night.tif")
writeRaster(bologna_LST_JLA_sd,"data/bologna_LST_JLA_sd_night.tif")
writeRaster(bologna_LST_MAM_sd,"data/bologna_LST_MAM_sd_night.tif")

################################################################################################
palermo_AATmap_LST=readRDS("data/palermo_AATmap_LST_night.rds")
palermo_slope_LST=palermo_AATmap_LST[["SlopeSEG1"]]
palermo_slope_LST=crop_raster(palermo_istat,palermo_slope_LST)
palermo_LST_DJF_mean=readRDS("data/palermo_LST_DJF_mean_night.rds")
palermo_LST_SON_mean=readRDS("data/palermo_LST_SON_mean_night.rds")
palermo_LST_MAM_mean=readRDS("data/palermo_LST_MAM_mean_night.rds")
palermo_LST_JLA_mean=readRDS("data/palermo_LST_JLA_mean_night.rds")
palermo_LST_DJF_sd=readRDS("data/palermo_LST_DJF_sd_night.rds")
palermo_LST_SON_sd=readRDS("data/palermo_LST_SON_sd_night.rds")
palermo_LST_MAM_sd=readRDS("data/palermo_LST_MAM_sd_night.rds")
palermo_LST_JLA_sd=readRDS("data/palermo_LST_JLA_sd_night.rds")

writeRaster(palermo_slope_LST,"data/palermo_slope_LST_night.tif")
writeRaster(palermo_LST_DJF_mean,"data/palermo_LST_DJF_mean_night.tif")
writeRaster(palermo_LST_SON_mean,"data/palermo_LST_SON_mean_night.tif")
writeRaster(palermo_LST_JLA_mean,"data/palermo_LST_JLA_mean_night.tif")
writeRaster(palermo_LST_MAM_mean,"data/palermo_LST_MAM_mean_night.tif")
writeRaster(palermo_LST_DJF_sd,"data/palermo_LST_DJF_sd_night.tif")
writeRaster(palermo_LST_SON_sd,"data/palermo_LST_SON_sd_night.tif")
writeRaster(palermo_LST_JLA_sd,"data/palermo_LST_JLA_sd_night.tif")
writeRaster(palermo_LST_MAM_sd,"data/palermo_LST_MAM_sd_night.tif")

################################################################################################

system("python rasterstats firenze_LST_poly.shp data/firenze_slope_LST_night.tif >firenze_slope_LST_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_DJF_mean_night.tif >firenze_LST_DJF_mean_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_SON_mean_night.tif >firenze_LST_SON_mean_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_MAM_mean_night.tif >firenze_LST_MAM_mean_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_JLA_mean_night.tif >firenze_LST_JLA_mean_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_DJF_sd_night.tif >firenze_LST_DJF_sd_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_SON_sd_night.tif >firenze_LST_SON_sd_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_MAM_sd_night.tif >firenze_LST_MAM_sd_night_stats.csv")
system("python rasterstats firenze_LST_poly.shp data/firenze_LST_JLA_sd_night.tif >firenze_LST_JLA_sd_night_stats.csv")


system("python rasterstats roma_LST_poly.shp data/roma_slope_LST_night.tif >roma_slope_LST_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_DJF_mean_night.tif >roma_LST_DJF_mean_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_SON_mean_night.tif >roma_LST_SON_mean_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_MAM_mean_night.tif >roma_LST_MAM_mean_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_JLA_mean_night.tif >roma_LST_JLA_mean_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_DJF_sd_night.tif >roma_LST_DJF_sd_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_SON_sd_night.tif >roma_LST_SON_sd_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_MAM_sd_night.tif >roma_LST_MAM_sd_night_stats.csv")
system("python rasterstats roma_LST_poly.shp data/roma_LST_JLA_sd_night.tif >roma_LST_JLA_sd_night_stats.csv")

system("python rasterstats milano_LST_poly.shp data/milano_slope_LST_night.tif >milano_slope_LST_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_DJF_mean_night.tif >milano_LST_DJF_mean_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_SON_mean_night.tif >milano_LST_SON_mean_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_MAM_mean_night.tif >milano_LST_MAM_mean_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_JLA_mean_night.tif >milano_LST_JLA_mean_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_DJF_sd_night.tif >milano_LST_DJF_sd_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_SON_sd_night.tif >milano_LST_SON_sd_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_MAM_sd_night.tif >milano_LST_MAM_sd_night_stats.csv")
system("python rasterstats milano_LST_poly.shp data/milano_LST_JLA_sd_night.tif >milano_LST_JLA_sd_night_stats.csv")


system("python rasterstats bologna_LST_poly.shp data/bologna_slope_LST_night.tif >bologna_slope_LST_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_DJF_mean_night.tif >bologna_LST_DJF_mean_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_SON_mean_night.tif >bologna_LST_SON_mean_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_MAM_mean_night.tif >bologna_LST_MAM_mean_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_JLA_mean_night.tif >bologna_LST_JLA_mean_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_DJF_sd_night.tif >bologna_LST_DJF_sd_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_SON_sd_night.tif >bologna_LST_SON_sd_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_MAM_sd_night.tif >bologna_LST_MAM_sd_night_stats.csv")
system("python rasterstats bologna_LST_poly.shp data/bologna_LST_JLA_sd_night.tif >bologna_LST_JLA_sd_night_stats.csv")


system("python rasterstats palermo_LST_poly.shp data/palermo_slope_LST_night.tif >palermo_slope_LST_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_DJF_mean_night.tif >palermo_LST_DJF_mean_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_SON_mean_night.tif >palermo_LST_SON_mean_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_MAM_mean_night.tif >palermo_LST_MAM_mean_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_JLA_mean_night.tif >palermo_LST_JLA_mean_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_DJF_sd_night.tif >palermo_LST_DJF_sd_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_SON_sd_night.tif >palermo_LST_SON_sd_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_MAM_sd_night.tif >palermo_LST_MAM_sd_night_stats.csv")
system("python rasterstats palermo_LST_poly.shp data/palermo_LST_JLA_sd_night.tif >palermo_LST_JLA_sd_night_stats.csv")

#######################################################################################################################





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
