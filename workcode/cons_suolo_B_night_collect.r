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
library(robust)

############################################################
# Setup directory

setwd("")

############################################################
# Load supplementary code

source("auxillary_functions.r")
file.copy(Sys.glob("*_csuolo*_stats.csv"),gsub("_stats.csv","_night_stats.csv",Sys.glob("*_csuolo*_stats.csv")))
stat_roma_list=lapply(Sys.glob("roma*_night_stats.csv"),read.csv)
names(stat_roma_list)<-gsub(".csv","",Sys.glob("roma*_night_stats.csv"))


stat_palermo_list=lapply(Sys.glob("palermo*_night_stats.csv"),read.csv)
names(stat_palermo_list)<-gsub(".csv","",Sys.glob("palermo*_night_stats.csv"))

stat_firenze_list=lapply(Sys.glob("firenze*_night_stats.csv"),read.csv)
names(stat_firenze_list)<-gsub(".csv","",Sys.glob("firenze*_night_stats.csv"))

stat_bologna_list=lapply(Sys.glob("bologna*_night_stats.csv"),read.csv)
names(stat_bologna_list)<-gsub(".csv","",Sys.glob("bologna*_night_stats.csv"))

stat_milano_list=lapply(Sys.glob("milano*_night_stats.csv"),read.csv)
names(stat_milano_list)<-gsub(".csv","",Sys.glob("milano*_night_stats.csv"))

# [1] "roma_csuolo_stats"           "roma_LST_DJF_mean_stats"    
# [3] "roma_LST_DJF_sd_stats"       "roma_LST_JLA_SON_mean_stats"
# [5] "roma_LST_JLA_SON_sd_stats"   "roma_LST_MAM_mean_stats"    
# [7] "roma_LST_MAM_sd_stats"       "roma_LST_SON_mean_stats"    
# [9] "roma_LST_SON_sd_stats"       "roma_slope_LST_stats"       

# [1] "X__fid__" "count"    "majority" "max"      "mean"     "median"   "min"     
# [8] "minority" "range"    "std"      "sum"      "unique"  


roma_df=data.frame(roma_sum_cs=stat_roma_list$roma_csuolo_night_stats$sum,
                   roma_mean_cs=stat_roma_list$roma_csuolo_night_stats$mean,
                   roma_std_cs=stat_roma_list$roma_csuolo_night_stats$std,
                   roma_LST_DJF_mean=stat_roma_list[[2]]$mean,
                   roma_LST_DJF_sd=stat_roma_list[[3]]$mean,
                   roma_LST_JLA_mean=stat_roma_list[[4]]$mean,
                   roma_LST_JLA_sd=stat_roma_list[[5]]$mean,
                   roma_LST_MAM_mean=stat_roma_list[[6]]$mean,
                   roma_LST_MAM_sd=stat_roma_list[[7]]$mean,
                   roma_LST_SON_mean=stat_roma_list[[8]]$mean,
                   roma_LST_SON_sd=stat_roma_list[[9]]$mean,
                   roma_LST_slope_sd=stat_roma_list[[10]]$mean
                   )
write.csv(roma_df,"roma_df_night.xls",row.names = FALSE)

milano_df=data.frame(milano_sum_cs=stat_milano_list$milano_csuolo_night_stats$sum,
                     milano_mean_cs=stat_milano_list$milano_csuolo_night_stats$mean,
                     milano_std_cs=stat_milano_list$milano_csuolo_night_stats$std,
                     milano_LST_DJF_mean=stat_milano_list[[2]]$mean,
                     milano_LST_DJF_sd=stat_milano_list[[3]]$mean,
                     milano_LST_JLA_mean=stat_milano_list[[4]]$mean,
                     milano_LST_JLA_sd=stat_milano_list[[5]]$mean,
                     milano_LST_MAM_mean=stat_milano_list[[6]]$mean,
                     milano_LST_MAM_sd=stat_milano_list[[7]]$mean,
                     milano_LST_SON_mean=stat_milano_list[[8]]$mean,
                     milano_LST_SON_sd=stat_milano_list[[9]]$mean,
                     milano_LST_slope_sd=stat_milano_list[[10]]$mean
)
write.csv(milano_df,"milano_df_night.xls",row.names = FALSE)


firenze_df=data.frame(firenze_sum_cs=stat_firenze_list$firenze_csuolo_night_stats$sum,
                      firenze_mean_cs=stat_firenze_list$firenze_csuolo_night_stats$mean,
                       firenze_std_cs=stat_firenze_list$firenze_csuolo_night_stats$std,
                      firenze_LST_DJF_mean=stat_firenze_list[[2]]$mean,
                      firenze_LST_DJF_sd=stat_firenze_list[[3]]$mean,
                      firenze_LST_JLA_mean=stat_firenze_list[[4]]$mean,
                      firenze_LST_JLA_sd=stat_firenze_list[[5]]$mean,
                      firenze_LST_MAM_mean=stat_firenze_list[[6]]$mean,
                      firenze_LST_MAM_sd=stat_firenze_list[[7]]$mean,
                      firenze_LST_SON_mean=stat_firenze_list[[8]]$mean,
                      firenze_LST_SON_sd=stat_firenze_list[[9]]$mean,
                      firenze_LST_slope_sd=stat_firenze_list[[10]]$mean
)
write.csv(firenze_df,"firenze_df_night.xls",row.names = FALSE)


bologna_df=data.frame(bologna_sum_cs=stat_bologna_list$bologna_csuolo_night_stats$sum,
                      bologna_mean_cs=stat_bologna_list$bologna_csuolo_night_stats$mean,
                      bologna_std_cs=stat_bologna_list$bologna_csuolo_night_stats$std,
                      bologna_LST_DJF_mean=stat_bologna_list[[2]]$mean,
                      bologna_LST_DJF_sd=stat_bologna_list[[3]]$mean,
                      bologna_LST_JLA_mean=stat_bologna_list[[4]]$mean,
                      bologna_LST_JLA_sd=stat_bologna_list[[5]]$mean,
                      bologna_LST_MAM_mean=stat_bologna_list[[6]]$mean,
                      bologna_LST_MAM_sd=stat_bologna_list[[7]]$mean,
                      bologna_LST_SON_mean=stat_bologna_list[[8]]$mean,
                      bologna_LST_SON_sd=stat_bologna_list[[9]]$mean,
                      bologna_LST_slope_sd=stat_bologna_list[[10]]$mean
)
write.csv(bologna_df,"bologna_df_night.xls",row.names = FALSE)


palermo_df=data.frame(palermo_sum_cs=stat_palermo_list$palermo_csuolo_night_stats$sum,
                      palermo_mean_cs=stat_palermo_list$palermo_csuolo_night_stats$mean,
                      palermo_std_cs=stat_palermo_list$palermo_csuolo_night_stats$std,
                      palermo_LST_DJF_mean=stat_palermo_list[[2]]$mean,
                      palermo_LST_DJF_sd=stat_palermo_list[[3]]$mean,
                      palermo_LST_JLA_mean=stat_palermo_list[[4]]$mean,
                      palermo_LST_JLA_sd=stat_palermo_list[[5]]$mean,
                      palermo_LST_MAM_mean=stat_palermo_list[[6]]$mean,
                      palermo_LST_MAM_sd=stat_palermo_list[[7]]$mean,
                      palermo_LST_SON_mean=stat_palermo_list[[8]]$mean,
                      palermo_LST_SON_sd=stat_palermo_list[[9]]$mean,
                      palermo_LST_slope_sd=stat_palermo_list[[10]]$mean
)
write.csv(palermo_df,"palermo_df_night.xls",row.names = FALSE)


nacs=which(!is.na(stat_roma_list$roma_csuolo_stats$sum))
numero_pixel_roma=max(stat_roma_list$roma_csuolo_stats$count)
nalst=which(stat_roma_list[[2]]$max>10)

plot(stat_roma_list$roma_csuolo_stats$sum[nalst],stat_roma_list[[2]]$max[nalst])
summary(lmRob(stat_roma_list[[2]]$max[nalst]~stat_roma_list$roma_csuolo_stats$mean[nalst]))

# require rasterstats python
max(stat_roma_list$roma_csuolo_stats$count)
max(stat_firenze_list$firenze_csuolo_stats$count)
max(stat_bologna_list$bologna_csuolo_stats$count)
max(stat_milano_list$milano_csuolo_stats$count)
max(stat_palermo_list$palermo_csuolo_stats$count)



1000000/43681


