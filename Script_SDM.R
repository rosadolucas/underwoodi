###  Script by Lucas Rosado Mendonça for master's analysis
###  March 2021 to August 2023
###  Mendonça et al. In Prep. Effects of the urban environment and climate change on a neotropical parthenogenetic lizard.


# Packages ----------------------------------------------------------------
{
  if(!require(devtools)) install.packages('devtools', dependencies = T)
  if(!require(Mapinguari)) install_github("gabrielhoc/Mapinguari")
  if(!require(magrittr)) install.packages('magrittr', dependencies = T)
  if(!require(tidyverse)) install.packages('tidyverse', dependencies = T)
  if(!require(ggplot2)) install.packages('ggplot2', dependencies = T)
  if(!require(ggmap)) install.packages('ggmap', dependencies = T)
  if(!require(raster)) install.packages('raster', dependencies = T)
  if(!require(geodata)) install.packages('geodata', dependencies = T)
  if(!require(mgcv)) install.packages('mgcv', dependencies = T)
  if(!require(progress)) install.packages('progress', dependencies = T)
  if(!require(rgdal)) install.packages('rgdal', dependencies = T)
  if(!require(sqldf)) install.packages('sqldf', dependencies = T)
  if(!require(maps)) install.packages('maps', dependencies = T)
  if(!require(testthat)) install.packages('testthat', dependencies = T)
  if(!require(roxygen2)) install.packages('roxygen2', dependencies = T)
  if(!require(sf)) install.packages('sf', dependencies = T)
} 

rm(list=ls())

# Dir's -------------------------------------------------------------------
var_dir <- 'D:/Variaveis_ENM'
main_dir <- 'D:/Lucas/Adequabilidade_G_underwoodi'

# Distribution records ----------------------------------------------------
underwoodiD <- read.csv(paste0(main_dir, "/Distribution_G_underwoodi.csv"), 
                        header = T, sep = ",")
head(underwoodiD)

underwoodiD <- underwoodiD %>% 
  mutate_if(is.character, as.factor) 

# Ploting data 
underwoodi_bbox <- ggmap::make_bbox(lat = lat,
                             lon = lon,
                             data =  underwoodiD,
                             f = 0.2)

underwoodi_big <- ggmap::get_map(location = underwoodi_bbox,
                          source = "stamen",
                          maptype = "toner")

ggmap(underwoodi_big)+
  geom_point(data = underwoodiD,
             mapping = aes(x = lon, y = lat),
             size = 2,
             colour = "red")

# Filter records by altitude. Removing points in sea
alt <- raster::brick(paste0(var_dir,'/Variables_crop/Current/elevation/elevation.tif')) # importing alt data.

underwoodiD_clean <- Mapinguari::clean_points(coord = underwoodiD[,-1],
                                  merge_dist = 10000, # eliminate points within 10km from each other
                                  filter_layer = !alt == 0) # !is.na(alt) keeps points only where alt is not NA. 

ggmap(underwoodi_big)+
  geom_point(data = underwoodiD_clean,
             mapping = aes(x = lon, y = lat),
             size = 2,
             colour = "red")

underwoodiD_clean2 <- underwoodiD_clean %>% 
  drop_na()

coordinates(underwoodiD_clean2) = underwoodiD_clean2[,c("lon","lat")]

# Saving clean records
saveRDS(underwoodiD_clean2, paste0(main_dir, '/underwoodiD_clean2.rds'))
write.table(underwoodiD_clean2, paste0(main_dir, '/UnderwoodiD_clean2.txt'), sep = '\t', row.names = F)

# Importing clean records
underwoodiD_clean2 <- read.table(paste0(main_dir, '/UnderwoodiD_clean2.txt'), header = T)
coordinates(underwoodiD_clean2) = underwoodiD_clean2[,c("lon","lat")]

# Creating biological rasters with Mapinguari ----------------------------------------------
### We did it following the vignette available on https://gabrielhoc.github.io/Mapinguari.html. 

setwd(var_dir)

# Fitting physiological models -------------------------#
TPC <- read.csv(paste0(main_dir, "/TPC.csv"), header = T, sep = ",")

head(TPC)
summary(TPC)

TPC_2 <- TPC %>% 
  mutate_if(is.character, as.factor) %>% 
  drop_na()

summary(TPC_2)

# Best model 
tpc_gamm <- mgcv::gamm(speed_mean ~ s(temp, by = distribution, bs = 'cs')
                       + distribution
                       + FAL_L
                       + FTL_L
                       + SVL,
                       random = list(id = ~ 1),
                       data = TPC_2)

pred_tpc <- Mapinguari::get_predict(tpc_gamm, type = 'response')

pred_tpc(distribution = 'Native', temp = 30, SVL = 40, FAL_L = 0.22862, FTL_L = 0.042468) # testing if its working

# Spatializing a physiological model ------------------#

### Present ------------------#

# Current variables 
underwoodi_climate_present <-
  Mapinguari::get_rasters(
    var = c('tmin', 'tmax', 'bio'), # variables names
    scenario = "present", 
    raster_path = paste0(var_dir, '/Variables_crop/Current/GRD'),
    margin = 5)

saveRDS(underwoodi_climate_present, 'underwoodi_climate_present.rds')

# importing .rds files saved
underwoodi_climate_present <- readRDS(paste0(main_dir,'./underwoodi_climate_present.rds'))

plot(underwoodi_climate_present$present$tmax_01)

# Spliting variables
present_tmax <- raster::subset(underwoodi_climate_present$present, grep('tmax', names(underwoodi_climate_present$present), value = T))

present_tmin <- raster::subset(underwoodi_climate_present$present, grep('tmin', names(underwoodi_climate_present$present), value = T))

plot(present_tmax)

# Creating model by month
underwoodi_present_performance_tmax <- 
  Mapinguari::transform_rasters(raster_stack = present_tmax, 
                                perf.RR.tmax.min = pred_tpc(temp = tmax, distribution = 'Native', 
                                                            SVL = min(TPC_2$SVL),
                                                            FAL_L = min(TPC_2$FAL_L),
                                                            FTL_L = min(TPC_2$FTL_L)),
                                perf.RR.tmax.max = pred_tpc(temp = tmax, distribution = 'Native', 
                                                            SVL = max(TPC_2$SVL),
                                                            FAL_L = max(TPC_2$FAL_L),
                                                            FTL_L = max(TPC_2$FTL_L)),
                                perf.RR.tmax.mean = pred_tpc(temp = tmax, distribution = 'Native', 
                                                             SVL = mean(TPC_2$SVL),
                                                             FAL_L = mean(TPC_2$FAL_L),
                                                             FTL_L = mean(TPC_2$FTL_L)),
                                perf.AM.tmax.min = pred_tpc(temp = tmax, distribution = 'Neonative',
                                                            SVL = min(TPC_2$SVL),
                                                            FAL_L = min(TPC_2$FAL_L),
                                                            FTL_L = min(TPC_2$FTL_L)),
                                perf.AM.tmax.max = pred_tpc(temp = tmax, distribution = 'Neonative',
                                                            SVL = max(TPC_2$SVL),
                                                            FAL_L = max(TPC_2$FAL_L),
                                                            FTL_L = max(TPC_2$FTL_L)),
                                perf.AM.tmax.mean = pred_tpc(temp = tmax, distribution = 'Neonative',
                                                             SVL = mean(TPC_2$SVL),
                                                             FAL_L = mean(TPC_2$FAL_L),
                                                             FTL_L = mean(TPC_2$FTL_L)),
                                ncores = 7) 

underwoodi_present_performance_tmin <- 
  Mapinguari::transform_rasters(raster_stack = present_tmin, 
                                perf.RR.tmin.min = pred_tpc(temp = tmin, distribution = 'Native',
                                                            SVL = min(TPC_2$SVL),
                                                            FAL_L = min(TPC_2$FAL_L),
                                                            FTL_L = min(TPC_2$FTL_L)),
                                perf.RR.tmin.max = pred_tpc(temp = tmin, distribution = 'Native',
                                                            SVL = max(TPC_2$SVL),
                                                            FAL_L = max(TPC_2$FAL_L),
                                                            FTL_L = max(TPC_2$FTL_L)),
                                perf.RR.tmin.mean = pred_tpc(temp = tmin, distribution = 'Native',
                                                             SVL = mean(TPC_2$SVL),
                                                             FAL_L = mean(TPC_2$FAL_L),
                                                             FTL_L = mean(TPC_2$FTL_L)),
                                perf.AM.tmin.min = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                            SVL = min(TPC_2$SVL),
                                                            FAL_L = min(TPC_2$FAL_L),
                                                            FTL_L = min(TPC_2$FTL_L)),
                                perf.AM.tmin.max = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                            SVL = max(TPC_2$SVL),
                                                            FAL_L = max(TPC_2$FAL_L),
                                                            FTL_L = max(TPC_2$FTL_L)),
                                perf.AM.tmin.mean = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                             SVL = mean(TPC_2$SVL),
                                                             FAL_L = mean(TPC_2$FAL_L),
                                                             FTL_L = mean(TPC_2$FTL_L)),
                                ncores = 7)

saveRDS(underwoodi_present_performance_tmin, paste0(main_dir, '/underwoodi_present_performance_limbs_tmin.rds'))
saveRDS(underwoodi_present_performance_tmax, paste0(main_dir, '/underwoodi_present_performance_limbs_tmax.rds'))

raster::writeRaster(underwoodi_present_performance_tmin, paste0(var_dir, "/Underwoodi_raster/Underwoodi_performance_present_tmin.grd")) # if wanna see

# importing .rds files 
underwoodi_present_performance_tmin <- readRDS(paste0(main_dir, './underwoodi_present_performance_tmin.rds'))
underwoodi_present_performance_tmax <- readRDS(paste0(main_dir, './underwoodi_present_performance_tmax.rds'))

# Summaries performance raster 
underwoodi_present_performance_tmax_summaries <- 
  transform_rasters(raster_stack = underwoodi_present_performance_tmax,
                    perf.RR.average.tmax.min = mean(perf.RR.tmax.min),
                    perf.RR.average.tmax.max = mean(perf.RR.tmax.max),
                    perf.RR.average.tmax.mean = mean(perf.RR.tmax.mean),
                    perf.AM.average.tmax.min = mean(perf.AM.tmax.min),
                    perf.AM.average.tmax.max = mean(perf.AM.tmax.max),
                    perf.AM.average.tmax.mean = mean(perf.AM.tmax.mean),
                    ncores = 7)

underwoodi_present_performance_tmin_summaries <- 
  transform_rasters(raster_stack = underwoodi_present_performance_tmin,
                    perf.RR.average.tmin.min = mean(perf.RR.tmin.min),
                    perf.RR.average.tmin.max = mean(perf.RR.tmin.max),
                    perf.RR.average.tmin.mean = mean(perf.RR.tmin.mean),
                    perf.AM.average.tmin.min = mean(perf.AM.tmin.min),
                    perf.AM.average.tmin.max = mean(perf.AM.tmin.max),
                    perf.AM.average.tmin.mean = mean(perf.AM.tmin.mean),
                    ncores = 7)

saveRDS(underwoodi_present_performance_tmax_summaries, 
        paste0(main_dir, '/underwoodi_present_performance_limbs_tmax_summaries.rds'))

saveRDS(underwoodi_present_performance_tmax_summaries, 
        paste0(var_dir, '/Underwoodi_raster/underwoodi_present_performance_limbs_tmax_summaries.rds'))

raster::writeRaster(underwoodi_present_performance_tmin_summaries, paste0(var_dir, "/Underwoodi_raster/present_performance_limbs_tmin_summaries.grd")) 

# importing .rds files 
summaries_present_performance_tmax <- 
  readRDS(paste0(main_dir, './underwoodi_present_performance_limbs_tmax_summaries.rds'))
summaries_present_performance_tmin <- 
  readRDS(paste0(main_dir, './underwoodi_present_performance_limbs_tmin_summaries.rds'))

# Future -----------------------------------#
# Importing future variables. tmin and tmax need to be import one at a time

setwd(paste0(var_dir, '/Variables_crop/Future/GRD'))

list <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
list.2 <- list[grep(patter = 'tmin', list)]; list.2
# list.2 <- list.2[!list.2 %in% c("./tmin_ipslcm6alr_370_100")] # if need to exclude some variable

i <- 1
for (i in 1:length(list.2)) {
  setwd(paste0(list.2[i])) # go into the GCM's folder
  
  list.3 <- list.files(path = ".")
  list.4 <- list.3[grep(patter = '.grd', list.3)]
  
  var <- raster::stack(list.4); var
  
  ordem_var <- c("tmin_01", "tmin_02", "tmin_03", "tmin_04", "tmin_05",
                 "tmin_06", "tmin_07", "tmin_08", "tmin_09", "tmin_10",
                 "tmin_11", "tmin_12" )
  
  var_names <- names(var)
  
  ifelse(var_names == ordem_var, print('Its OK'), stop('ERROR'))
  
  print('----------------------------------------------------------------')
  
  obj_name <- substring(list.2[i], 3, nchar(list.2[i]))
  
  assign(paste0(obj_name), var)
  
  print('_________________________________________________________________________')
  print(paste0('############ Done: ', i, '-', length(list.2), '  #################'))
  print('_________________________________________________________________________')
  
  setwd('./..')
} # importing future variables

# Spatializing future physiological model
var_future_tmax <- c(tmax_ipslcm6alr_370_100, tmax_ipslcm6alr_370_60, 
                     tmax_ipslcm6alr_585_100, tmax_ipslcm6alr_585_60,
                     tmax_miroc6_370_100, tmax_miroc6_370_60, 
                     tmax_miroc6_585_100, tmax_miroc6_585_60,
                     tmax_mriesm20_370_100, tmax_mriesm20_370_60, 
                     tmax_mriesm20_585_100, tmax_mriesm20_585_60,
                     tmax_ukesm10ll_370_100, tmax_ukesm10ll_370_60, 
                     tmax_ukesm10ll_585_100, tmax_ukesm10ll_585_60)

var_future_tmin <- c(tmin_ipslcm6alr_370_100, tmin_ipslcm6alr_370_60, 
                     tmin_ipslcm6alr_585_100, tmin_ipslcm6alr_585_60,
                     tmin_miroc6_370_100, tmin_miroc6_370_60, 
                     tmin_miroc6_585_100, tmin_miroc6_585_60,
                     tmin_mriesm20_370_100, tmin_mriesm20_370_60, 
                     tmin_mriesm20_585_100, tmin_mriesm20_585_60,
                     tmin_ukesm10ll_370_100, tmin_ukesm10ll_370_60, 
                     tmin_ukesm10ll_585_100, tmin_ukesm10ll_585_60)

summary(tmax_ukesm10ll_370_100) == summary(var_future_tmax[[13]]) # assessing if the c() is in the same order of var_future

# progress bar to spatializing future ecophysiological raster. This will take some time, go to read a book. 

n_iter <- length(list.2)

pb <-
  progress_bar$new(format='(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimate time remaining :  :eta]',
                   total = n_iter,
                   complete = '=',
                   incomplete = '-',
                   current = '>',
                   clear = F,
                   width = 100)

# Model spatializing. Run with Tmin and Tmax
j <- 1
for (j in 1:n_iter) {
  
  pb$tick()
  
  performance_future <- 
    Mapinguari::transform_rasters(raster_stack = var_future_tmin[[j]],
                                  perf.fut.RR.tmin.min = pred_tpc(temp = tmin, distribution = 'Native', 
                                                                  SVL = min(TPC_2$SVL),
                                                                  FAL_L = min(TPC_2$FAL_L),
                                                                  FTL_L = min(TPC_2$FTL_L)),
                                  perf.fut.RR.tmin.max = pred_tpc(temp = tmin, distribution = 'Native', 
                                                                  SVL = max(TPC_2$SVL),
                                                                  FAL_L = max(TPC_2$FAL_L),
                                                                  FTL_L = max(TPC_2$FTL_L)),
                                  perf.fut.RR.tmin.mean = pred_tpc(temp = tmin, distribution = 'Native', 
                                                                   SVL = mean(TPC_2$SVL),
                                                                   FAL_L = mean(TPC_2$FAL_L),
                                                                   FTL_L = mean(TPC_2$FTL_L)),
                                  perf.fut.AM.tmin.min = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                                  SVL = min(TPC_2$SVL),
                                                                  FAL_L = min(TPC_2$FAL_L),
                                                                  FTL_L = min(TPC_2$FTL_L)),
                                  perf.fut.AM.tmin.max = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                                  SVL = max(TPC_2$SVL),
                                                                  FAL_L = max(TPC_2$FAL_L),
                                                                  FTL_L = max(TPC_2$FTL_L)),
                                  perf.fut.AM.tmin.mean = pred_tpc(temp = tmin, distribution = 'Neonative',
                                                                   SVL = mean(TPC_2$SVL),
                                                                   FAL_L = mean(TPC_2$FAL_L),
                                                                   FTL_L = mean(TPC_2$FTL_L)),
                                  ncores = 7)
  
  assign(paste0('future_performance_', substring(list.2[j], 3, nchar(list.2[j]))), 
         performance_future)
  
  saveRDS(performance_future, paste0(main_dir,'/future_performance_limbs_', substring(list.2[j], 3, nchar(list.2[j])), '.rds'))
  
  print(paste0('Done: ', j, ' - ', length(list.2)))
  
}

# importing .rds files. Run with Tmin and Tmax
setwd(main_dir)

list <- list.files(path = ".", full.names = TRUE, recursive = FALSE)
list.2 <- list[grep(patter = 'future_performance_limbs_tmax', list)]; list.2
list.3 <- list.2 # [!list.2 %in% c("./summaries_future_performance_limbs_tmax_ipslcm6alr_370_100.rds")]; list.3 # if need to exclude some variable from list

for (r in 1:length(list.3)) {
  
  performance_raster <- readRDS(paste0(substring(list.2[r], 3, nchar(list.2[r]))))
  assign(paste0(substring(list.2[r], 3, nchar(list.2[r])-4)), performance_raster)
  print(paste0('Done: ', r, '-', length(list.3)))  
  
}

# Summaring future rasters
performance_var_future <- c(future_performance_limbs_tmax_ipslcm6alr_370_100, future_performance_limbs_tmax_ipslcm6alr_370_60,
                            future_performance_limbs_tmax_ipslcm6alr_585_100, future_performance_limbs_tmax_ipslcm6alr_585_60,
                            future_performance_limbs_tmax_miroc6_370_100, future_performance_limbs_tmax_miroc6_370_60,
                            future_performance_limbs_tmax_miroc6_585_100, future_performance_limbs_tmax_miroc6_585_60,
                            future_performance_limbs_tmax_mriesm20_370_100, future_performance_limbs_tmax_mriesm20_370_60,
                            future_performance_limbs_tmax_mriesm20_585_100, future_performance_limbs_tmax_mriesm20_585_60,
                            future_performance_limbs_tmax_ukesm10ll_370_100, future_performance_limbs_tmax_ukesm10ll_370_60,
                            future_performance_limbs_tmax_ukesm10ll_585_100, future_performance_limbs_tmax_ukesm10ll_585_60)

summary(future_performance_tmax_ukesm10ll_585_100) == summary(performance_var_future[[15]])

# progress bar
n_iter <- length(performance_var_future)
pb <-
  progress_bar$new(format='(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimate time remaining :  :eta]',
                   total = n_iter,
                   complete = '=',
                   incomplete = '-',
                   current = '>',
                   clear = F,
                   width = 100)

t <- 1
for (t in 1:length(performance_var_future)) {
  
  pb$tick()
  
  summaries_future <- 
    Mapinguari::transform_rasters(raster_stack = performance_var_future[[t]],
                                  perf.fut.RR.average.tmax.min = mean(perf.fut.RR.tmax.min),
                                  perf.fut.RR.average.tmax.max = mean(perf.fut.RR.tmax.max),
                                  perf.fut.RR.average.tmax.mean = mean(perf.fut.RR.tmax.mean),
                                  perf.fut.AM.average.tmax.min = mean(perf.fut.AM.tmax.min),
                                  perf.fut.AM.average.tmax.max = mean(perf.fut.AM.tmax.max),
                                  perf.fut.AM.average.tmax.mean = mean(perf.fut.AM.tmax.mean),
                                  ncores = 7)
  
  #assign(paste0(substring(list.3[t], 3, nchar(list.3[t])-4), '_summaries'), summaries_future)
  
  writeRaster(summaries_future, 
              paste0(var_dir, '/Underwoodi_raster/summaries_', substring(list.3[t], 3, nchar(list.3[t])-4), '.grd'))
  
  saveRDS(summaries_future, paste0(main_dir, '/summaries_', substring(list.3[t], 3, nchar(list.3[t])-4), '.rds'))
  
  print(paste0('Done: ', t , '|',length(list.3)))
  
} # summaring rasters

# If the summarise above do not work try this
for (k in 1:16) {
  
  var_dir <- 'D:/Variaveis_ENM'
  main_dir <- 'D:/Lucas/Trabalho - Desktop/Trabalho/Mestrado 2021/UFAM/Projeto/Dados/R_script/Adequabilidade_G_underwoodi'
  
  setwd(main_dir)
  
  list <- list.files(path = ".", full.names = TRUE, recursive = FALSE)
  list.2 <- list[grep(patter = 'future_performance_limbs_tmin', list)]; list.2
  list.3 <- list.2
  
  r <- k
  t <- k
  
  performance_var_future <- readRDS(paste0(substring(list.2[r], 3, nchar(list.2[r]))))
  
  gc()
  
  summaries_future <- 
    Mapinguari::transform_rasters(raster_stack = performance_var_future,
                                  perf.fut.RR.average.tmin.min = mean(perf.fut.RR.tmin.min),
                                  perf.fut.RR.average.tmin.max = mean(perf.fut.RR.tmin.max),
                                  perf.fut.RR.average.tmin.mean = mean(perf.fut.RR.tmin.mean),
                                  perf.fut.AM.average.tmin.min = mean(perf.fut.AM.tmin.min),
                                  perf.fut.AM.average.tmin.max = mean(perf.fut.AM.tmin.max),
                                  perf.fut.AM.average.tmin.mean = mean(perf.fut.AM.tmin.mean),
                                  ncores = 7)
  
  #assign(paste0(substring(list.3[t], 3, nchar(list.3[t])-4), '_summaries'), summaries_future)
  
  writeRaster(summaries_future, 
              paste0(var_dir, '/Underwoodi_raster/summaries_', substring(list.3[t], 3, nchar(list.3[t])-4), '.grd'))
  
  saveRDS(summaries_future, paste0(main_dir, '/summaries_', substring(list.3[t], 3, nchar(list.3[t])-4), '.rds'))
  
  print(paste0('Done: ', t , '|', '16'))
  
} 

# importing summaries .rds files 
setwd(main_dir)

list <- list.files(path = ".", full.names = TRUE, recursive = FALSE)
list.2 <- list[grep(patter = 'summaries_future_performance_limbs', list)]; list.2
list.3 <- list.2 #[!list.2 %in% c("./underwoodi_present_performance_tmax_summaries.rds",
#              "./underwoodi_present_performance_tmin_summaries.rds")]; list.3 

s <- 1
for (s in 1:length(list.3)) {
  
  performance_raster_future <- readRDS(paste0(substring(list.3[s], 3, nchar(list.3[s]))))
  assign(paste0(substring(list.3[s], 3, nchar(list.3[s])-4)), performance_raster_future)
  
  print(paste0('Done: ', s, '-', length(list.3)))  
}

plot(summaries_future_performance_limbs_tmin_ukesm10ll_585_100)

### Time of Activity -----------------------------------------------------------#
# Tpref gradient
underwoodiGradient <- read.csv(paste0(main_dir, "/Tpref_SDM.csv"), header = T, sep = ",", fileEncoding="latin1")
head(underwoodiGradient)

# remove outliers 
outliers_gr <- boxplot(underwoodiGradient$temp)$out
underwoodiGradient_no <- underwoodiGradient[-which(underwoodiGradient$temp %in% outliers_gr),]
boxplot(underwoodiGradient_no$temp)$out

# calcule VTmax (95% percentil) and VTmin (5% percentil)
vtmin_AM <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'AM'], 0.05)
vtmax_AM <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'AM'], 0.95)
vtmin_RR <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'RR'], 0.05)
vtmax_RR <- quantile(underwoodiGradient_no$temp[underwoodiGradient_no$site == 'RR'], 0.95)

# Microclimate method ---------------------------------------------------------------------------#
micro_dir <- 'D:/Variaveis_ENM/Microclim'

head(underwoodiD_clean2) # see processing records data

underwoodi_microclim <- 
  Mapinguari::multi_extract(raster_path = paste0(micro_dir),
                            coord = underwoodiD_clean2[-1],
                            layers = 6:18)

underwoodi_microclim <- underwoodi_microclim %>% 
  drop_na()

saveRDS(underwoodi_microclim, paste0(main_dir, '/underwoodi_microclim.rds'))

underwoodi_microclim <- readRDS(paste0(main_dir, '/underwoodi_microclim.rds'))

head(underwoodi_microclim)

hvtFUN_AM <- function(x) ifelse(x > vtmin_AM & x < vtmax_AM, 1, 0)
hvtFUN_RR <- function(x) ifelse(x > vtmin_RR & x < vtmax_RR, 1, 0)

mc_hour_AM <- 
  underwoodi_microclim %>% 
  dplyr::mutate(hvt_h.AM = hvtFUN_AM(value))

mc_hour_RR <- 
  underwoodi_microclim %>% 
  dplyr::mutate(hvt_h.RR = hvtFUN_RR(value))


mc_mh_AM <-
  mc_hour_AM %>% 
  dplyr::group_by(Lon, Lat, file_ind, layer) %>% 
  dplyr::summarise(hvt_tr_mh.AM = max(hvt_h.AM),
                   t_air_max_mh.AM = max(value))

mc_mh_RR <-
  mc_hour_RR %>% 
  dplyr::group_by(Lon, Lat, file_ind, layer) %>% 
  dplyr::summarise(hvt_tr_mh.RR = max(hvt_h.RR),
                   t_air_max_mh.RR = max(value))

mc_day_AM <- 
  mc_mh_AM %>% 
  dplyr::group_by(Lon, Lat, file_ind) %>% 
  dplyr::summarise(hvt_tr.AM = sum(hvt_tr_mh.AM),
                   t_air_tp.AM = mean(t_air_max_mh.AM))

mc_day_RR <- 
  mc_mh_RR %>% 
  dplyr::group_by(Lon, Lat, file_ind) %>% 
  dplyr::summarise(hvt_tr.RR = sum(hvt_tr_mh.RR),
                   t_air_tp.RR = mean(t_air_max_mh.RR))


hvt_tr_logistic.AM <- nls(hvt_tr.AM + 0.00001 ~ SSlogis(t_air_tp.AM, Asym, xmid, scal), data = mc_day_AM)

hvt_tr_logistic.RR <- nls(hvt_tr.RR + 0.00001 ~ SSlogis(t_air_tp.RR, Asym, xmid, scal), data = mc_day_RR)

pred_ha.AM <- 
  Mapinguari::get_predict(list(hvt_trFUN.AM = hvt_tr_logistic.AM))

pred_ha.RR <- 
  Mapinguari::get_predict(list(hvt_trFUN.RR = hvt_tr_logistic.RR))

ha_present_MicroClim_AM <- 
  Mapinguari::transform_rasters(tmax_tmin_present,
                                ha.vt.tr.AM = pred_ha.AM$hvt_trFUN.AM(t_air_tp.AM = tmax),
                                ncores = 5)

ha_present_MicroClim_RR <- 
  Mapinguari::transform_rasters(tmax_tmin_present,
                                ha.vt.tr.RR = pred_ha.RR$hvt_trFUN.RR(t_air_tp.RR = tmax),
                                ncores = 5)


plot(ha_present_MicroClim_AM$ha.vt.tr.AM_01)
plot(underwoodiD_clean2, pch=16, col = "black", cex = 1, add=TRUE)

saveRDS(ha_present_MicroClim_AM, paste0(main_dir, '/ha_present_microclim_AM.rds'))
saveRDS(ha_present_MicroClim_RR, paste0(main_dir, '/ha_present_microclim_RR.rds'))

# Summarise 
ha_present_microclim_AM_summaries <- 
  Mapinguari::transform_rasters(raster_stack = ha_present_MicroClim_AM,
                                ha_vt_tr_AM_average = mean(ha.vt.tr.AM),
                                ha_vt_tr_AM_sum = sum(ha.vt.tr.AM),
                                ncores = 5)

ha_present_microclim_RR_summaries <- 
  Mapinguari::transform_rasters(raster_stack = ha_present_MicroClim_RR,
                                ha_vt_tr_RR_average = mean(ha.vt.tr.RR),
                                ha_vt_tr_RR_sum = sum(ha.vt.tr.RR),
                                ncores = 5)

plot(ha_present_microclim_RR_summaries)

# saving .rds
saveRDS(ha_present_microclim_AM_summaries, paste0(main_dir, '/ha_present_microclim_AM_summaries.rds'))
saveRDS(ha_present_microclim_RR_summaries, paste0(main_dir, '/ha_present_microclim_RR_summaries.rds'))

ha_present_wc <- raster::stack(ha_present_MicroClim_AM, ha_present_MicroClim_RR)
ha_present_wc_summaries <- raster::stack(ha_present_microclim_AM_summaries, ha_present_microclim_RR_summaries) 

saveRDS(ha_present_wc, paste0(main_dir, '/ha_present_wc.rds'))
saveRDS(ha_present_wc_summaries, paste0(main_dir, '/ha_present_wc_summaries.rds'))

writeRaster(ha_present_wc_summaries, paste0(var_dir, '/Underwoodi_raster/summaries_ha_present_wc.grd'))

# importing .rds
ha_present_wc <- readRDS(paste0(main_dir, '/ha_present_wc.rds'))
ha_present_wc_summaries <- readRDS(paste0(main_dir, '/ha_present_wc_summaries.rds'))

plot(ha_present_wc_summaries)

# Future --------------------------------------#
var_future_mc <- c(tmax_ipslcm6alr_370_100, tmax_ipslcm6alr_370_60, 
                   tmax_ipslcm6alr_585_100, tmax_ipslcm6alr_585_60,
                   tmax_miroc6_370_100, tmax_miroc6_370_60, 
                   tmax_miroc6_585_100, tmax_miroc6_585_60,
                   tmax_mriesm20_370_100, tmax_mriesm20_370_60, 
                   tmax_mriesm20_585_100, tmax_mriesm20_585_60,
                   tmax_ukesm10ll_370_100, tmax_ukesm10ll_370_60, 
                   tmax_ukesm10ll_585_100, tmax_ukesm10ll_585_60) # go back to line 233 if you start here and import future variables

names(var_future_mc) <- c('tmax_ipslcm6alr_370_100', 'tmax_ipslcm6alr_370_60', 
                          'tmax_ipslcm6alr_585_100', 'tmax_ipslcm6alr_585_60',
                          'tmax_miroc6_370_100', 'tmax_miroc6_370_60', 
                          'tmax_miroc6_585_100', 'tmax_miroc6_585_60',
                          'tmax_mriesm20_370_100', 'tmax_mriesm20_370_60', 
                          'tmax_mriesm20_585_100', 'tmax_mriesm20_585_60',
                          'tmax_ukesm10ll_370_100', 'tmax_ukesm10ll_370_60', 
                          'tmax_ukesm10ll_585_100', 'tmax_ukesm10ll_585_60')

# Progress bar 
n_iter <- length(var_future_mc)
pb <-
  progress_bar$new(format ='(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimate time remaining :  :eta]',
                   total = n_iter,
                   complete = '=',
                   incomplete = '-',
                   current = '>',
                   clear = F,
                   width = 100)

h <- 1
for (h in 1:n_iter) {
  pb$tick()
  
  out <- 
    Mapinguari::transform_rasters(var_future_mc[[h]],
                                  ha.vt.tr.AM = pred_ha.AM$hvt_trFUN.AM(t_air_tp.AM = tmax),
                                  ha.vt.tr.RR = pred_ha.RR$hvt_trFUN.RR(t_air_tp.RR = tmax),
                                  ncores = 7)
  
  saveRDS(out, paste0(main_dir, '/ha_mc_', 
                      substring(names(var_future_mc)[h], 6, nchar(names(var_future_mc)[h])), '.rds'))  
  
  print(paste0('Done: ', h, '|', n_iter))
  
} # Spatializing ha_mc

# importing .rds
list_5 <- list.files(main_dir, pattern = '.rds'); list_5
list_6 <- list_5[grep(patter = 'ha_mc', list_5)]; list_6

y <- 1
for (y in 1:length(list_6)) {
  out <- readRDS(paste0(main_dir, '/', list_6[y]))
  
  assign(substring(list_6[y], 1, nchar(list_6[y])-4), out)
  
  print(paste0('Done: ', list_6[y], ' -- ', y, ':', length(list_6)))
} 

ha_mc_futures <- c(ha_mc_ipslcm6alr_370_100, ha_mc_ipslcm6alr_370_60,
                   ha_mc_ipslcm6alr_585_100, ha_mc_ipslcm6alr_585_60,
                   ha_mc_miroc6_370_100, ha_mc_miroc6_370_60,
                   ha_mc_miroc6_585_100, ha_mc_miroc6_585_60,
                   ha_mc_mriesm20_370_100, ha_mc_mriesm20_370_60,
                   ha_mc_mriesm20_585_100, ha_mc_mriesm20_585_60,
                   ha_mc_ukesm10ll_370_100, ha_mc_ukesm10ll_370_60,
                   ha_mc_ukesm10ll_585_100, ha_mc_ukesm10ll_585_60)

names(ha_mc_futures) <- c('ha_mc_ipslcm6alr_370_100', 'ha_mc_ipslcm6alr_370_60',
                          'ha_mc_ipslcm6alr_585_100', 'ha_mc_ipslcm6alr_585_60',
                          'ha_mc_miroc6_370_100', 'ha_mc_miroc6_370_60',
                          'ha_mc_miroc6_585_100', 'ha_mc_miroc6_585_60',
                          'ha_mc_mriesm20_370_100', 'ha_mc_mriesm20_370_60',
                          'ha_mc_mriesm20_585_100', 'ha_mc_mriesm20_585_60',
                          'ha_mc_ukesm10ll_370_100', 'ha_mc_ukesm10ll_370_60',
                          'ha_mc_ukesm10ll_585_100', 'ha_mc_ukesm10ll_585_60')

# Summaring future rasters
x <- 1
for (x in 1:length(ha_mc_futures)) {
  summaries_future_mc <- 
    Mapinguari::transform_rasters(raster_stack = ha_mc_futures[[x]],
                                  ha_vt_tr_mc_RR_average = mean(ha.vt.tr.RR),
                                  ha_vt_tr_mc_RR_sum = sum(ha.vt.tr.RR),
                                  ha_vt_tr_mc_AM_average = mean(ha.vt.tr.AM),
                                  ha_vt_tr_mc_AM_sum = sum(ha.vt.tr.AM),
                                  ncores = 4)
  
  raster::writeRaster(summaries_future_mc, paste0(var_dir, '/Underwoodi_raster/summaries_ha_mc_', 
                                                  substring(names(ha_mc_futures)[x], 7, nchar(names(var_future_mc)[x])), '.grd'))
  
  saveRDS(summaries_future_mc, paste0(main_dir, '/summaries_ha_mc_', 
                                      substring(names(ha_mc_futures)[x], 7, nchar(names(var_future_mc)[x])), '.rds')) 
  
  
  print(paste0('Done: ', names(ha_mc_futures)[x], ' ==>> ', x, '|', length(ha_mc_futures)))
  
} 

# importing .rds
list_7 <- list.files(main_dir, pattern = '.rds'); list_7
list_8 <- list_7[grep(patter = 'summaries', list_7)]; list_8

y <- 1
for (y in 1:length(list_8)) {
  out <- readRDS(paste0(main_dir, '/', list_8[y]))
  
  assign(substring(list_8[y], 1, nchar(list_8[y])-4), out)
  
  print(paste0('Done: ', list_8[y], ' -- ', y, ':', length(list_8)))
} 

plot(raster::stack(paste0(var_dir, '/Underwoodi_raster/summaries_ha_mc_miroc6_585_10.grd')))


# First Variable selection ---------------------------------------------------
# Run this for BC variables firt, and after do with ecophysiological raster

# VIF
set.seed(22111962)
var <- raster::stack(list.files('D:/Variaveis_ENM/Variables_ENMTML/Present', 
                                full.names = T, pattern = '.tif')) 

# creating 1000 random points from mask
mask <- var$layer.5
rnd.points <- dismo::randomPoints(mask, 10000)

# VIF of bc
vif_var <- usdm::vifstep(var, th = 10)

# Subset environmental stack
env.selected <- usdm::exclude(var, vif_var)

plot(env.selected)  


# Distribution records bias - EnvSample ---------------------------------------------------------------
envSample <- function(coord, filters, res, do.plot=TRUE){
  
  n<- length (filters)
  pot_points<- list ()
  for (i in 1:n){
    k<- filters [[i]] [!is.na(filters[[i]])]
    ext1<- range (k)
    ext1 [1]<- ext1[1]- 1
    x<- seq(ext1[1],ext1[2], by=res[[i]])
    pot_points[[i]]<- x
  }
  pot_p<- expand.grid(pot_points)
  
  ends<- NULL
  for (i in 1:n){
    fin<- pot_p [,i] + res[[i]]
    ends<- cbind (ends, fin)
  }
  
  pot_pp<- data.frame (pot_p, ends)
  pot_pp<- data.frame (pot_pp, groupID=c(1:nrow (pot_pp)))
  rows<- length (filters[[1]])
  filter<- data.frame(matrix(unlist(filters), nrow=rows))
  real_p<- data.frame (coord, filter)
  
  names_real<- c("lon", "lat")
  names_pot_st<- NULL
  names_pot_end<- NULL
  sql1<- NULL
  for (i in 1:n){
    names_real<- c(names_real, paste ("filter", i, sep=""))
    names_pot_st<- c(names_pot_st, paste ("start_f", i, sep=""))
    names_pot_end<- c(names_pot_end, paste ("end_f", i, sep=""))
    sql1<- paste (sql1, paste ("real_p.filter", i, sep=""), sep=", ")   
  }
  
  names (real_p)<- names_real
  names (pot_pp)<- c(names_pot_st, names_pot_end, "groupID")
  
  conditions<- paste ("(real_p.filter", 1, "<= pot_pp.end_f", 1,") and (real_p.filter", 1, "> pot_pp.start_f", 1, ")", sep="")
  for (i in 2:n){
    conditions<- paste (conditions, 
                        paste ("(real_p.filter", i, "<= pot_pp.end_f", i,") and (real_p.filter", i, "> pot_pp.start_f", i, ")", sep=""), 
                        sep="and")
  }
  
  selection_NA<- sqldf(paste ("select real_p.lon, real_p.lat, pot_pp.groupID",   
                              sql1, "from pot_pp left join real_p on", conditions, sep=" "))
  
  selection<- selection_NA [complete.cases(selection_NA),]
  
  final_points<- selection[!duplicated(selection$groupID), ]
  coord_filter<- data.frame (final_points$lon, final_points$lat) 
  names (coord_filter)<- c("lon", "lat")
  
  if (do.plot==TRUE){
    par (mfrow=c(1,2), mar=c(4,4,0,0.5))
    plot (filters[[1]], filters[[2]], pch=19, 
          col="grey50", xlab="Filter 1", ylab="Filter 2")
    points (final_points$filter1, final_points$filter2, 
            pch=19, col="#88000090")
    plot (coord, pch=19, col="grey50")
    map(add=T)
    points (coord_filter, pch=19, col="#88000090")
    
  }
  coord_filter
}

coord <- read.csv(paste0(main_dir, "/Distribution_G_underwoodi.csv"), 
                        header = T, sep = ",")

setwd ("D:/Variaveis_ENM/Variables_ENMTML/Present") # where are the variables selected to model (BC + ecophysiological + landcover + slope)
var <- list.files(pattern=".tif")
wc <- stack(var)
data <- extract(wc, coord)
data <- as.data.frame(data)

sample <- envSample(coord, filters = list(data$layer.2, data$layer.4,
                                          data$layer.5, data$layer.6,
                                          data$layer.7, data$layer.10,
                                          data$layer.11), 
                    res=list(10, 10, 10, 10, 10, 10, 10), do.plot=TRUE)

write.table(sample, './EnvSample_cleaned.txt', sep = ",")

# Creating the SDM with ENMTML - Andrade et al 2020 ---------------------------------------------
if(!require(ENMTML)) devtools::install_github("andrefaa/ENMTML", force=T, dependence = T)  
if(!require(viridis)) install.packages('viridis', dependencies = T)
if(!require(tmap)) install.packages('tmap', dependencies = T)

ENMTML(pred_dir = paste0(var_dir, '/Variables_ENMTML/Present'),
       occ_file = paste0(main_dir, '/EnvSample_cleaned.txt'),
       result_dir = paste0(main_dir, '/Script_ENM/Result'),
       proj_dir = paste0(var_dir, '/Variables_ENMTML/Future'),
       sp = "sp",
       x = "lon",
       y = "lat",
       min_occ = 20,
       thin_occ = c(method = 'USER-DEFINED', distance = '10'),
       colin_var = NULL,
       imp_var = TRUE,
       sp_accessible_area = c(method = 'MASK', 
                              filepath = paste0(main_dir, '/minimum_convex_poly/polygon_500km.shp')),
       pres_abs_ratio = 1,
       pseudoabs_method = c(method='GEO_CONST', width='50'),
       part = c(method = 'BLOCK'),
       save_part = F, save_final = T,
       algorithm = c('MXD', "RDF", "GAU", "SVM"),
       thr = c(type ='SORENSEN'),
       msdm= NULL,
       ensemble = c(method = 'W_MEAN', metric = 'Sorensen'),
       extrapolation = F,
       cores = 6)

plot(raster::brick(paste0(main_dir, '/Script_ENM/Result/Ensemble/W_MEAN/G_underwoodi.tif')), 
     col = viridis(10))

### Predictions ensemble -----------#
spp370_60 <- mean(raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ipslcm6alr_370_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/miroc6_370_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/mriesm20_370_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ukesm10ll_370_60/Ensemble/W_MEAN/G_underwoodi.tif')))

spp370_100 <- mean(raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ipslcm6alr_370_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/miroc6_370_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/mriesm20_370_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ukesm10ll_370_100/Ensemble/W_MEAN/G_underwoodi.tif')))

spp585_60 <- mean(raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ipslcm6alr_585_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/miroc6_585_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/mriesm20_585_60/Ensemble/W_MEAN/G_underwoodi.tif')),
                  raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ukesm10ll_585_60/Ensemble/W_MEAN/G_underwoodi.tif')))

spp585_100 <- mean(raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ipslcm6alr_585_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/miroc6_585_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/mriesm20_585_100/Ensemble/W_MEAN/G_underwoodi.tif')),
                   raster::brick(paste0(main_dir, '/Script_ENM/Result/Projection/ukesm10ll_585_100/Ensemble/W_MEAN/G_underwoodi.tif')))

par(mfrow = c(2,2))
plot(spp370_60, col = viridis(10))
plot(spp370_100, col = viridis(10))
plot(spp585_60, col = viridis(10))
plot(spp585_100, col = viridis(10))

writeRaster(spp585_60, paste0(main_dir, '/Script_ENM/Result/spp585_60.tif'))

# Binary Plot -------------------------------------------------------------
# by Josué Azevedo and Isabela Oliveira <------------------------------------##

# Objects that you will need 
# Projection from present
log_pred <- raster::brick(paste0(main_dir, '/Script_ENM/Result/Ensemble/W_MEAN/G_underwoodi.tif')) 

# Projection from future
SPP370_2060_pred <- raster::brick(paste0(main_dir, '/Script_ENM/Result/spp370_60.tif'))
SPP585_2060_pred <- raster::brick(paste0(main_dir, '/Script_ENM/Result/spp585_60.tif'))
SPP370_2100_pred <- raster::brick(paste0(main_dir, '/Script_ENM/Result/spp370_100.tif'))
SPP585_2100_pred <- raster::brick(paste0(main_dir, '/Script_ENM/Result/spp585_100.tif'))

# Threshold value
ev2 <- 0.32

# Distribution records
occ_raw <- read.table(paste0(main_dir, '/Distribution_G_underwoodi.txt'), header = T)
occ_raw <- drop_na(occ_raw)

occ_xy <- occ_raw[,-1]
coordinates(occ_raw) = occ_raw[,c("lon","lat")]
occs.xy.z <- read.table(paste0(main_dir, '/Distribution_G_underwoodi.txt'), header = T)

# Species 
genero <- 'Gymnophthalmus'
especie <- 'underwoodi'

# Neotropico shapefile
Neotropico <- readOGR(dsn = 'D:/qgis/Shapes/World_Countries_(Generalized)', 
                      layer = 'World_Countries__Generalized_')

fortify(Neotropico) # if you need to see the col names

buffer500km <- readOGR(paste0(main_dir, '/minimum_convex_poly/polygon_500km.shp'))
limite <- c(-90, -45, -10, 25)

# Current binary plot
pa.pa <- raster(log_pred) # creating an empty raster
pa.pa[] <- ifelse(log_pred[] >= ev2, 1, NA)
plot(pa.pa, col = c('white', 'blue'), legend = FALSE)
points(occ_xy, pch = 16)
maps::map('world', add = TRUE)
title(paste0("Thresholded_map_current"))  


# Detect continous clamps in the rasters
clumped <- clump(pa.pa, directions = 8)
sp_buf <- buffer500km # here we use a buffer remove suitable areas that are not continuous with the species records. Check if it makes sense for your species
inter <- raster::extract(clumped, sp_buf, na.rm = TRUE) 
inters <- na.exclude(inter[[1]])
raster1 <- match(clumped, inters)

current_bin <- pa.pa # if wanna use full map binary

# Rasters whithout extrapolations
current_bin <- raster(raster1) # creating an empty raster
current_bin[] <- ifelse(raster1[] >= 1, 1, NA)
current_bin <- raster::mask(current_bin, buffer500km)
plot(current_bin, col = c('white', 'blue'), legend = FALSE)
points(occ_raw, pch = 16)
maps::map('world', add = TRUE)
title(paste0("Thresholded_map_current"))  

# Save binary raster
outputname = paste(paste(genero,"_", especie, "_threshold" , '.tif', sep=""))
writeRaster(current_bin, outputname, bylayer = TRUE, options = c("COMPRESS=DEFLATE"), format="GTiff", overwrite=TRUE)

### Future binary plot
binary_fut = function(future_pred){
  pr.pa_fut <- raster(future_pred) # creating an empty raster
  pr.pa_fut[] <- ifelse(future_pred[] >= ev2, 1, NA)
  clumped <- clump(pr.pa_fut, directions = 8)
  sp_buf <- buffer500km ## here we use a buffer remove suitable areas that are not continuous with the species records
  inter <- raster::extract(clumped, sp_buf, na.rm = TRUE) 
  inters <- na.exclude(inter[[1]])
  raster1 <- match(clumped, inters)
  #raster binário sem extrapolação
  current_bin_fut <- raster(raster1) # creating an empty raster
  current_bin_fut[] <- ifelse(raster1[] >= 1, 1, NA)
  current_bin_fut <- raster::mask(current_bin_fut, buffer500km)
  plot(current_bin_fut)
  points(occ_raw, pch = 16)
  maps::map('world',add = TRUE)
  title(paste0("Thresholded_map"))  
  return(pr.pa_fut)
}

#Layers futuro binario
SPP370_2060_bin = binary_fut(future_pred = SPP370_2060_pred)
SPP585_2060_bin = binary_fut(future_pred = SPP585_2060_pred)
SPP370_2100_bin = binary_fut(future_pred = SPP370_2100_pred)
SPP585_2100_bin = binary_fut(future_pred = SPP585_2100_pred)

plot(current_bin,col = "blue", main = "Current", legend = FALSE)
plot(SPP585_2060_bin, col = "red", main = "2060", add = TRUE, legend = FALSE)
plot(SPP585_2100_bin, col = "yellow", main = "2100", add = TRUE, legend = FALSE)


# Calculate area
current_area = sum(na.exclude(getValues(raster::area(current_bin, na.rm = TRUE))))
SPP370_2060_area = sum(na.exclude(getValues(raster::area(SPP370_2060_bin, na.rm = TRUE))))
SPP585_2060_area = sum(na.exclude(getValues(raster::area(SPP585_2060_bin, na.rm = TRUE))))
SPP370_2100_area = sum(na.exclude(getValues(raster::area(SPP370_2100_bin, na.rm = TRUE))))
SPP585_2100_area = sum(na.exclude(getValues(raster::area(SPP585_2100_bin, na.rm = TRUE))))

write.csv(cbind.data.frame(especie, current_area, 
                           SPP370_2060_area, 
                           SPP585_2060_area, 
                           SPP370_2100_area, 
                           SPP585_2100_area),
          paste0(main_dir, "/Script_ENM/Result/Areas_distribruicao_binario.csv"), 
          append = TRUE, col.names = TRUE, row.names = FALSE)

# Loss, gain, stable and expansion 
raster_data_1 = current_bin
raster_data_2 = SPP370_2060_bin
raster_data_3 = SPP370_2100_bin
raster_data_4 = raster::mask(raster_data_3,raster_data_1, maskvalue = 1)
period = ''

#raster_data_1 <- raster::mask(raster_data_1, buffer500km)
#raster_data_2 <- raster::mask(raster_data_2, buffer500km)
#raster_data_3 <- raster::mask(raster_data_3, buffer500km)
#raster_data_4 <- raster::mask(raster_data_4, buffer500km)

# Raster to points for ggplot
TAB_RDA1 <- as.data.frame(rasterToPoints(raster_data_1))
colnames(TAB_RDA1)[3] <- "value"  
TAB_RDA1$variable <- factor(rep("Current", nrow(TAB_RDA1)), 
                            levels = c("Current"))
head(TAB_RDA1)
tail(TAB_RDA1)        

# Reverse values to fit the legend
TAB_RDA1$value = 1 - TAB_RDA1$value

# Raster to points for ggplot
TAB_RDA2 <- as.data.frame(rasterToPoints(raster_data_2))
colnames(TAB_RDA2)[3] <- "value"  
TAB_RDA2$variable <- factor(rep("Current", nrow(TAB_RDA2)), 
                            levels = c("Current"))
head(TAB_RDA2)
tail(TAB_RDA2)        
# Reverse values to fit the legend
TAB_RDA2$value = 1 - TAB_RDA2$value

# Raster to points for ggplot
TAB_RDA3 <- as.data.frame(rasterToPoints(raster_data_3))
colnames(TAB_RDA3)[3] <- "value"   
TAB_RDA3$variable <- factor(rep("Current", nrow(TAB_RDA3)), 
                            levels = c("Current"))
head(TAB_RDA3)
tail(TAB_RDA3)        
# Reverse values to fit the legend
TAB_RDA3$value = 1 - TAB_RDA3$value

# Raster to points for ggplot
TAB_RDA4 <- as.data.frame(rasterToPoints(raster_data_4))
colnames(TAB_RDA4)[3] <- "value" 
TAB_RDA4$variable <- factor(rep("Current", nrow(TAB_RDA4)), 
                            levels = c("Current"))
head(TAB_RDA4)
tail(TAB_RDA4)        
# Reverse values to fit the legend
TAB_RDA4$value = 1 - TAB_RDA4$value

# Plot 
my_ggplot_2 = ggplot() +
  geom_polygon(data = Neotropico, 
               aes(long, lat, group = group), fill = NA, col = NA, size = .2)+
  annotate(geom = "raster", x = TAB_RDA1$x, y = TAB_RDA1$y, fill = "#d73027")+ # loss
  annotate(geom = "raster", x = TAB_RDA2$x, y = TAB_RDA2$y, fill = "#fee090")+ # gain in 2060 loss in 2100
  annotate(geom = "raster", x = TAB_RDA3$x, y = TAB_RDA3$y, fill = "#4575b4")+ # stable
  annotate(geom = "raster", x = TAB_RDA4$x, y = TAB_RDA4$y, fill = "#984ea3")+ # gain until 2100
  geom_polygon(data = Neotropico, 
               aes(long, lat, group = group), 
               fill = NA, col = "black", size = .55)+
  geom_point(data = occs.xy.z, 
             aes(x = lon, y = lat), col = "black", fill = "black", size = 2.5, shape = 1, stroke = 0.75)+
  coord_cartesian(xlim = extent_plot[1:2], ylim = extent_plot[3:4])+
  ggsn::scalebar(x.min = extent_plot[1], x.max = extent_plot[2], 
                 y.min = extent_plot[3], y.max = extent_plot[4],
                 dist = 200, dist_unit = "km", model = 'WGS84', st.size = 2.5,
                 location = "bottomleft", transform=T, border.size=.5)+
  ggsn::north(anchor = c("x" = -45, "y" = 25), x.min = extent_plot[1], x.max = extent_plot[2], 
              y.min = extent_plot[3], y.max = extent_plot[4])+
  xlab(NULL) + ylab(NULL) +
  ggtitle(bquote(italic(.(genero)~.(especie))~.(period)))+
  guides(fill = guide_legend(title = "Suitability")) +
  theme_bw(base_size = 14, base_family = "Times") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), 
        panel.background = element_blank(), strip.text = element_text(size=14))

my_ggplot_2

pdf(paste(paste(genero,"_", especie, "Range_Shift_ssp585" , '.pdf', sep="")))
print(my_ggplot_2)
dev.off()

# converting to raster to save as .tif to edit in Qgis
{
  coordinates(TAB_RDA1) <- ~ x + y
  gridded(TAB_RDA1) <- T
  rasterTABRDA1 <- raster(TAB_RDA1)
  
  writeRaster(rasterTABRDA1, paste0(main_dir, '/Script_ENM/Result/rasterTABRDA1.tif'))
}

# Variables importance --------------------------------------------
variableImportance <- read.table(paste0(main_dir, '/Script_ENM/Result/Algorithm/RDF/Response Curves & Variable Importance/VariableImportance.txt'), header = T)
head(variableImportance)

layer_var <- c('layer.1', 'layer.2', 'layer.3', 'layer.4', 'layer.5',
               'layer.6', 'layer.7', 'layer.8', 'layer.9', 'layer.10', 'layer.11')  

names_layer <- c('Isothermality (bio3)', 'Temperature Seasonality (bio4)',
                 'Max Temperature of Warmest Month (bio5)', 'Precipitation of Wettest Month (bio13)',
                 'Precipitation Seasonality (bio15)', 'Clay content (0.0 m depth)', 'Landcover',
                 'Slope', 'havt_tr_AM_sum', 'perf_AM_tmin_min', 'perf_RR_tmax_mean')

layer_datframe <- cbind.data.frame(layer_var, names_layer)

variableImportance$names_var <- plyr::mapvalues(variableImportance$Variables, 
                                                to = layer_datframe$names_layer, 
                                                from = layer_datframe$layer_var)
variableImportance <- 
  data.frame(importance = variableImportance$IncNodePurity)

ggplot(variableImportance,
       mapping = aes(x = reorder(names_layer, importance), y = importance, fill = importance)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  theme(line = element_blank()) +
  scale_fill_gradient(low = "light blue", high = "blue", guide = "none") +
  scale_x_discrete("Variables") + scale_y_continuous("Increase in Node Impurity")

# Converting .grd to .tif - ENMTML do not work with grd -------------------------------------------------
setwd(paste0(var_dir, '/Underwoodi_raster'))
list_11 <- list.files(path = paste0(var_dir, '/Underwoodi_raster'), pattern = '.grd'); list_11
list_12 <- list_11[grep(patter = 'summaries_future_performance_limbs_tmax', list_11)]; list_12

r <- 1
for (r in 1:length(list_12)) {
  out <- raster::stack(list_12[r])
  names(out)
  
  #names(out) <- c("summary_ha_vt_RR_average", "summary_ha_vt_tr_RR_average", 
  #                "summary_ha_vt_RR_sum", "summary_ha_vt_tr_RR_sum",    
  #                "summary_ha_vt_AM_average", "summary_ha_vt_tr_AM_average", 
  #                "summary_ha_vt_AM_sum", "summary_ha_vt_tr_AM_sum") 
  
  names(out) <- c("summary_perf_RR_average_tmax_min", "summary_perf_RR_average_tmax_max", 
                  "summary_perf_RR_average_tmax_mean", "summary_perf_AM_average_tmax_min",
                  "summary_perf_AM_average_tmax_max", "summary_perf_AM_average_tmax_mean")
  
  folder_name <- substring(list_12[r], 41, nchar(list_12[r])-4)
  setwd(paste0(var_dir, '/Variables_ENMTML/Future/', folder_name))
  #dir.create('limbs_var')
  setwd('./limbs_var')
  raster::writeRaster(out, filename = names(out), bylayer = T, format = 'Gtiff')
  
  setwd(paste0(var_dir, '/Underwoodi_raster'))
  
  print(paste0('Done: ', r, '|', length(list_12)))
}

plot(out$summary_perf_RR_average_tmax_mean)
plot(raster::brick(
  paste0(var_dir, '/Variables_ENMTML/Present/summary_ha_vt_AM_average.tif')))

# Creating a polygon ------------------------------------------------------
if(!require(sf)) install.packages('sf')
if(!require(tidygraph)) install.packages('tidygraph')

sf_use_s2(T)

data <- read.csv(paste0(main_dir, '/Distribution_G_underwoodi.csv'), header = T)
dat_sf <- st_as_sf(data, coords = c('lon', 'lat'), crs = 4326)

conv <- st_simplify(st_buffer(st_convex_hull(st_union(st_geometry(dat_sf))), dist = 500000), dTolerance = 2)

plot(conv)

setwd(paste0(main_dir, '/minimum_convex_poly'))
conv %>% 
  as_tibble() %>% 
  st_as_sf() %>% 
  st_write('polygon_500km.shp')



