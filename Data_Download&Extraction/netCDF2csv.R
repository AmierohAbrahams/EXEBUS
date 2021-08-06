# netCDF2csv
# The purpose of this script is to...
# The steps taken are:
# 1: Setup environment
# 2: Extract the OISST data
# 3: Extract the ERA5 wind data
# 4: Extract the ERA5  SLP data

# 1: Setup environment -------------------------------------------------------------------------------------
library(ncdf4) # library for processing netCDFs
library(ncdf4)
library(data.table)
library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)
library(doMC); doMC::registerDoMC(cores = 4)
library(heatwaveR)
library(plyr)
# bbox were determined using the parameters in the following paper
# Reduced Nearshore Warming Associated With Eastern Boundary Upwelling Systems 
# Rui Seabra 1 , Rubén Varela 2 , António M. Santos 1,3 , Moncho Gómez-Gesteira 2 ,
# Claudia Meneghesso 1,3 , David S. Wethey 4 and Fernando P. Lima 1 *
  
# Under Pressure: Climate Change, Upwelling, and Eastern Boundary  Upwelling Ecosystems
# Marisol García-Reyes 1 *, William J. Sydeman 1 , David S. Schoeman 2 ,
# Ryan R. Rykaczewski 3 , Bryan A. Black 4 , Albertus J. Smit 5 and Steven J. Bograd 6

# 2: Extract the OISST data ---------------------------------------------------------------------------------
  bbox <- data.frame(BC = c(-35, -15, 10, 20), # Benguela Current
                   CC = c(15, 45, 340, 350), # Canary Current
                   CalC = c(25, 45, 230, 250), # California Current
                   HC = c(-45.5, -7.5, 280, 290), # Humboldt Current
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))


nc.dir <- "/home/amieroh/Documents/Data/Datasets/AVHRR/OISST"
csv.dir <- "/media/amieroh/Amieroh/spatial"


read_nc <- function(ncFile, location = location, csv.dir = csv.dir) {
  coords <- bbox[,location]
  nc <- nc_open(ncFile)
  pathLen <- nchar(nc.dir) + 1 
  name.stem <-
    substr(ncFile, pathLen + 1, pathLen + 13)
  date.stamp <- substr(ncFile, pathLen + 15, pathLen + 22)
  LatIdx <- which(nc$dim$lat$vals > coords[1] & nc$dim$lat$vals < coords[2])
  LonIdx <- which(nc$dim$lon$vals > coords[3] & nc$dim$lon$vals < coords[4])
  sst <- ncvar_get(nc,
                   varid = "sst",
                   start = c(LonIdx[1], LatIdx[1], 1, 1),
                   count = c(length(LonIdx), length(LatIdx), 1, 1)) %>%
    round(4) 
  dimnames(sst) <- list(lon = nc$dim$lon$vals[LonIdx],
                        lat = nc$dim$lat$vals[LatIdx])
  nc_close(nc)
  sst <-
    as.data.table(melt(sst, value.name = "temp"), row.names = NULL) %>%
    mutate(t = ymd(date.stamp)) %>%
    na.omit()
  fwrite(sst,
         file = paste(csv.dir, "/", location, "-", name.stem, ".", strtDate, "-", endDate, ".csv", sep = ""),
         append = TRUE, col.names = FALSE)
  rm(sst)
}

# the list of files
ncList <- list.files(path = nc.dir, pattern = "*.nc", full.names = TRUE, include.dirs = TRUE)
strtDate <- str_sub(ncList[1], start = 15, end = 22)
endDate <- str_sub(ncList[length(ncList)], start = 15, end = 22)
dplyr::
## apply the function
system.time(llply(ncList, read_nc, location = "BC", csv.dir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, location = "CC", csv.dir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, location = "CalC", csv.dir = csv.dir, .parallel = TRUE))
system.time(llply(ncList, read_nc, location = "HC", csv.dir = csv.dir, .parallel = TRUE))

BC <- BC-avhrr-only-v2.Document-Document.csv
CC <- CC-avhrr-only-v2.Document-Document.csv
CalC <- CalC-avhrr-only-v2.Document-Document.csv
HC <- HC-avhrr-only-v2.Document-Document.csv


# 3: Extract the ERA5 wind data ---------------------------------------------------------------------------------
# https://cds.climate.copernicus.eu/api-how-to#use-the-cds-api-client-for-data-access
# Extracting ERA 5 wind u and v variables
# N/W/S/E coordinate format when downloading
# -26.00/15.00/-36.00/22.00
# This extraction is for the Benguela Current the same is repeated for the Canary, California and Humboldt current
# New downloading ERA5 u and v variables
# Extracting the data
# To obtain ERA 5 wind data one needs to download u and v wind variables
# These variables are collected every hour for every day
# Here I extract u and v wind variables and then convert them to daily data in order to match the OISST daily temperature data
# ERA 5 does not have daily wind variables only hourly

library(ncdf4)
library(ncdump)
library(tidync) # tidync is unfortunately not loading on the computer with more RAM for some reason
library(tidyverse)
library(reshape2)
library(lubridate)
library(stringr)
library(circular)
library(doParallel); doParallel::registerDoParallel(cores = 8) 
# RWS: I find doParallel to be more stable than doMC after recent updates to the tidyverse
# The AJ way ---------------------------------------------------------------------------------------------------------------------------------
nc <- nc_open("/home/amieroh/Downloads/data1.nc")
u10 <- ncvar_get(nc, varid = "u10")
v10 <- ncvar_get(nc, varid = "v10")
dimnames(v10) <- list(lon = nc$dim$longitude$vals,
                      lat = nc$dim$latitude$vals,
                      time = nc$dim$time$vals)

t_origin <- ncatt_get(nc, "time", "units")$value
t_origin

as.ymd <- function(x) {
  as.Date(as.POSIXct(x * 3600, origin = "1900-01-01 00:00:00.0"),
          "GMT",
          "%Y-%m-%d")
}

u10_df <- as_tibble(reshape2::melt(u10, value.name = "u10"), row.names = NULL) %>%
  mutate(time = as.ymd(time)) %>%
  na.omit()

v10_df <- as_tibble(reshape2::melt(v10, value.name = "v10"), row.names = NULL) %>%
  mutate(time = as.ymd(time)) %>%
  na.omit()

v10_df <- v10_df %>% 
  dplyr::select(-lat,-lon,-time)
wind1982 <- cbind(u10_df,v10_df)4

###################################################################################################
# Load with tidync --------------------------------------------------------
# Downloading ERA5 data from the requesting site and subset it.
# Extract and bind with current data

rm(BC_match);gc()
rm(BC_match1982);gc()
rm(BC_wind_fin);gc()

rm(wind_daily);gc()

#rm(wind_daily_2)

wind_daily <- tidync("/home/amieroh/Downloads/HC_sub/data1.nc") %>%
  hyper_tibble() %>% 
  dplyr::select(longitude, latitude, time, v10, u10) %>% 
  rename(lon = longitude,
         lat = latitude) %>% 
  mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1900-01-01"))) %>% 
  dplyr::select(-time)
combined <- rbind(combined,wind_daily)
# BC_wind <- combined
# wind_daily_2 <- tidync("/home/amieroh/Downloads/CalC_sub/data.nc") %>%
#   hyper_tibble() %>%
#   dplyr::select(longitude, latitude, time, v10, u10) %>%
#   rename(lon = longitude,
#          lat = latitude) %>%
#   mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1900-01-01"))) %>%
#  dplyr::select(-time)

#combined <- rbind(wind_daily,wind_daily_2)

# 4: Extract the ERA5 SLP data ---------------------------------------------------------------------------------

# Viewing netCDF properties
# ncin <- nc_open("/home/amieroh/Downloads/CC_test.nc")
# print(ncin)

SLP_HC <- tidync("/home/amieroh/Downloads/HC_SLP.nc") %>%
  hyper_tibble() %>% 
  dplyr::select(lon, lat, time, slp) %>% 
  mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1800-01-01 00:00:00.0"))) %>% 
  dplyr::select(-time)

SLP_BC <- tidync("/home/amieroh/Downloads/BC_SLP.nc") %>%
  hyper_tibble() %>% 
  dplyr::select(lon, lat, time, slp) %>% 
  mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1800-01-01 00:00:00.0"))) %>% 
  dplyr::select(-time)

SLP_CalC <- tidync("/home/amieroh/Downloads/CalC_SLP.nc") %>%
  hyper_tibble() %>% 
  dplyr::select(lon, lat, time, slp) %>% 
  mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1800-01-01 00:00:00.0"))) %>% 
  dplyr::select(-time)


SLP_CC <- tidync("/home/amieroh/Downloads/CC_SLP.nc") %>%
  hyper_tibble() %>% 
  dplyr::select(lon, lat, time, slp) %>% 
  mutate(t = as.Date(as.POSIXct(time * 60 * 60, origin = "1800-01-01 00:00:00.0"))) %>% 
  dplyr::select(-time)

save(SLP_BC , file = "data/SLP_BC.RData")
save(SLP_CC , file = "data/SLP_CC.RData")
save(SLP_CalC , file = "data/SLP_CalC.RData")
save(SLP_HC , file = "data/SLP_HC.RData")
