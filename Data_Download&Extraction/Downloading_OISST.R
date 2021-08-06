# NB: The packages only need to be installed from GitHub once
# install.packages("dplyr")
# install.packages("lubridate")
# install.packages("ggplot2")
# install.packages("tidync")
# install.packages("doParallel")
# install.packages("plyr")

# Load the packages once they have been downloaded and installed
# The packages we will use
library(dplyr) # A staple for most modern data management in R
library(lubridate) # Useful functions for dealing with dates
library(ggplot2) # The preferred library for data visualisation
library(tidync) # For easily dealing with NetCDF data
library(doParallel) # For parallel processing

# First we tell R where the data are on the interwebs
OISST_base_url <- "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/"

OISST_dates <- data.frame(t = seq(as.Date("1981-09-01"), as.Date("1981-10-02"), by = "day"))

# To finish up this step we add some text to those dates so they match the OISST file names
OISST_files <- OISST_dates %>%
  mutate(t_day = gsub("-", "", t),
         t_month = substr(t_day, 1, 6),
         t_year = year(t),
         file_name = paste0(OISST_base_url, t_month, "/", "oisst-avhrr-v02r01.", t_day ,".nc"))

OISST_url_daily_dl <- function(target_URL){
  dir.create("~/home/amieroh/Documents/EBUS/data/OISST", showWarnings = F)
  file_name <- paste0("~/home/amieroh/Documents/EBUS/data/OISST/",sapply(strsplit(target_URL, split = "/"), "[[", 10))
  if(!file.exists(file_name)) download.file(url = target_URL, method = "libcurl", destfile = file_name)
}

doParallel::registerDoParallel(cores = 4)

# And with that we are clear for take off
system.time(plyr::l_ply(OISST_files$file_name, .fun = OISST_url_daily_dl, .parallel = T)) # ~15 seconds


OISST_load <- function(file_name, lon1, lon2, lat1, lat2){
  OISST_dat <- tidync(file_name) %>%
    hyper_filter(lon = between(lon, lon1, lon2),
                 lat = between(lat, lat1, lat2)) %>%
    hyper_tibble() %>%
    select(lon, lat, time, sst) %>%
    dplyr::rename(t = time, temp = sst) %>%
    mutate(t = as.Date(t, origin = "1978-01-01"))
  return(OISST_dat)
}

# Locate the files that will be loaded
OISST_files <- dir("/home/amieroh/Documents/EBUS/data/OISST", full.names = T)

rm(OISST_dat);gc()

# Load the data in parallel
OISST_BC <- plyr::ldply(.data = OISST_files, .fun = OISST_load, .parallel = T,
                        lon1 = 20, lon2 = 25, lat1 = -27.5, lat2 = -35.5)

save(OISST_BC , file = "data_complete/OISST_BC.RData")

OISST_HC <- plyr::ldply(.data = OISST_files, .fun = OISST_load, .parallel = T,
                        lon1 = 280 , lon2 = 290 , lat1 = -45.5, lat2 = -7.5)
save(OISST_HC , file = "data_complete/OISST_HC.RData")

OISST_CC <- plyr::ldply(.data = OISST_files, .fun = OISST_load, .parallel = T,
                        lon1 = 340 , lon2 = 350 , lat1 = 15, lat2 = 45)
save(OISST_CC , file = "data_complete/OISST_CC.RData")

OISST_CalC <- plyr::ldply(.data = OISST_files, .fun = OISST_load, .parallel = T,
                          lon1 = 230, lon2 = 250, lat1 = 25, lat2 = 45)
save(OISST_CalC , file = "data_complete/OISST_CalC.RData")

# the netCDF file extracts OISST data from 1982 -2018 see below loading the CalC
CalC_avhrr_only_v2_Document_Document <- read_csv("~/Documents/CalC-avhrr-only-v2.Document-Document.csv",col_names = c("lon","lat","date","temp"))

##Then combining the two datsets with rbind
CalC <- rbind(CalC,OISST_CalC)
#save(CalC , file = "data/CalC.RData")

# For the next 2020 updated download rbind the new data with the data created on line 82, which is now the full dataset
OISST_global <- tidync("~/Documents/EBUS/data/OISST/oisst-avhrr-v02r01.20180701.nc") %>%
  # hyper_filter(lon = between(lon, lon1, lon2),
  #              lat = between(lat, lat1, lat2)) %>%
  hyper_tibble() %>%
  dplyr::select(lon, lat, time, sst) %>%
  dplyr::rename(t = time, temp = sst) %>%
  mutate(t = as.Date(t, origin = "1978-01-01"))

save(OISST_global, file = "data/OISST_global.RData")

