# 2_upwelling_Identification.R.R
# The purpose of this script is to...
# The steps taken are:
# 1: Setup environment
# 2: Find the coastal pixels
# 3: Calculate upwelling and the metrics

# 1: Setup environment ----------------------------------------------------

# Loading Libraries
# library(circular)
library(gridExtra)
library(geosphere)
library(tidyverse)
library(heatwaveR)
## devtools::install_github("robwschlegel/coastR")
library(coastR)
library(FNN)
source("functions/theme.R")

# Load data
# Data loaded here was created in "1_Temp_wind_data.R"
load("data/CC_complete.RData")
load("data/CalC_complete.RData")
load("data/HC_complete.RData")
load("data/BC_complete.RData")

CC_complete <- CC_complete %>% 
  mutate(lon = lon - 360) # Fix all use of a static -360

CalC_complete <- CalC_complete %>% 
  mutate(lon = lon - 360)

HC_complete <- HC_complete %>% 
  mutate(lon = lon - 360)

# Loading the temperature data this is the OISST data extracted to the regions (See netCDF2CSV script in the data extraction folder)
load("~/Documents/EBUS/data/BC.RData")
load("~/Documents/EBUS/data/HC.RData")
load("~/Documents/EBUS/data/CC.RData")
load("~/Documents/EBUS/data/CalC.RData")

# 2: Find the coastal pixels ----------------------------------------------
# Isolate the unique pixel coordinates

coord_func <- function(df){
  coords <- df %>% 
    dplyr::select(lon, lat) %>%
    mutate(lon = lon - 360) %>% # Fix
    unique()
}

CC_coords <- coord_func(df = CC)
CalC_coords <- coord_func(df = CalC)
HC_coords <- coord_func(df = HC)

BC_coords <- BC %>% 
  dplyr::select(lon, lat) %>% 
  unique()

# Take coastal coordinates from a global map
coastline_func <- function(df){
  coastline <- fortify(maps::map(xlim = c(min(df$lon, na.rm = T), 
                                          max(df$lon, na.rm = T)), 
                                 ylim = c(min(df$lat, na.rm = T), 
                                          max(df$lat, na.rm = T)), 
                                 plot = F, interior = F, fill = T, lforce = "e", map = "world"))
}

CC_coastline <- coastline_func(df = CC_coords)
CalC_coastline <- coastline_func(df = CalC_coords)
HC_coastline <- coastline_func(df = HC_coords)
BC_coastline <- coastline_func(df = BC_coords)

plot_func <- function(df){
  ggplot(data = df, aes(x = long, y = lat)) +
    geom_polygon(colour = "black", fill = "grey70", aes(group = group))
}

plot_func(df = CC_coastline)
plot_func(df = CalC_coastline)
plot_func(df = HC_coastline)
plot_func(df = BC_coastline)

# Find which of the EBUS pixels match closest to the coastal pixels

coastal_index_func <- function(coord_df,coastline_df){
  BC_coastal_index <- as.vector(knnx.index(as.matrix(coord_df[, c("lon", "lat")]),
                                           as.matrix(coastline_df[ ,c("long", "lat")]), k = 1))
}

CC_coastal_index <- coastal_index_func(coord_df = CC_coords, coastline_df = CC_coastline)
CalC_coastal_index <- coastal_index_func(coord_df = CalC_coords, coastline_df = CalC_coastline)
HC_coastal_index <- coastal_index_func(coord_df = HC_coords, coastline_df = HC_coastline)
BC_coastal_index <- coastal_index_func(coord_df = BC_coords, coastline_df = BC_coastline)

BC_coastal_coords <- unique(BC_coords[BC_coastal_index,]) # Rather use "distinct" instead of "unique"
CC_coastal_coords <- unique(CC_coords[CC_coastal_index,])
CalC_coastal_coords <- unique(CalC_coords[CalC_coastal_index,])
HC_coastal_coords <- unique(HC_coords[HC_coastal_index,])

# Find the coastal angle for each point

BC_transects <- transects(BC_coastal_coords, spread = 2, alongshore = T) %>% 
  dplyr::rename(coastal_angle = heading) %>% 
  mutate(coastal_angle = round(coastal_angle),
         coastal_angle = coastal_angle+180,
         coastal_angle = ifelse(coastal_angle > 360, coastal_angle-360, coastal_angle))

transect_func <- function(df){
  transects <- transects(df, spread = 2) %>% 
    dplyr::rename(coastal_angle = heading) %>% 
    mutate(coastal_angle = round(coastal_angle))
}

CC_transects <- transect_func(df = CC_coastal_coords)
CalC_transects <- transect_func(df = CalC_coastal_coords)
HC_transects <- transect_func(df = HC_coastal_coords)

# Bind it all together

CC_coastal <- left_join(CC_coastal_coords, CC_complete, by = c("lon", "lat")) %>% 
  left_join(CC_transects, by = c("lon", "lat"))
# save(CC_coastal, file = "data/CC_coastal.RData")
rm(CC_complete, CC_temp); gc()

CalC_coastal <- left_join(CalC_coastal_coords, CalC_complete, by = c("lon", "lat")) %>% 
  left_join(CalC_transects, by = c("lon", "lat"))
# save(CalC_coastal, file = "data/CalC_coastal.RData")
rm(CalC_complete, CalC_temp); gc()

HC_coastal <- left_join(HC_coastal_coords, HC_complete, by = c("lon", "lat")) %>% 
  left_join(HC_transects, by = c("lon", "lat"))
# save(HC_coastal, file = "data/HC_coastal.RData")
rm(HC_complete, HC_temp); gc()

BC_coastal <- left_join(BC_coastal_coords, BC_complete, by = c("lon", "lat")) %>% 
  left_join(BC_transects, by = c("lon", "lat"))
# save(BC_coastal, file = "data/BC_coastal.RData")


# 3: Calculate upwelling and the metrics ----------------------------------

# Determining the upwelling index per coastal pixel
upwelling_func <- function(df){
  UI <- df %>%  
    mutate(ui = wind_spd * (cos(wind_dir_from - coastal_angle)), 
           ui_TF = ifelse(ui > 0, TRUE, FALSE)) 
}

CC_UI <- upwelling_func(df = CC_coastal) %>% 
  dplyr::rename(t = date)
# save(CC_UI, file = "data/CC_UI.RData")

CalC_UI <- upwelling_func(df = CalC_coastal) %>% 
  dplyr::rename(t = date)
# save(CalC_UI, file = "data/CalC_UI.RData")

HC_UI <- upwelling_func(df = HC_coastal) %>% 
  dplyr::rename(t = date)
# save(HC_UI, file = "data/HC_UI.RData")


BC_UI <- upwelling_func(df = BC_coastal) %>% 
  dplyr::rename(t = date)
# save(BC_UI, file = "data/BC_UI.RData")

# The custom function for detecting upwelling and extracting only the metrics
detect_event_custom <- function(df){
  res <- detect_event(df, threshClim2 = df$ui_TF, minDuration = 1, coldSpells = T)$event
  return(res)
}

# Calculate the upwelling event metrics
clim_func <- function(df){
  clim <- df %>% 
    group_by(lon, lat) %>% 
    nest() %>% 
    mutate(clim = purrr::map(data, ts2clm, pctile = 25, climatologyPeriod = c("1982-01-01", "2011-12-31"))) %>% 
    select(-data) %>%
    unnest(cols = clim) %>% 
    ungroup()
}

CC_clim <- clim_func(df = CC_UI)
# save(CC_clim, file = "data/CC_clim.RData")

CalC_clim <- clim_func(df = CalC_UI)
# save(CalC_clim, file = "data/CalC_clim.RData")

<<<<<<< HEAD
HC_clim <- clim_func(df = HC_UI) 
=======
HC_clim <- clim_func(df = HC_UI)
>>>>>>> 09fe0ed976130e21c8d543760db451cd72e16a1e
# save(HC_clim, file = "data/HC_clim.RData")

BC_clim <- clim_func(df = BC_UI) 
# save(BC_clim, file = "data/BC_clim.RData")

# Calculate the upwelling metrics
UI_metrics_func <- function(df,clim_df){
  UI_metrics <- df %>% 
    left_join(clim_df, by = c("lon", "lat", "t", "temp")) %>% 
    group_by(lon, lat) %>% 
    nest() %>% 
    mutate(event = purrr::map(data, detect_event_custom)) %>% 
    select(-data) %>%
    unnest(cols = event) %>% 
    ungroup()
}

CC_UI_metrics <- UI_metrics_func(df = CC_UI, clim_df = CC_clim)
CalC_UI_metrics <- UI_metrics_func(df = CalC_UI, clim_df = CalC_clim)
HC_UI_metrics <- UI_metrics_func(df = HC_UI, clim_df = HC_clim)
BC_UI_metrics <- UI_metrics_func(df = BC_UI, clim_df = BC_clim)

# save(CC_UI_metrics, file = "data/CC_UI_metrics.RData")
# save(CalC_UI_metrics, file = "data/CalC_UI_metrics.RData")
# save(HC_UI_metrics, file = "data/HC_UI_metrics.RData")
# save(BC_UI_metrics, file = "data/BC_UI_metrics.RData")
