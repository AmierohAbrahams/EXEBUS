# 4_site_map\
# The steps taken in this script are:
# 1: Setup environment
# 2: Plotting the EBUS
# 3: Using plotdap to plot EBUSlibrary(tidyverse) # Base suite of functions
library(raster)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(PBSmapping)
library(viridis)
library(plotdap)
library(heatwaveR) # For detecting MHWs
# cat(paste0("heatwaveR version = ", packageDescription("heatwaveR")$Version))
library(FNN) # For fastest nearest neighbour searches
library(tidync) # For a more tidy approach to managing NetCDF data
# library(SDMTools) # For finding points within polygons
library(lubridate)

# 4: Plotting Robs suggestion
# 5: Plotting paper discription

# 1: Setup environment ----------------------------------------------------------------------------------------------------


# 2: Plotting EBUS ----------------------------------------------------------------------------------------------------------
# Default CRS
# +proj=longlat +datum=WGS84 +no_defs
# Load the Global Self-consistent, Hierarchical, High-resolution Geography Database
# Use the full resolution version
gshhsDir <- "/home/amieroh/Documents/Data/Datasets/gshhg-bin-2.3.7"
b <- c(5, 10, 15, 20, 25,30)
colors <- c('#EFF573', '#EFF573', '#EFF573', '#3DA394', '#3A828C', '#456075')


bbox <- data.frame(BC = c(-35, -25, 15, 20), # Benguela Current
                   CC = c(25, 35, 340, 355), # Canary Current
                   CalC = c(35, 45, 225, 240), # California Current
                   HC = c(-17.5, -7.5, 275, 290), # Humboldt Current
                   row.names = c("latmin", "latmax", "lonmin", "lonmax"))

# Make a coastline for the world in sf format
coastline <- importGSHHS(paste0(gshhsDir, "/gshhs_l.b"),
                         xlim = c(0, 360), ylim = c(-90, 90), maxLevel = 1, useWest = FALSE)
polygon <- coastline %>%
  group_by(PID) %>%
  st_as_sf(coords = c("X", "Y"), crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON")

# a simple theme for the map
theme_opts <- list(theme(panel.border = element_rect(colour = "black", size = 0.4, fill = NA),
                         axis.text = element_text(size = 16),
                         axis.line = element_line(colour = "black", size = 0.2),
                         panel.background = element_rect(colour = "black", fill = NA),
                         panel.grid.minor = element_blank(),
                         # panel.grid.major = element_blank(),
                         panel.grid.major = element_line(colour = "black", linetype = "dashed", size = 0.2),
                         axis.ticks = element_line(colour = "black")))


BC <- importGSHHS(paste0(gshhsDir, "/gshhs_f.b"),
                  xlim = c(15.00, 20.00), ylim = c(-35.00, -25.00), maxLevel = 1, useWest = FALSE)
(BC_map <- ggplot() +
    geom_raster(data = BC_complete, aes(x = lon, y = lat, fill = temp)) +
    geom_polygon(data = BC, aes(x = X, y = Y, group = PID), col = "black", fill = "grey60", size = 0.2) +
    coord_fixed(ratio = 1, expand = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "", accuracy = .01),
                       breaks = c(39.00, 41.00, 43.00)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = "", accuracy = .01)) +
    labs(title = NULL, x = NULL, y = NULL) +
  scale_fill_gradientn(colors = colors,breaks = b)+
    theme_opts +
    theme(axis.text = element_text(size = 12)))

HC <- importGSHHS(paste0(gshhsDir, "/gshhs_f.b"),
                  xlim = c(275.00, 290.00), ylim = c(-17.5, -7.5), maxLevel = 1, useWest = FALSE)
(HC_map <- ggplot() +
    geom_raster(data = HC_complete, aes(x = lon, y = lat, fill = temp)) +
    geom_polygon(data = HC, aes(x = X, y = Y, group = PID), col = "black", fill = "grey60", size = 0.2) +
    coord_fixed(ratio = 1, expand = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "", accuracy = .01),
                       breaks = c(39.00, 41.00, 43.00)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = "", accuracy = .01)) +
    labs(title = NULL, x = NULL, y = NULL) +
    scale_fill_gradientn(colors = colors,breaks = b)+
    theme_opts +
    theme(axis.text = element_text(size = 12)))


CC <- importGSHHS(paste0(gshhsDir, "/gshhs_f.b"),
                  xlim = c(340.00, 355.00), ylim = c(25.00, 35.00), maxLevel = 1, useWest = TRUE)
(CC_map <- ggplot() +
    geom_polygon(data = CC, aes(x = X, y = Y, group = PID), col = "black", fill = "grey60", size = 0.2) +
    coord_fixed(ratio = 1, expand = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "", accuracy = .01)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = "", accuracy = .01),
                       breaks = seq(from = -6.40, to = -5.80, by = 0.20)) +
    labs(title = NULL, x = NULL, y = NULL) +
    scale_fill_gradientn(colors = colors,breaks = b)+
    theme_opts +
    theme(axis.text = element_text(size = 12)))



CalC <- importGSHHS(paste0(gshhsDir, "/gshhs_f.b"),
                    xlim = c(225.00, 240.00), ylim = c(35.00, 45.00), maxLevel = 1, useWest = FALSE)
(CalC_map <- ggplot() +
    geom_raster(data = CalC_complete, aes(x = lon, y = lat, fill = temp)) +
    geom_polygon(data = CalC, aes(x = X, y = Y, group = PID), col = "black", fill = "grey60", size = 0.2) +
    coord_fixed(ratio = 1, expand = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°E", sep = "", accuracy = .01)) +
    scale_y_continuous(expand = c(0, 0), labels = scales::unit_format(unit = "°S", sep = "", accuracy = .01),
                       breaks = seq(from = -6.40, to = -5.80, by = 0.20)) +
    labs(title = NULL, x = NULL, y = NULL) +
    scale_fill_gradientn(colors = colors,breaks = b)+
    theme_opts +
    theme(axis.text = element_text(size = 12)))
  
# 3: Using plotdap to plot EBUS ----------------------------------------------------------------------------------------------------------

library(plotdap)
plotdap()
plotdap("base")
sstInfo <- rerddap::info('erdVHsstaWS3day')
# get latest 3-day composite sst
viirsSST <- rerddap::griddap(sstInfo, 
                             latitude = c(41., 31.), 
                             longitude = c(-128., -115), 
                             time = c('last','last'), 
                             fields = 'sst')
xpos <- c(135.25, 240.25)
ypos <- c(20.25, 60.25)
zpos <- c(70.02, 70.02)
remove <- c("UK:Great Britain", "France", "Spain", "Algeria", "Mali", "
            Burkina Faso", "Ghana", "Togo")
#subset world2Hires with those countries removed
w <- map("world2Hires", plot = FALSE, fill = TRUE, ylim = ypos, xlim = xpos)
w <- map("world2Hires", regions = w$names[!(w$names %in% remove)], 
         plot = FALSE, fill = TRUE, ylim = ypos, xlim = xpos)
# plot result
plotdap(mapData = w)
# write plot to disk using the Cairo package
library(sf)
#> Linking to GEOS 3.7.2, GDAL 2.4.2, PROJ 5.2.0
library(mapdata)
plotdap(mapTitle = "Grid over Land") %>%
  add_griddap(
    viirsSST, 
    ~sst, 
    fill = "thermal"
  )

library(cshapes)

world <- cshp(date=as.Date("2008-1-1"))
world.points <- fortify(world, region='COWCODE')

# Matching world points with the cordinates in OISST data

OISST_dat_remove <- OISST_dat %>% 
  mutate(lon = lon - 360)

ggplot(world.points, aes(long,lat,group=group)) + 
  geom_polygon() +
  ggplot(OISST_dat, aes(x = lon, y = lat))+
  geom_raster(aes(fill = temp))  

# 4: Plotting Robs suggesion----------------------------------------------------------------------------------------------------------------------------------

load("data/OISST_global.RData") # Created in Data_extraction folder/Downloading_OISST.R
site_squares <- read_csv("data/site_squares.csv")

OISST_global <- OISST_global %>%
  mutate(lon = ifelse(lon > 180, lon - 360, lon))

map_base <- ggplot2::fortify(maps::map(fill = TRUE, col = "grey80", plot = FALSE)) %>%
  dplyr::rename(lon = long) %>%
  mutate(group = ifelse(lon > 180, group+9999, group),
         lon = ifelse(lon > 180, lon-360, lon))

load("data_official/south_BC.RData")
load("data_official/north_BC.RData")
load("data_official/chile.RData")
load("data_official/peru.RData")
load("data_official/south_CalC.RData")
load("data_official/north_CalC.RData")
load("data_official/Canary_current.RData")

south_BC <- south_BC %>% 
  dplyr::select(lat,lon) %>% 
  distinct()
north_BC <- north_BC %>% 
  dplyr::select(lat,lon) %>% 
  distinct()

chile <- chile %>% 
  dplyr::select(lat,lon) %>% 
  distinct() 

write_csv(chile, path = "chile.csv")

peru <- peru %>% 
  dplyr::select(lat,lon) %>% 
  distinct()

Canary_current <- Canary_current %>% 
  dplyr::select(lat,lon) %>% 
  distinct()

write_csv(Canary_current, path = "Canary_current.csv")

south_CalC <- south_CalC %>% 
  dplyr::select(lat,lon) %>% 
  distinct()
north_CalC <- north_CalC %>% 
  dplyr::select(lat,lon) %>% 
  distinct()
# ggplot(map_base, aes(lon, lat, group = group)) + 
#   geom_polygon()
# world <- ne_countries(scale = "medium", returnclass = "sf")
# class(world)
# ggplot(data = world) + 
#   geom_sf(color = "black", fill = "black")

# SA_map <- ggplot() +
#   geom_sf(data = world_ne, col = "black", fill = NA, size = 0.6) +
#   # geom_sf(data = SA_prov_coast, col = "black", fill = "grey80", size = 0.4) +
#   geom_sf(data = SA_prov_coast, col = "black", fill = NA, size = 0.3) +
#   coord_sf(xlim = c(xmin, xmax),
#            ylim = c(ymin, ymax),
#            expand = FALSE) +
#   scale_x_continuous(breaks = seq(-15, 35, by = 5)) +
#   scale_y_continuous(breaks = seq(-35, -25, by = 5)) +
#   labs(x = NULL, y = NULL) +
#   theme_map2()

Map <- ggplot() +
  # geom_point(data = south_BC, aes(x = lon, y = lat), colour = "Red") +
  geom_raster(data = OISST_global, aes(x = lon, y = lat, fill = temp), show.legend = TRUE) +
  geom_polygon(data = map_base, aes(x = lon, y = lat, group = group)) +
  geom_rect(data = site_squares, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              colour = "red", fill = NA, size = 0.75) +
  annotate("text", label = "D", x = 4.0, y = -25.0,
           size =3, angle = 0, colour = "black") +
  annotate("text", label = "C", x = -25.0, y = 29.5,
           size =3, angle = 0, colour = "black") +
  annotate("text", label = "B", x = -88, y = -30,
           size =3, angle = 0, colour = "black") +
  annotate("text", label = "A", x = -133, y = 40.0,
           size =3, angle = 0, colour = "black") +
  scale_x_continuous(breaks = seq(-150, 150, by = 50),
                     labels  = c("150°W", "100°W", "50°W", "0", "50°E", "100°E", "150°E")) +
  scale_y_continuous(breaks = seq(-60, 60, by = 30),
                     labels = c("60°S", "30°S", "0", "30°N", "60°N")) +
  coord_quickmap(ylim = c(-78.5, 90), expand = FALSE) +
  scale_fill_gradientn("SST (°C)", values = scales::rescale(c(-1, 7,19,26)),
                       colors = c("lightcyan1", "orchid1", "skyblue", "blue3")) +
  geom_point(data = south_BC, aes(x = lon, y = lat), colour = "Red", size = 1)+
  geom_point(data = north_BC, aes(x = lon, y = lat), colour = "Blue", size = 1)+
  geom_point(data = chile, aes(x = lon, y = lat), colour = "Red", size = 1)+
  geom_point(data = peru, aes(x = lon, y = lat), colour = "Blue", size = 1)+
  geom_point(data = south_CalC, aes(x = lon, y = lat), colour = "Red", size = 1)+
  geom_point(data = north_CalC, aes(x = lon, y = lat), colour = "Blue", size = 1)+
  geom_point(data = Canary_current, aes(x = lon, y = lat), colour = "Red", size = 1)+
  labs(x = NULL, y = NULL) +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1),
        axis.text = element_text(colour = "black", size = 13, family = "Palatino"),
        axis.title = element_text(colour = "black", size = 20, family = "Palatino"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        legend.text=element_text(size=5, family = "Palatino"),
        legend.title = element_text(size = 6, family = "Palatino"),
        legend.key.size = unit(0.5, "cm"))

ggsave(filename = "Map.jpg", plot = last_plot(), width=180, height = 200, units = "mm",dpi = 300, device = "jpg", path = "figures/")


  # 5: Plotting paper description------------------------------------------------------------------------------------------------------------------------------------

# Regions surrounding the Benguela EBUS, are dominated by anticyclonic high-pressure cells with quasi-stationary positions, resulting in abundant southerly 
# and south easterly winds (Risien et al., 2004; Hagen et al., 2009). The South Atlantic Ocean High is situated along the west, drawing cool, dry air onto the
# west of the subcontinent (Van Heerden and Hurry, 1998). Solar heating during summer may result in the development of low-pressure cells, known as heat lows, 
# which are absent during the winter (Tyson and Preston-Whyte, 2000). 

load("data/south_africa_coast.RData")
load("data/BC_wind_plot.RData")

south_africa_coast <- south_africa_coast %>% 
  dplyr::rename(lon = long)
# Continental data

test <- ggplot() + 
  geom_raster(data = BC_wind_plot, aes(x = lon, y = lat, fill = wind_spd)) +
  geom_polygon(data = south_africa_coast, aes(x = lon, y = lat,group = group),
               fill = NA, colour = "black", size = 1.0, show.legend = FALSE) +
  scale_fill_distiller(palette = "RdYlGn") +
  labs(x = NULL, y = NULL,  fill = "Wind\nSpeed") +
  xlab("") + ylab("") +
  annotate("text", label = "ATLANTIC\nOCEAN", x = 15.0, y = -32.0,
           size = 5, angle = 0, colour = "black") + # RWS: I don't know that this label is necessary
  guides(fill = guide_colourbar()) +
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(14, 18, 2),
                     labels = scales::unit_format(unit = "°E", sep = "")) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(-34, -31), labels = c("34°S", "33°S", "32°S", "31°S")) +
                     # labels = scales::unit_format(unit = "°S", sep = "")) + # RWS: If you label the latitude with degrees South, the numbers should no longer be negative
  coord_fixed(ratio = 1, xlim = c(12.5, 20), ylim = c(-35, -30),  expand = TRUE) +
  theme(axis.text = element_text(size = rel(1.2), colour = "black", family = "Palatino"),
        axis.text.x = element_text(vjust = 1, family = "Palatino"),
        axis.text.y = element_text(hjust = 1, family = "Palatino"),
        axis.title.x = element_text(vjust = 0, family = "Palatino"),
        axis.title.y = element_text(angle = 90, vjust = 0.3, family = "Palatino"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.4, "cm"),
        panel.background = element_rect(fill = "white", colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 2),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = rel(1.0), family = "Palatino"),
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text.x = element_text(),
        strip.text.y = element_text(angle = -90),        
        plot.background = element_rect(colour = "white"),
        plot.title = element_text(size = rel(1.2), vjust = 1),
        plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"),
        legend.position = "top",
        legend.title = element_text(family = "Palatino"),
        legend.text = element_text(family = "Palatino"))

