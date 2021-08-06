# 3_Analyses
# The purpose of this script is to...
# The steps taken are:
# 1: Setup environment
# 2: Observe how wind patterns change over time
# 3: ANOVA analyses
# 4: Linear models

# climate change as a result of global warning resulted in changes in wind patterns and
# ultimately lead to changes in the duration and intensity of upwelling events overtime.
# Changing upwelling region boundaries for each current

# 1: Setup environment --------------------------------------------------------------------------------------------------------------------------------------------
library(gridExtra)
library(geosphere)
library(tidyverse)
library(lubridate)
library(heatwaveR)
## devtools::install_github("robwschlegel/coastR")
library(coastR)
library(mgcv) # for gam
library(FNN)
library(broom)
library(circular)
library(grid)
source("functions/theme.R")
options(scipen=999) 

# New facet label names
supp.labs <- c("Benguela current", "Humboldt current", "California current", "Canary current")
names(supp.labs) <- c("BC","HC","CalC","CC")
my.formula <- y ~ x

# 2: Wind pattern observation --------------------------------------------------------------------------------------------------------------------------------------
# Analyses done to compare how the wind blown in a SE direction during summer months varied over a 30 year period

# The datasets used here were created in script "5_SLP.R"
load("data/CC_coastal_SLP.RData") 
load("data/CalC_coastal_SLP.RData")
load("data/BC_coastal_SLP.RData")
load("data/HC_coastal_SLP.RData")

BC <- BC_coastal_SLP%>% 
  mutate(current = "BC") 
HC <- HC_coastal_SLP %>% 
  mutate(current = "HC")
CC <- CC_coastal_SLP %>% 
  mutate(current = "CC")
CalC <- CalC_coastal_SLP %>% 
  mutate(current = "CalC")

current_winds <- rbind(BC,HC,CC,CalC)
rm(BC,BC_coastal_SLP,CalC,CalC_coastal_SLP,CC_coastal_SLP,CC, HC_coastal_SLP,HC);gc()
  # save(current_winds, file = "data/current_winds.RData")

# Then create different temporal results
# First filter out only the SE data
SE_renamed <- current_winds %>% 
  filter(wind_dir_from >= 180, wind_dir_from <= 270) %>% 
  unique()
rm(current_winds);gc()
#save(SE_renamed, file = "data/SE_renamed.RData")

# Then create diifferent temporal results
SE_summer <- SE_renamed %>% 
  filter(season == "Summer") %>% 
  group_by(current, year, season) %>% 
  summarise(count = n(),
            circ_dir = mean.circular(circular(wind_dir_from, units = "degrees")),
            circ_wspd = mean.circular(circular(wind_spd, units = "degrees")),
            mean_temp = mean(temp, na.rm = T),
            mean_SLP = mean(slp, na.rm = T))

SE_monthly <- SE_renamed %>% 
  filter(season == "Summer") %>% 
  group_by(current, year, season, month) %>% 
  summarise(count = n(),
            circ_dir = mean.circular(circular(wind_dir_from, units = "degrees")),
            circ_wspd = mean.circular(circular(wind_spd, units = "degrees")),
            mean_temp = mean(temp, na.rm = T),
            mean_SLP = mean(slp, na.rm = T))

# Wind-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Changes in the number of SE wind blown
SE_monthly_BC <- SE_monthly %>% 
  filter(current == "BC") %>% 
  mutate(no_SE = count/73)

SE_monthly_CC <- SE_monthly %>% 
  filter(current == "CC") %>% 
  mutate(no_SE = count/106)


SE_monthly_CalC <- SE_monthly %>% 
  filter(current == "CalC") %>% 
  mutate(no_SE = count/185)

SE_monthly_HC <- SE_monthly %>% 
  filter(current == "HC") %>% 
  mutate(no_SE = count/193)

SE_winds <- rbind(SE_monthly_BC,SE_monthly_CalC,SE_monthly_CC,SE_monthly_HC)

  plotA <- ggplot(data = SE_winds, aes(x = year, y = no_SE)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "SE wind events 
(count)")+
  theme_bw() +
  labs(colour = "Month") +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=14, family = "Palatino"),
    axis.title = element_text(size = 20, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 20, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 15, hjust = 0),
    legend.title = element_text(size = 14, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
    legend.key = element_rect(size = 1, colour = NA),
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank())

ggsave(filename = "PlotA.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")


anova_func <- function(df){
  sites_aov <- aov(no_SE ~ current * year * month, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = SE_winds))

noSE_year <- SE_winds %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(no_SE ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
noSE_month <- SE_winds %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(no_SE ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

ggplot(data = SE_monthly, aes(x = year, y = mean_SLP)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs)) +
  labs(x = "Year", y = "Sea level pressure") +
  theme(strip.text = element_text(face="bold", size=12)) +
  theme_Publication()


# Determine wind duration-------------------------------------------------------------------------------------------------------------------

load("data/CC_coastal_SLP.RData") 
load("data/CalC_coastal_SLP.RData")
load("data/BC_coastal_SLP.RData")
load("data/HC_coastal_SLP.RData")

HC_coastal_SLP <- HC_coastal_SLP %>% 
  arrange(date)

wind_dur_func <- function(df){
  wind<- df %>% 
  select(date, wind_dir_from, lat, lon) %>% 
  rename(date = date,
         wind_se = wind_dir_from)
}

BC_wind_dur <- wind_dur_func(df = BC_coastal_SLP)
CC_wind_dur <- wind_dur_func(df = CC_coastal_SLP)
CalC_wind_dur <- wind_dur_func(df = CalC_coastal_SLP)
HC_wind_dur <- wind_dur_func(df = HC_coastal_SLP)

SE_wind_func <- function(df){
  SE<- df %>% 
  filter(wind_dir_from >= 180, wind_dir_from <= 270) %>% 
  unique() %>% 
  select(date,wind_dir_from,lat,lon)
}

BC_SE <- SE_wind_func(df = BC_coastal_SLP)
CC_SE <- SE_wind_func(df = CC_coastal_SLP)
CalC_SE <- SE_wind_func(df = CalC_coastal_SLP)
HC_SE<- SE_wind_func(df = HC_coastal_SLP)


BC_prep <- right_join(BC_SE, BC_wind_dur)
BC_prep[is.na(BC_prep)] <- 0

CalC_prep <- right_join(CalC_SE, CalC_wind_dur)
CalC_prep[is.na(CalC_prep)] <- 0

CC_prep <- right_join(CC_SE, CC_wind_dur)
CC_prep[is.na(CC_prep)] <- 0

HC_prep <- right_join(HC_SE, HC_wind_dur)
HC_prep[is.na(HC_prep)] <- 0


dur_prep <- function(df){
  dur <- df %>% 
  rename(temp = wind_dir_from,
         t = date) %>% 
  select(temp,t)
}

BC_dur <- dur_prep(df = BC_prep)
CC_dur <- dur_prep(df = CC_prep)
CalC_dur <- dur_prep(df = CalC_prep)
HC_dur <- dur_prep(df = HC_prep)

HC_dur <- HC_dur %>% 
  arrange(t)

exc_BC <- exceedance(BC_dur, minDuration = 1, threshold = 0)
exc_CC <- exceedance(CC_dur, minDuration = 1, threshold = 0)
exc_CalC <- exceedance(CalC_dur, minDuration = 1, threshold = 0)
exc_HC <- exceedance(HC_dur, minDuration = 1, threshold = 0)

wind_func <- function(df){
  wind_duration <- df$exceedance %>%
  ungroup() %>%
  select(exceedance_no, duration, date_start, date_peak, intensity_max, intensity_cumulative) 
}

BC_dur <- wind_func(df = exc_BC)
CC_dur <- wind_func(df = exc_CC)
CalC_dur <- wind_func(df = exc_CalC)
HC_dur <- wind_func(df = exc_HC)

# Seasons for the southern hemisphere
seasons_S_func <- function(df){
  df_seasons <- df %>% 
    mutate(month = month(date_start, abbr = T, label = T),
           year = year(date_start)) %>% 
    mutate(season = case_when(month %in% c("Dec", "Jan", "Feb") ~ "Summer", 
                              month %in% c("Mar", "Apr", "May") ~ "Autumn",
                              month %in% c("Jun", "Jul", "Aug") ~ "Winter",
                              month %in% c("Sep", "Oct", "Nov") ~ "Spring"))
}

BC_wind<- seasons_S_func(df = BC_dur)
BC_wind <- BC_wind %>% 
  mutate(current = "BC")
HC_wind <- seasons_S_func(df = HC_dur)
HC_wind <- HC_wind %>% 
  mutate(current = "HC")

# Seasons for the Northern Hemisphere
seasons_N_func <- function(df){
  df_seasons <- df %>% 
    mutate(month = month(date_start, abbr = T, label = T),
           year = year(date_start)) %>% 
    mutate(season = case_when(month %in% c("Dec", "Jan", "Feb") ~ "Winter", 
                              month %in% c("Mar", "Apr", "May") ~ "Spring",
                              month %in% c("Jun", "Jul", "Aug") ~ "Summer",
                              month %in% c("Sep", "Oct", "Nov") ~ "Autumn"))
}

CC_wind <- seasons_N_func(df = CC_dur)
CC_wind <- CC_wind %>% 
  mutate(current = "CC")
CalC_wind <- seasons_N_func(df = CalC_dur)
CalC_wind <- CalC_wind %>% 
  mutate(current = "CalC")

duration_wind_currents <- rbind(BC_wind,CC_wind,CalC_wind,HC_wind)
# save(duration_wind_currents, file = "data/duration_wind_currents.RData")
load("data/duration_wind_currents.RData")

wind_currents <- duration_wind_currents %>% 
  filter(season == "Summer") %>% 
  group_by(year, month, current) %>% 
  summarise(mean_dur = mean(duration))

wind_currents <- as.data.frame(wind_currents)


wind_currents$month <- as.factor(wind_currents$month)
wind_currents$month <- factor(wind_currents$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
wind_currents$current = factor(wind_currents$current, levels=c('BC','HC','CalC','CC'))

plotB <- ggplot(data = wind_currents, aes(x = year, y = mean_dur)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "Duration of SE winds
(Days)")+
  theme_bw() +
  labs(colour = "Month") +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=14, family = "Palatino"),
    axis.title = element_text(size = 20, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 20, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 15, hjust = 0),
    legend.title = element_text(size = 14, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
   legend.key = element_rect(size = 1, colour = NA),
   legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank())

ggsave(filename = "PlotB.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")

anova_func <- function(df){
  sites_aov <- aov(mean_dur ~ current * year * month, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = wind_currents))

dur_year <- wind_currents %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_dur ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
dur_month <- wind_currents %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_dur ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Determining the number of pixels within each current -----------------------------------------------------------------------------------
# BC_pixels <- SE_renamed %>% 
#   filter(current == "BC") %>% 
#   dplyr::select(lon, lat) %>% 
#   unique()
# 
# CC_pixels <- SE_renamed %>% 
#   filter(current == "CC") %>% 
#   dplyr::select(lon, lat) %>% 
#   unique()
# 
# HC_pixels <- SE_renamed %>% 
#   filter(current == "HC") %>% 
#   dplyr::select(lon, lat) %>% 
#   unique()
# 
# CalC_pixels <- SE_renamed %>% 
#   filter(current == "CalC") %>% 
#   dplyr::select(lon, lat) %>% 
#   unique()

SE_monthly$month <- as.factor(SE_monthly$month)
SE_monthly$month <- factor(SE_monthly$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
SE_monthly$current = factor(SE_monthly$current, levels=c('BC','HC','CalC','CC'))
# Monthly mean temperature
ggplot(data = SE_monthly, aes(x = year, y = mean_temp)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "SST (°C)")+
  theme_bw() +
  labs(colour = "Month") +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=14, family = "Palatino"),
    axis.title = element_text(size = 20, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 20, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 15, hjust = 0),
    legend.title = element_text(size = 14, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
    legend.key = element_rect(size = 1, colour = NA),
    legend.key.size = unit(0.8, "cm"),
    legend.background = element_blank())


anova_func <- function(df){
  sites_aov <- aov(mean_temp ~ current * year * month, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = SE_monthly))

## Regression
# Slope in summer signals over time, years+months
count_year <- SE_monthly %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_temp ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
count_month <- SE_monthly %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_temp ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

ggplot(data = SE_monthly, aes(x = year, y = mean_SLP)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs)) +
  labs(x = "Year", y = "Sea level pressure") +
  theme(strip.text = element_text(face="bold", size=12)) +
  theme_Publication()

# No warming changes in SST over a 30 year periodwhen comparing summer seasons.

# For each current and not each pixel
CC_wind <- SE_monthly %>% 
  filter(current == "CC") %>% 
  mutate(signal = count / 106)

BC_wind <- SE_monthly %>% 
  filter(current == "BC") %>% 
  mutate(signal = count / 73)

CalC_wind <- SE_monthly %>% 
  filter(current == "CalC") %>%
  mutate(signal = count / 185)

HC_wind <- SE_monthly %>% 
  filter(current == "HC") %>% 
  mutate(signal = count / 193)

complete_wind <- rbind(CC_wind,BC_wind,CalC_wind,HC_wind)
# save(complete_wind, file = "data/complete_wind.RData")

load("data/complete_wind.RData")

# Plots
## Annual count of SE wind in Summer
## Summer month count of SE winds

# Number of SE wind
ggplot(data = complete_wind, aes(x = year, y = signal)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
 facet_wrap(~current,  labeller = labeller(current = supp.labs)) +
  labs(x = "Year", y = "Count") +
  theme(strip.text = element_text(face="bold", size=12)) +
  theme_Publication()

complete_wind$month <- as.factor(complete_wind$month)
complete_wind$month <- factor(complete_wind$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
complete_wind$current = factor(complete_wind$current, levels=c('BC','HC','CalC','CC'))

# Wind intensity
plotC <- ggplot(data = complete_wind, aes(x = year, y = circ_wspd, colour = Month)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  #ylab(expression("SE wind intensity")) +
  labs(x = "", y = "SE wind intensity")+
# (m.s^-1)") +
  # geom_smooth(aes(colour = month), method = "lm", se=FALSE, formula = my.formula) +
  # stat_poly_eq(formula = my.formula,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=14, family = "Palatino"),
    axis.title = element_text(size = 20, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 20, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 15, hjust = 0),
    legend.title = element_text(size = 14, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
   legend.key = element_rect(size = 1, colour = NA),
  legend.key.size = unit(0.8, "cm"),
  legend.background = element_blank())

ggsave(filename = "plotC.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")

combined <- ggarrange(plotA, plotB, plotC,labels = c("A", "B", "C"), ncol = 1, common.legend = TRUE,legend = "right")

ggsave(filename = "combined.jpg", plot = combined, width=180, height = 200, units = "mm",dpi = 300,  path = "figures/")

# ANOVA analyses on wind speed

anova_func <- function(df){
  sites_aov <- aov(circ_wspd ~ current * year, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = complete_wind))
#### Regression

circwspd_year <- complete_wind %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(circ_wspd ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

  # Slope in upwelling signals per month over years
circwspd_month <- complete_wind %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(circ_wspd ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Linear model -------------------------------------------------------------------------------------------------------

slope_calc <- function(df){
  df %>% 
    mutate(row_num = 1:n()) %>% 
    do(mod1 = lm(count ~ row_num, data = .),
       mod2 = lm(mean_temp ~ row_num, data = .),
       mod3 = lm(mean_temp ~ count, data = .),
       mod4 = cor(.$mean_temp, .$count, method = "pearson", use = "complete.obs")[1]) %>% 
    mutate(wind_slope = summary(mod1)$coeff[2],
           temp_slope = summary(mod2)$coeff[2],
           temp_wind_slope = summary(mod3)$coeff[2],
           temp_wind_r = mod4[1],
           temp_wind_r2 = glance(mod3)$adj.r.squared) %>%
    select(-mod1, -mod2, -mod3, -mod4) %>% 
    mutate_if(is.numeric, round, 2)
}
# glance(lm(mean_temp ~ count, data = SE_annual))
# Summer stats
complete_wind%>% 
  group_by(current, season) %>% 
  slope_calc()

# Monthly summer stats
complete_wind%>% 
  group_by(current, season, month) %>% 
  slope_calc()

# 3: ANOVA analyses ------------------------------------------------------------------------------------------------
# ANOVA analyses comparing is the number of gc(signals detected each year and each season varied over time

load("data/BC_UI_metrics_SLP.RData")
load("data/HC_UI_metrics_SLP.RData")
load("data/CC_UI_metrics_SLP.RData")
load("data/CalC_UI_metrics_SLP.RData")

# Seasons for the southern hemisphere
seasons_S_func <- function(df){
  df_seasons <- df %>% 
    mutate(month = month(date_start, abbr = T, label = T),
           year = year(date_start)) %>% 
    mutate(season = case_when(month %in% c("Dec", "Jan", "Feb") ~ "Summer", 
                              month %in% c("Mar", "Apr", "May") ~ "Autumn",
                              month %in% c("Jun", "Jul", "Aug") ~ "Winter",
                              month %in% c("Sep", "Oct", "Nov") ~ "Spring"))
}

BC_UI_metrics <- seasons_S_func(df = BC_UI_metrics_SLP)
HC_UI_metrics <- seasons_S_func(df = HC_UI_metrics_SLP)

# Seasons for the Northern Hemisphere
seasons_N_func <- function(df){
  df_seasons <- df %>% 
    mutate(month = month(date_start, abbr = T, label = T),
           year = year(date_start)) %>% 
    mutate(season = case_when(month %in% c("Dec", "Jan", "Feb") ~ "Winter", 
                              month %in% c("Mar", "Apr", "May") ~ "Spring",
                              month %in% c("Jun", "Jul", "Aug") ~ "Summer",
                              month %in% c("Sep", "Oct", "Nov") ~ "Autumn"))
}

CC_UI_metrics <- seasons_N_func(df = CC_UI_metrics_SLP)
CalC_UI_metrics <- seasons_N_func(df = CalC_UI_metrics_SLP)

BC_UI_metrics <- BC_UI_metrics %>% 
  mutate(current = "BC") 
HC_UI_metrics <- HC_UI_metrics %>% 
  mutate(current = "HC")
CC_UI_metrics <- CC_UI_metrics %>% 
  mutate(current = "CC")
CalC_UI_metrics <- CalC_UI_metrics %>% 
  mutate(current = "CalC")

combined_products <- rbind(BC_UI_metrics,HC_UI_metrics,CC_UI_metrics,CalC_UI_metrics)
#save(combined_products, file = "data/combined_products.RData")

# Total signals at each pixel
total_signals <- combined_products %>%
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year) %>% 
  group_by(lat, lon) %>% 
  summarise(y = n()) %>% 
  rename(count = y) %>% 
  data.frame() 

CC_signals <- CC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 106)
  

BC_signals <- BC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 73)

CalC_signals <- CalC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 185)


HC_signals <- HC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 193)
    
complete_signal <- rbind(CC_signals,BC_signals,CalC_signals,HC_signals)
# save(complete_signal, file = "data/complete_signal.RData")
load("data/complete_signal.RData")

summer_signal <- complete_signal %>% 
  filter(season == "Summer") #%>% 
  # group_by(year, current, month)

summer_signal$month <- as.factor(summer_signal$month)
summer_signal$month <- factor(summer_signal$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
summer_signal$current = factor(summer_signal$current, levels=c('BC','HC','CalC','CC'))

plot_1 <- ggplot(data = summer_signal, aes(x = year, y = signal, colour = Month)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "Upwelling events
(count)")+
  # geom_smooth(aes(colour = month), method = "lm", se=FALSE, formula = my.formula) +
  # stat_poly_eq(formula = my.formula,
  #              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
  #              parse = TRUE) +
  theme_set(theme_grey()) +
  theme_grey() +
  theme(panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
        panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
        #panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
        # panel.grid.major = element_line(size = 0.2, linetype = 2),
        # panel.grid.minor = element_line(colour = NA),
        strip.text = element_text(size=8, family = "Palatino"),
        axis.title = element_text(size = 9, face = "bold", family = "Palatino"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 8, colour = "black", family = "Palatino"),
        plot.title = element_text(size = 18, hjust = 0),
        legend.title = element_text(size = 10, family = "Palatino"),
        legend.text = element_text(size = 9, family = "Palatino"),
        legend.key.size = unit(0.2, "cm"))
        
ggsave(filename = "Plot1.pdf", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "pdf", path = "figures/")

ggsave(filename = "Plot1.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")


# Anova analyses to test whether or not a significant difference exist in the amount of 
# signals detected by each of the currents for each year and season

anova_func <- function(df){
  sites_aov <- aov(signal ~ current * year * month, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = summer_signal))

## Regression
# Slope in summer signals over time, years+months
summer_signal_year <- summer_signal %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(signal ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>%
  dplyr::rename(slope = estimate) %>%
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
summer_signal_month <- summer_signal %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(signal ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::glance)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  # filter(term == "year") %>%  3To use broom function comment this out
  # dplyr::rename(slope = estimate) %>% 
  # mutate(p.value = round(p.value, 4))

# 3: Linear models -------------------------------------------------------------------------------------------------------------------------------------------------
# ANOVA Analyses testing if there is a significant difference in the duration/mean intensity etc,
# between currents and seasons over a 30 year period

load("data/combined_products.RData")

# Calculate all of the linear models
lm_coeff <- function(df){
  res <- lm(formula = val ~ date_peak, data = df)
  res_coeff <- as.numeric(res$coefficient[2])
}

lm_metrics <- combined_products %>% 
  ungroup() %>% 
  dplyr::select(-c(lon:index_end), -c(month)) %>% 
  pivot_longer(cols = c(duration, intensity_mean:rate_decline), 
               names_to = "var", values_to = "val") %>% 
  group_by(current, season, year, var) %>% 
  nest() %>% 
  mutate(slope = purrr::map(data, lm_coeff)) %>% 
  dplyr::select(-data) %>% 
  unnest(cols = slope) %>% 
  # convert from daily to decadal values
  mutate(slope = round((slope*365.25*10), 2)) %>% 
  ungroup()

#save(lm_metrics, file = "data/lm_metrics.RData")


lm_metrics_wide <- pivot_wider(lm_metrics, 
                               id_cols = current:season, 
                               names_from = var, values_from = slope,
                               values_fn = mean) %>% 
  filter(season == "Summer")


summary(aov(duration ~ current + season, data = lm_metrics_wide))
summary(aov(intensity_mean ~ current + season, data = lm_metrics_wide))
summary(aov(intensity_max ~ current  + season, data = lm_metrics_wide))
summary(aov(intensity_cumulative ~ current + season, data = lm_metrics_wide))


summary(aov(duration ~ current + year, data = lm_metrics_wide))
  summary(aov(intensity_mean ~ current + year, data = lm_metrics_wide))
summary(aov(intensity_max ~ current  + year, data = lm_metrics_wide))
summary(aov(intensity_cumulative ~ year + season, data = lm_metrics_wide))

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Combined products to view for changes in the intensity of upwelling signals...

CC_int <- CC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(mean_intensity = mean(intensity_mean)) #%>% 
  #mutate(signal = y / 106)

HC_int <- HC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(mean_intensity = mean(intensity_mean)) #%>% 
  #mutate(mean_int = mean_intensity / 193)

CalC_int <- CalC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(mean_intensity = mean(intensity_mean)) #%>% 
 # mutate(signal = y / 185)

BC_int <- BC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(mean_intensity = mean(intensity_mean)) #%>% 
 # mutate(signal = y / 73)

mean_int <- rbind(BC_int,CC_int, CalC_int, HC_int)
mean_int <- mean_int %>% 
  filter(season == "Summer")

mean_int$month <- as.factor(mean_int$month)
mean_int$month <- factor(mean_int$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
mean_int$current = factor(mean_int$current, levels=c('BC','HC','CalC','CC'))

plot_2 <- ggplot(data = mean_int, aes(x = year, y = mean_intensity, colour = Month)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "Mean intensity of 
signals")+
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=8, family = "Palatino"),
    axis.title = element_text(size = 9, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 8, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 18, hjust = 0),
    legend.title = element_text(size = 10, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
    #legend.key = element_rect(size = 0.8, colour = NA),
    legend.key.size = unit(0.2, "cm"),
    legend.background = element_blank())

ggsave(filename = "Plot2.pdf", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "pdf", path = "figures/")

ggsave(filename = "Plot2.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")


anova_func <- function(df){
  sites_aov <- aov(mean_intensity ~ current * year * month, data = df)
  return(sites_aov)
}
summary(count_aov <- anova_func(df = mean_int))

mean_int_year <- mean_int %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_intensity ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
mean_int_month <- mean_int %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_intensity ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# To get R2 value

mean_int_month <- mean_int %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_intensity ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::glance)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") 


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

CC_int <- CC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(cum_intensity = mean(intensity_cumulative)) 

HC_int <- HC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(cum_intensity = mean(intensity_cumulative))

CalC_int <- CalC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(cum_intensity = mean(intensity_cumulative)) 

BC_int <- BC_UI_metrics %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(cum_intensity = mean(intensity_cumulative)) 

cum_int <- rbind(BC_int,CC_int, CalC_int, HC_int)
cum_int <- cum_int %>% 
  filter(season == "Summer")

cum_int$month <- as.factor(cum_int$month)
cum_int$month <- factor(cum_int$month, levels = c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
cum_int$current = factor(cum_int$current, levels=c('BC','HC','CalC','CC'))

plot_3 <- ggplot(data = cum_int, aes(x = year, y = cum_intensity, colour = Month)) +
  geom_line(aes(colour = month)) +
  geom_smooth(aes(colour = month), method = "lm") +
  facet_wrap(~current,  labeller = labeller(current = supp.labs), ncol = 4) +
  labs(x = "", y = "Cumulative intensity
of signals")+
  theme_set(theme_grey()) +
  theme_grey() +
  theme(#panel.border = element_rect(colour = "black", fill = NA, size = 1.0),
    # panel.grid.major = element_line(size = 0.2, linetype = 2),
    # panel.grid.minor = element_line(colour = NA),
    strip.text = element_text(size=8, family = "Palatino"),
    axis.title = element_text(size = 9, face = "bold", family = "Palatino"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major = element_line("grey70", linetype = "dashed", size = 0.2),
    panel.grid.minor = element_line("grey70", linetype = "dashed", size = 0.2),
    axis.text = element_text(size = 8, colour = "black", family = "Palatino"),
    plot.title = element_text(size = 8, hjust = 0),
    legend.title = element_text(size = 10, family = "Palatino"),
    legend.text = element_text(size = 9, family = "Palatino"),
    #legend.key = element_rect(size = 0.8, colour = NA),
    legend.key.size = unit(0.2, "cm"),
    legend.background = element_blank())

ggsave(filename = "Plot3.jpg", plot = last_plot(), width=180, height = 50,units = "mm",dpi = 300, device = "jpg", path = "figures/")


combined_upwellng <- ggarrange(plot_1, plot_2, plot_3,labels = c("A", "B", "C"), ncol = 1, common.legend = TRUE,legend = "right")
ggsave(filename = "combined_upwellng.jpg", plot = combined_upwellng, width=180, height = 200, units = "mm",dpi = 300,  path = "figures/")


anova_func <- function(df){
  sites_aov <- aov(cum_intensity ~ current * year * month, data = df)
  return(sites_aov)
}
summary(count_aov <- anova_func(df = cum_int))

cum_int_year <- cum_int %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(cum_intensity ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Slope in upwelling signals per month over years
cum_int_month <- cum_int %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(cum_intensity ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

