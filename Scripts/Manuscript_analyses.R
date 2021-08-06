# Statistical analyses
# The purpose of this script is to...
# The steps taken are:
# 1: Setup environment
# 2: Loading the datasets produced in Wind_data.R and Upwelling&SST.R script

# 1: Setup environment ------------------------------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(heatwaveR)
## devtools::install_github("robwschlegel/coastR")
library(mgcv) # for gam
library(doParallel); registerDoParallel(cores = 4)
source("functions/theme.R")
options(scipen = 999) 

# 2: Loading the data --------------------------------------------------------------------------------------------------------

# Loading the wind data
load("data_official/winds.RData")

wind_means <- winds %>% 
  group_by(current) %>% 
  summarise(mean_duration = mean(duration_mean))

# Loading the temperature data4
load("data_official/temp_monthly.RData")
  
anova_func <- function(df){
  sites_aov <- aov(mean_temp ~ current * year * month, data = df)
  return(sites_aov)
}

summary(count_aov <- anova_func(df = temp_monthly))

temp <- temp_monthly %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_temp ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% #change glance to rtidy
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

count_month <- temp_monthly %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(mean_temp ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Loading the upwelling metrics
load("data_official/upwell_south_BC.RData")
load("data_official/upwell_north_BC.RData")
load("data_official/upwell_Canary_current.RData")
load("data_official/upwell_chile.RData")
load("data_official/upwell_peru.RData")
load("data_official/upwell_south_CalC.RData")
load("data_official/upwell_north_CalC.RData")

current_upwelling <- rbind(upwell_south_BC, upwell_north_BC, upwell_Canary_current,
                           upwell_chile, upwell_peru, upwell_south_CalC, upwell_north_CalC)

upwelling_metrics <- current_upwelling %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(mean_intensity = mean(intensity_mean),
            cum_intensity = mean(intensity_cumulative)) 

# Observing changes in upwelling metrics, wind metrics and temperature overtime. 

# 3: Some analyses ------------------------------------------------------------------------------------------------------------
# Duration
anova_func <- function(df){
  sites_aov <- aov(duration_mean ~ current * year * month, data = df)
  return(sites_aov)
}

summary(dur_aov <- anova_func(df = winds))

wind_dur <- winds %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(duration_mean ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Replace glance with tidy
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))
# - lowering in value means more intense
winds_month <- winds %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(duration_mean ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::glance)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Upwelling favourable wind event count
anova_func <- function(df){
  sites_aov <- aov(event_count ~ current * year * month, data = df)
  return(sites_aov)
}

summary(event_aov <- anova_func(df = winds))

wind_event_count <- winds %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(event_count ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

winds_month <- winds %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(event_count ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Upwelling favourable wind intensity
anova_func <- function(df){
  sites_aov <- aov(intensity_mean ~ current * year * month, data = df)
  return(sites_aov)
}

summary(intensity_aov <- anova_func(df = winds))

wind_intensity <- winds %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(intensity_mean ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

winds_month <- winds %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(intensity_mean ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::glance)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

# Analysing trends in upwelling metrics

load("data_official/current_upwelling.RData")

current_upwelling <- as.data.frame(current_upwelling)

# Calculate all of the linear models
lm_coeff <- function(df){
  res <- lm(formula = val ~ date_peak, data = df)
  res_coeff <- as.numeric(res$coefficient[2])
}

lm_metrics <- current_upwelling %>% 
  ungroup() %>% 
  drop_na() %>% 
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


##### Upwelling metrics
load("data_official/current_upwelling.RData")

current_upwelling <- current_upwelling %>% 
  filter(season == "Summer")

BC_S_signals <- upwell_south_BC %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 45) # Here we are dividing by the number of pixels within this region

BC_N_signals <- upwell_north_BC %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 16)

CC_signals <- upwell_Canary_current %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 82)

CalC_S_signals <- upwell_south_CalC %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 24)

CalC_N_signals <- upwell_north_CalC %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 20)

Peru_signals <- upwell_peru %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 25)

Chile_signals <- upwell_chile %>% 
  mutate(year = year(date_start)) %>% 
  group_by(current, season,year, month) %>% 
  summarise(y = n()) %>% 
  mutate(signal = y / 99)

complete_signal <- rbind(BC_S_signals, BC_N_signals, CC_signals, CalC_N_signals,
                         CalC_S_signals, Peru_signals, Chile_signals)

anova_func <- function(df){
  sites_aov <- aov(intensity_mean ~ current * year * month, data = df)
  return(sites_aov)
}

summary(intensity_aov <- anova_func(df = current_upwelling))

upwell_intensity <- current_upwelling %>% 
  group_by(current) %>% 
  mutate(year_month = 1:n()) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(intensity_mean ~ year_month, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year_month") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))

upwell_month <- current_upwelling %>% 
  group_by(current, month) %>% 
  nest() %>% 
  mutate(model_out = purrr::map(data, ~lm(intensity_mean ~ year, data = .)),
         model_tidy = purrr::map(model_out, broom::tidy)) %>% # Using 'broom::glance' will provide the R2 value 
  dplyr::select(-data, -model_out) %>% 
  unnest(cols = "model_tidy") %>% 
  filter(term == "year") %>% 
  dplyr::rename(slope = estimate) %>% 
  mutate(p.value = round(p.value, 4))
  
summer_signal_year <- summer_signal_new %>% 
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
