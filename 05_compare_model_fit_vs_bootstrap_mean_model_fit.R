#======================================================================
#clean work space
#======================================================================
rm(list = ls())
graphics.off()
gc()
#======================================================================
#load libraries
#======================================================================
library(sf)
library(tidyverse)
library(spData)
library(raster)
library(lubridate)
library(colorspace)
library(lme4)
library(piecewiseSEM)
library(lmerTest)
library(broom.mixed)
library(patchwork)
source("00_load_functions.R")
#======================================================================
#load helper functions
#======================================================================
load_summary_boot_coef_data <- function(filename) {
  coef_summary <-
    list.files(dir_model,
               filename ,
               full.names = TRUE) %>%
    readRDS()  %>%
   # filter(coef != "(Intercept)") %>%
    drop_na() %>%
    group_by(group_id, coef) %>%
    summarise(
      mean = mean(value),
      median = median(value),
      percentile_5th = quantile(value, 0.05),
      percentile_95th = quantile(value, 0.95)
    ) %>%
    mutate(significant = as.integer(sign(percentile_5th) == sign(percentile_95th))) %>%
    separate(group_id, c("State", "County"), sep = "_")
  
  return(coef_summary)
  
}
#======================================================================
#load helper function
#======================================================================
#Define the prediction function based on the model equation
calc_yield_predictions <- function(data) {
  
  coef_intercept <- data$`coef_(Intercept)`
  coef_tmx_spring <- data$coef_mean_tmax_spring
  coef_tmx_summer <- data$coef_mean_tmax_summer
  coef_tmx_interaction <- data$`coef_mean_tmax_spring:mean_tmax_summer`
  coef_sm_spring <- data$coef_mean_sm_spring
  coef_sm_summer <- data$coef_mean_sm_summer
  coef_sm_interaction <- data$`coef_mean_sm_spring:mean_sm_summer`
  
  
  data$predicted_yield <- coef_intercept +
    coef_tmx_spring * data$mean_tmax_spring  +
    coef_tmx_summer * data$mean_tmax_summer +
    coef_tmx_interaction * data$mean_tmax_spring * data$mean_tmax_summer +
    coef_sm_spring * data$mean_sm_spring  +
    coef_sm_summer * data$mean_sm_summer +
    coef_sm_interaction * data$mean_sm_spring * data$mean_sm_summer 
  
  return(data)
}
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
dir_data <- file.path(dir_root, "data")
dir_echam <- file.path(dir_data, "echam_climatology")
dir_storyline <- file.path(dir_data, "echam_storylines")
dir_obs_climate <- file.path(dir_data, "obs_climate_data")
dir_crop <- file.path(dir_data, "crop_data")
dir_shp <- file.path(dir_data, "shapefiles")
dir_model <- file.path(dir_data, "model_data")
dir_figures <- file.path(dir_root, "figures")
#======================================================================
#Load spatial data
#======================================================================
#spatial info from crop data
spatial_info <-
  list.files(dir_model,
             "spatial_county_group_id_ref.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  tibble() %>%
  dplyr::select(-group_id) %>%
  distinct() 
#----------------------------------------------------------------------
#state spatial info US
us_state_shp <- list.files(dir_shp, ".json", full.names = TRUE) %>%
  map_dfr(st_read) %>%
  mutate(NAME_1  = NAME_1 %>%  str_to_upper) %>%
  filter(!NAME_1 %in%  c("ALASKA", "HAWAII"))
#----------------------------------------------------------------------
#state spatial info EU
eu_state_shp <-
  list.files(dir_shp, "NUTS_RG_20M_2021_3035.shp", full.names = TRUE) %>%
  st_read() %>%
  filter(LEVL_CODE == 2) %>%
  st_transform(., crs = st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>%
  filter(CNTR_CODE %in% spatial_info$State)
#======================================================================
#load bootstrap coefficient summary per crop incl. CIs (Tx in degree C)
#======================================================================
#load coefficient bootstrap for corn model
boot_corn_model <-
  load_summary_boot_coef_data("coef_ci_corn.rds")
#----------------------------------------------------------------------
#load coefficient bootstrap for soy model
boot_soy_model <- load_summary_boot_coef_data("coef_ci_soy.rds")
#----------------------------------------------------------------------
#process bootstrap corn coefficients
boot_corn_coefs <-boot_corn_model %>% 
  dplyr::select(State,County,coef,mean) %>% 
  pivot_wider(names_from = coef, values_from = mean)  %>%
  rename_with(~ paste0("coef_", .),-all_of(c("State","County"))) %>% 
  mutate(crop = "corn")
#----------------------------------------------------------------------
#process bootstrap corn coefficients
boot_soy_coefs <-boot_soy_model %>% 
  dplyr::select(State,County,coef,mean) %>% 
  pivot_wider(names_from = coef, values_from = mean)  %>%
  rename_with(~ paste0("coef_", .),-all_of(c("State","County"))) %>% 
  mutate(crop = "soy")
#----------------------------------------------------------------------
boot_coefs <-bind_rows(boot_corn_coefs,boot_soy_coefs)
#======================================================================
#load model ceofficients
#======================================================================
#load model incl. all predictors
model_original_coefs <-
  list.files(dir_model,
             "_full" ,
             full.names = TRUE) %>%
  map(~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_)(.*?)(?=\\_full.rds)")
    data <- readRDS(.)
    tibble(crop = crop, model_fit = list(data))
  }) %>%
  bind_rows() %>%
  mutate(model_coefs = map(
    model_fit,
    ~ coef(.) %>%
      .$group_id %>%
      rownames_to_column("group_id")
    
  )) %>%
  dplyr::select(-model_fit) %>%
  unnest(c(model_coefs))  %>%
  rename_with(~ paste0("coef_", .),-all_of(c("crop", "group_id"))) %>% 
  separate(group_id, c("State", "County"), sep = "_")
#===============================================================================
#load model original dataset
#===============================================================================
model_data <-
  list.files(dir_model, "model_data_", full.names = TRUE) %>%
  map_dfr( ~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_data_)(.*?)(?=\\.rds)")
    data <- readRDS(.)
    tibble(crop = crop, data)
  }) %>%
  dplyr::select(-geometry, -group_id) %>%
  drop_na() %>%
  filter(n >= 15) %>%
  group_by(County, State, zone, crop) %>%
  mutate_at(
    vars(
      crop_yield,
      mean_sm_spring,
      mean_sm_summer,
      mean_tmax_spring,
      mean_tmax_summer
    ),
    ~ scale(., center = TRUE, scale = FALSE)
  )  %>% 
  mutate_at(
    vars(
      # crop_yield,
      mean_sm_spring,
      mean_sm_summer,
      # mean_tmax_spring,
      #mean_tmax_summer
    ),
    ~ scale(., center = FALSE, scale = TRUE)
  ) 
#===============================================================================
#Calculate R-square based on both model coefficients and bootsrap mean coefs
#===============================================================================
county_scale_rsq_original_model <-inner_join(model_data,
                                             model_original_coefs) %>% 
  group_by(crop,State,County) %>% 
  nest() %>%
  summarise(data = map(data, calc_yield_predictions)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct() %>% 
  group_by(crop,State,County) %>% 
  summarise(rsq = cor(crop_yield,predicted_yield)^2) %>% 
  ungroup() %>%  
  distinct() %>% 
  inner_join(spatial_info)
#----------------------------------------------------------------------
county_scale_rsq_bootstrap_model <-inner_join(model_data,
                                              boot_coefs) %>% 
  group_by(crop,State,County) %>% 
  nest() %>%
  summarise(data = map(data, calc_yield_predictions)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct() %>% 
  group_by(crop,State,County) %>% 
  summarise(rsq = cor(crop_yield,predicted_yield)^2) %>% 
  ungroup() %>%  
  distinct() %>% 
  inner_join(spatial_info)
#----------------------------------------------------------------------
#Plot R-square
county_scale_rsq_original_model %>%
    ungroup() %>% 
    st_as_sf() %>%
    ggplot(aes(fill = rsq)) +
    geom_sf() +
    geom_sf(
      data = us_state_shp,
      color = "black",
      size = 0.00000001,
      fill = "transparent"
    ) +
    theme_void(base_size = 5) +
    scale_fill_continuous_sequential(name = "", "reds 3", limit = c(0,1)) +
    facet_grid(~crop)+
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1, 'cm'),
      legend.title.align = 0.5,
      legend.title = element_text(size = 15, angle = 0),
      legend.key.width = unit(5, "cm"),
      legend.box = "horizontal",
      plot.title = element_text(size = 20, hjust = 0.5),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      legend.text = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
      
    ) +
    guides(fill = guide_colorbar(title.position = "bottom"))
  