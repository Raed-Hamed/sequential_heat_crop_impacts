#======================================================================
#clean work space
#======================================================================
rm(list = ls())
graphics.off()
gc()
#======================================================================
#load libraries
#======================================================================
library(tidyverse)
library(sf)
library(raster)
library(colorspace)
library(spData)
library(tools)
library(patchwork)
source("00_load_functions.R")
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_model <- file.path(root_dir, "data/model_data")
dir_cmip6 <- file.path(root_dir, "data/cmip6/")  
#======================================================================
#load helper function
#======================================================================
# Define the prediction function based on the model equation
calc_yield_predictions <- function(data) {
  
  coef_intercept <- data$`(Intercept)`
  coef_tmx_spring <- data$mean_tmax_spring
  coef_tmx_summer <- data$mean_tmax_summer
  coef_sequential_heat <- data$`mean_tmax_spring:mean_tmax_summer`
  
  data$predicted_yield <- coef_intercept +
    coef_tmx_spring * data$delta_spring_tx  +
    coef_tmx_summer * data$delta_summer_tx  +
    coef_sequential_heat * data$delta_spring_tx * data$delta_summer_tx
  
  data$tx_spring_effect <- 
    coef_tmx_spring * data$delta_spring_tx
  
  data$tx_summer_effect <- 
    coef_tmx_summer * data$delta_summer_tx
  
  data$sequential_heat_effect <- 
    coef_sequential_heat  * data$delta_spring_tx * data$delta_summer_tx
  
  return(data)
}
#----------------------------------------------------------------------
#load summary bootstrap coefficients
load_summary_boot_coef_data <- function(filename) {
  coef_summary <-
    list.files(dir_model,
               filename ,
               full.names = TRUE) %>%
    readRDS()  %>%
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
#load CMIP6 data
#======================================================================
#load TX delta changes
cmip6_tx_delta_change_per_model <-
  list.files(dir_cmip6,
             "counties_difference_per_model.shp",
             full.names = TRUE) %>%
  read_sf() %>%
  tibble() %>%
  dplyr::select(-geometry) %>% 
  pivot_longer("X245_su_M0":"X119_sp_mea") %>%
  dplyr::select(NAME_0,NAME_1,NAME_2,name,value) %>% 
  separate(name, c("model_scenario", "season","model_ref")) %>%
  mutate(across('season', str_replace, 'su', 'delta_summer_tx')) %>%
  mutate(across('season', str_replace, 'sp', 'delta_spring_tx')) %>%
  filter(!NAME_1  %in% c("Alaska", "Hawaii"))%>%
  mutate_at(vars(NAME_1, NAME_2), str_to_upper) %>% 
  pivot_wider(names_from = season, values_from = value) %>% 
  unite("group_id",NAME_1 :NAME_2) %>% 
  mutate(delta_spring_sm = 0,
         delta_summer_sm= 0)
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
  mutate(crop = "corn")
#----------------------------------------------------------------------
#process bootstrap corn coefficients
boot_soy_coefs <-boot_soy_model %>% 
  dplyr::select(State,County,coef,mean) %>% 
  pivot_wider(names_from = coef, values_from = mean)  %>%
  mutate(crop = "soy")
#----------------------------------------------------------------------
boot_coefs <-bind_rows(boot_corn_coefs,boot_soy_coefs) %>% 
  unite(group_id, c("State","County"))
#======================================================================
#calculate yield projections per temperature driver (tx predictor)
#======================================================================
#calculate crop yield projections
calc_yield_predictions_df <-
  inner_join(cmip6_tx_delta_change_per_model,
             boot_coefs) %>%
  dplyr::select(-NAME_0) %>%
  group_by(model_scenario, model_ref, group_id, crop) %>%
  nest() %>%
  summarise(data = map(data, calc_yield_predictions)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct()
#----------------------------------------------------------------------
#save crop yield projections for all scenarios and models per scenario
saveRDS(calc_yield_predictions_df,
        file.path(dir_cmip6,
                  "crop_yield_projections.rds"))
