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
library(broom.mixed)
source("00_load_functions.R")
library(piecewiseSEM)
library(patchwork)
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
dir_data <- file.path(dir_root, "data")
dir_model <- file.path(dir_data, "model_data")
dir_shp <- file.path(dir_data, "shapefiles")
dir_figures <- file.path(dir_root, "figures")
#======================================================================
#Get quantile equivalents for variables of interest
#======================================================================
#load model data 
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
  ) %>% 
  unite("group_id",State:County) %>% 
  unite(crop, c(crop, zone), sep = "_") %>%
  mutate(
    crop = case_when(
      crop == "corn_US" ~ "corn",
      crop == "soy_US" ~ "soy",
      crop == "wheat_US" ~ "wheat_us",
      crop == "wheat_EU" ~ "wheat_eu",
      TRUE ~ crop
    )
  )
#----------------------------------------------------------------------
#store spring tx quantile values
spring_tx_quantile_values <-
  model_data %>%
  group_by(crop) %>%
  summarise_at(vars(c(mean_tmax_spring)),
               ~ data.frame(
                 q05 = quantile(., 0.05),
                 q50 = quantile(., 0.50),
                 q95 = quantile(., 0.95)
               )) %>%
  mutate(q05 = mean_tmax_spring$q05,
         q50 = mean_tmax_spring$q50,
         q95 = mean_tmax_spring$q95) %>%
  dplyr::select(-mean_tmax_spring)
#----------------------------------------------------------------------
#store summer tx quantile values
summer_tx_quantile_values <-
  model_data %>%
  group_by(crop) %>%
  summarise_at(vars(c(mean_tmax_summer)),
               ~ data.frame(
                 q05 = quantile(., 0.05),
                 q50 = quantile(., 0.50),
                 q95 = quantile(., 0.95)
               )) %>%
  mutate(q05 = mean_tmax_summer$q05,
         q50 = mean_tmax_summer$q50,
         q95 = mean_tmax_summer$q95) %>%
  dplyr::select(-mean_tmax_summer)
#======================================================================
#load model fit and extract key features
#======================================================================
#load model incl. all predictors
model_fit_full <-
  list.files(dir_model,
             "_full" ,
             full.names = TRUE) %>%
  map(~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_)(.*?)(?=\\_full.rds)")
    data <- readRDS(.)
    tibble(crop = crop, model_fit = list(data))
  }) %>%
  bind_rows() 
#----------------------------------------------------------------------
#get model r-squared values for full model
model_fit_full %>% 
  mutate(map_dfr(model_fit, ~ rsquared(.)))
#----------------------------------------------------------------------
#load model incl. all predictors minus sequential interaction terms (to compare)
model_fit_simple <-
  list.files(dir_model,
             "_interaction" ,
             full.names = TRUE) %>%
  map( ~ {
    filename <- basename(.)
    crop <-
      str_extract(filename,
                  "(?<=model_)(.*?)(?=\\_no_sequential_interaction.rds)")
    data <- readRDS(.)
    tibble(crop = crop, model_fit = list(data))
  }) %>%
  bind_rows() 
#----------------------------------------------------------------------
#get model r-squared values for simple model (to compare)
model_fit_simple %>% 
  mutate(map_dfr(model_fit, ~ rsquared(.)))
#----------------------------------------------------------------------
#check model fixed coefficients
model_fixed_coefs <- model_fit_full %>%
  mutate(coefs = map(model_fit,
                     ~ (.) %>%
                       tidy())) %>%
  unnest(coefs) %>%
  filter(
    term %in% c(
      "mean_tmax_spring",
      "mean_tmax_summer",
      "mean_tmax_spring:mean_tmax_summer"
    )
  ) %>% 
  dplyr::select(crop,term,estimate) %>% 
  pivot_wider(names_from = crop, values_from =estimate)
#----------------------------------------------------------------------
#store tx coefficients of interest
model_tx_coefs <- model_fit_full %>% 
  mutate(summer_coef = map(
    model_fit,
    ~ (.) %>%
      tidy() %>%
      filter(term == "mean_tmax_summer") %>%
      pull(estimate)
  ))%>% 
  mutate(spring_coef = map(
    model_fit,
    ~ (.) %>%
      tidy() %>%
      filter(term == "mean_tmax_spring") %>%
      pull(estimate)
  ))%>% 
  mutate(interaction_coef = map(
    model_fit,
    ~ (.) %>%
      tidy() %>%
      filter(term == "mean_tmax_spring:mean_tmax_summer") %>%
      pull(estimate)
  )) %>% 
  dplyr::select(-model_fit) %>% 
  unnest(c(spring_coef,summer_coef,interaction_coef))
#======================================================================
#join model coefficients and quantile spring TX values
#======================================================================
calc_tx_spring_effect_on_summer_tx_impact <-
  inner_join(spring_tx_quantile_values, model_tx_coefs) %>%
  mutate(
    summer_effect_cold_spring = summer_coef + q05 * interaction_coef,
    summer_effect_hot_spring = summer_coef + q95 * interaction_coef,
    summer_effect_mild_spring = summer_coef + q50 * interaction_coef
  ) %>%
  mutate(additional_hot_spring_impacts =
           summer_effect_hot_spring / summer_effect_mild_spring)
#======================================================================
#join model coefficients and quantile spring TX values
#======================================================================
calc_tx_summer_effect_on_spring_tx_impact <-
  inner_join(summer_tx_quantile_values, model_tx_coefs) %>%
  mutate(
    spring_effect_hot_summer = spring_coef + q95 * interaction_coef,
    spring_effect_mild_summer = spring_coef + q50 * interaction_coef
  ) %>%
  mutate(additional_hot_summer_impacts =
           spring_effect_hot_summer / spring_effect_mild_summer)

