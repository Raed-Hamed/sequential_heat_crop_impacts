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
library(lme4)
library(parallel)
library(sf)
library(piecewiseSEM)
#======================================================================
#load directories and key helper functions
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_crop <- file.path(root_dir, "data/usda-nass")
dir_eurostat <- file.path(root_dir, "data/eurostat")
dir_model <- file.path(root_dir, "data/model_data")
#======================================================================
model_df <-
  list.files(dir_model, "model_data_", full.names = TRUE) %>%
  map_dfr(~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_data_)(.*?)(?=\\.rds)")
    data <- readRDS(.)
    tibble(crop = crop, data)
  }) %>%
  dplyr::select(-geometry,-group_id) %>%
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
  ) %>% 
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
  unite("group_id",State:County) 
#======================================================================
#store spatial info separately for easy plotting later on
spatial_info <- model_df %>%
  dplyr::select(crop, zone, State, County, geometry, group_id) %>%
  distinct()
#======================================================================
model_soy <- model_df %>%
  ungroup() %>%
  filter(crop == "soy")

model_corn <- model_df %>%
  ungroup() %>%
  filter(crop == "corn")

model_wheat_us <- model_df %>%
  ungroup() %>%
  filter(crop == "wheat") %>%
  filter(zone == "US")

model_wheat_eu <- model_df %>%
  ungroup() %>%
  filter(crop == "wheat") %>%
  filter(zone == "EU")
#======================================================================
mixed_model_soy <- lmer(
  crop_yield ~
    mean_sm_spring + mean_sm_summer +
    mean_tmax_spring + mean_tmax_summer +
    mean_tmax_spring:mean_sm_spring+
    mean_sm_spring:mean_sm_summer +
    mean_tmax_spring:mean_tmax_summer +
    mean_tmax_summer:mean_sm_summer +
     (1 + mean_sm_spring +
       mean_sm_summer +
       mean_tmax_spring +
       mean_tmax_summer +
       mean_tmax_spring:mean_sm_spring+
       mean_sm_spring:mean_sm_summer +
       mean_tmax_spring:mean_tmax_summer +
       mean_tmax_summer:mean_sm_summer
     | group_id
    )
  ,
  data = model_soy
  ,
  control = lmerControl(optimizer = "bobyqa",
                          optCtrl = list(maxfun = 1e7))
)
#======================================================================
mixed_model_corn <- lmer(
  crop_yield ~
    mean_sm_spring + mean_sm_summer +
    mean_tmax_spring + mean_tmax_summer +
    mean_tmax_spring:mean_sm_spring+
    mean_sm_spring:mean_sm_summer +
    mean_tmax_spring:mean_tmax_summer +
    mean_tmax_summer:mean_sm_summer +
    (1 + mean_sm_spring +
       mean_sm_summer +
       mean_tmax_spring +
       mean_tmax_summer +
       mean_tmax_spring:mean_sm_spring+
       mean_sm_spring:mean_sm_summer +
       mean_tmax_spring:mean_tmax_summer +
       mean_tmax_summer:mean_sm_summer
     | group_id
    ),
  
  data = model_corn
  ,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 1e7))
)
#======================================================================
mixed_model_wheat_us <- lmer(
  crop_yield ~
    mean_sm_spring + mean_sm_summer +
    mean_tmax_spring + mean_tmax_summer +
    mean_tmax_spring:mean_sm_spring+
    mean_sm_spring:mean_sm_summer +
    mean_tmax_spring:mean_tmax_summer +
    mean_tmax_summer:mean_sm_summer +
    (1 + mean_sm_spring +
       mean_sm_summer +
       mean_tmax_spring +
       mean_tmax_summer +
       mean_tmax_spring:mean_sm_spring+
       mean_sm_spring:mean_sm_summer +
       mean_tmax_spring:mean_tmax_summer +
       mean_tmax_summer:mean_sm_summer
     | group_id
    ),
  
  data = model_wheat_us
  ,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 1e7))
)
#======================================================================
mixed_model_wheat_eu <- lmer(
  crop_yield ~
    mean_sm_spring + mean_sm_summer +
    mean_tmax_spring + mean_tmax_summer +
    mean_tmax_spring:mean_sm_spring+
    mean_sm_spring:mean_sm_summer +
    mean_tmax_spring:mean_tmax_summer +
    mean_tmax_summer:mean_sm_summer +
    (1 + mean_sm_spring +
       mean_sm_summer +
       mean_tmax_spring +
       mean_tmax_summer +
       mean_tmax_spring:mean_sm_spring+
       mean_sm_spring:mean_sm_summer +
       mean_tmax_spring:mean_tmax_summer +
       mean_tmax_summer:mean_sm_summer
     | group_id
    ),
  
  data = model_wheat_eu
  ,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 1e7))
)
#======================================================================
#check model assumptions
#======================================================================
#get marginal and conditional r-squared per model
# mixed_model_corn %>% rsquared(.)
# mixed_model_soy  %>% rsquared(.)
# #store fitted and residual corn model values
# fitted_vals_corn  <- fitted.values(mixed_model_fit_corn)
# residuals_corn = residuals(mixed_model_fit_corn)
# #store fitted and residual soy model values
# fitted_vals_soy  <- fitted.values(mixed_model_soy)
# residuals_soy = residuals(mixed_model_soy)
# #======================================================================
# #create corn diagnostic df
# diagnostics_corn <-
#   tibble(fitted_vals =
#            fitted_vals_corn,
#          residuals  = residuals_corn)
# #create soy diagnostic df
# diagnostics_soy <-
#   tibble(fitted_vals =
#            fitted_vals_soy,
#          residuals  = residuals_soy)
# #======================================================================
# #check homogeneity of variances
# diagnostics_soy %>%
#   ggplot(aes(x = fitted_vals, y = residuals)) +
#   geom_point(shape = ".") +
#   theme_bw()
# #plot residual histogram (check normality)
# diagnostics_corn %>% 
#   ggplot(aes(x=residuals))+
#   geom_density(fill = "red", alpha = 0.7)
#======================================================================
#save models
#======================================================================
#save corn model
saveRDS(
  mixed_model_corn,
  file.path(
    dir_model,
    "model_corn_full.rds"
  )
)
#save soy model
saveRDS(
  mixed_model_soy,
  file.path(
    dir_model,
    "model_soy_full.rds"
  ))

#save wheat us model
saveRDS(
  mixed_model_wheat_us,
  file.path(
    dir_model,
    "model_wheat_us_full.rds"
  ))

#save wheat eu model
saveRDS(
  mixed_model_wheat_eu,
  file.path(
    dir_model,
    "model_wheat_eu_full.rds"
  ))

#store county group ID reference 
saveRDS(spatial_info,
        file.path(dir_model, "spatial_county_group_id_ref.rds"))