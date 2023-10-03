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
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
#======================================================================
#load key data sets
#======================================================================
model_df <-
  list.files(dir_root, "model_data_soy.rds", full.names = TRUE) %>%
  readRDS() %>% 
  dplyr::select(-geometry,-group_id) %>%
  drop_na() %>%
  filter(n >= 15) %>%
  group_by(County, State, zone) %>%
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
      mean_sm_spring,
      mean_sm_summer
    ),
    ~ scale(., center = FALSE, scale = TRUE)
  ) %>% 
  unite("group_id",State:County) 
#======================================================================
#single fit mixed effect model 
#======================================================================
#southeast south america model fit
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
    ),
  
  data = model_soy
  ,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 1e7))
)
#======================================================================
#load CMIP6 data
#======================================================================
#load TX delta changes
cmip6_tx_delta_change_data <-
  list.files(dir_root,
             "counties_difference_per_model.shp",
             full.names = TRUE) %>%
  read_sf() %>%
  tibble() %>%
  dplyr::select(-geometry) %>% 
  pivot_longer("245_su_M0":"119_sp_mea") %>%
  dplyr::select(NAME_0,NAME_1,NAME_2,name,value) %>% 
  separate(name, c("model_scenario", "season","model_ref")) %>%
  mutate(across('season', str_replace, 'su', 'delta_summer_tx')) %>%
  mutate(across('season', str_replace, 'sp', 'delta_spring_tx')) %>%
  filter(!NAME_1  %in% c("Alaska", "Hawaii")) %>% 
  filter(model_ref == "mea") %>%
  mutate_at(vars(NAME_1, NAME_2), str_to_upper) %>% 
  pivot_wider(names_from = season, values_from = value) %>% 
  unite("group_id",NAME_1 :NAME_2) %>% 
  mutate(delta_spring_sm = 0,
         delta_summer_sm= 0)
#======================================================================
#subset CMIP6 scenarios and run bootsrap predictions
#======================================================================
#below 1.5 degree
cmip6_119 <-
  cmip6_tx_delta_change_data %>% 
  filter(model_scenario  == 119)

#predict function
predFun_ssp119 <- function(fit) {
  predict(fit,cmip6_119)
}

#bootstrap function
bc119 <-
  bootMer(
    model_fit_us,
    nsim = 1000,
    FUN = predFun_ssp119,
    seed = 101,
    use.u = TRUE,
    parallel = "multicore",
    ncpus = 10
  )

#store predictions
predictions_119 <- as.data.frame(t(bc119$t)) %>% 
  mutate(group_id = cmip6_119$group_id) 

#save file
saveRDS(predictions_119, file.path(dir_root,"prediction_soy_119.rds"))