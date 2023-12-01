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
dir_figures <- file.path(root_dir, "figures")
#======================================================================
#load helper function
#======================================================================
# Define the prediction function based on the complex model equation
calc_yield_pred_complex <- function(data) {
  coef_intercept <- data$`coef_(Intercept)`
  coef_tmx_spring <- data$coef_mean_tmax_spring
  coef_tmx_summer <- data$coef_mean_tmax_summer
  coef_sequential_heat <-data$`coef_mean_tmax_spring:mean_tmax_summer`
  
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
#---------------------------------------------------------------------
# Define the prediction function based on the simple model equation
calc_yield_pred_simple <- function(data) {
  
  coef_intercept <- data$`coef_(Intercept)`
  coef_tmx_spring <- data$coef_mean_tmax_spring
  coef_tmx_summer <- data$coef_mean_tmax_summer

    data$predicted_yield <- coef_intercept +
    coef_tmx_spring * data$delta_spring_tx  +
    coef_tmx_summer * data$delta_summer_tx 
  
  data$tx_spring_effect <- 
    coef_tmx_spring * data$delta_spring_tx
  
  data$tx_summer_effect <- 
    coef_tmx_summer * data$delta_summer_tx
  
  return(data)
}
#======================================================================
#Load recent yield and calculate weighted harvested areas
#======================================================================
#load USDA harvested soy data
global_harvested_area <-
  list.files(dir_model,
             "global_crop_harvested_area",
             full.names = TRUE) %>%
  readRDS()
#----------------------------------------------------------------------
#load average soy recent yield
global_absolute_yield <-
  list.files(dir_model,
             "global_absolute_recent_yield.rds" ,
             full.names = TRUE) %>%
  readRDS() %>% 
  group_by(State,County,crop) %>% 
  summarise_at(vars(absolute_yield),mean) %>% 
  ungroup()
#----------------------------------------------------------------------
#calculate recent production area (2010-2021)
recent_production_area <-
  inner_join(global_harvested_area, global_absolute_yield) %>%
  mutate(production_total = absolute_yield * Value) 
#======================================================================
#load model coefficients
#======================================================================
#load complex model incl. all predictors
model_complex_coefs <-
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
  separate(group_id, c("State", "County"), sep = "_") %>% 
  filter(crop %in% c("corn", "soy"))
#---------------------------------------------------------------------
#load simple model incl. all predictors
model_simple_coefs <-
  list.files(dir_model,
             "_no_sequential_interaction.rds" ,
             full.names = TRUE) %>%
  map( ~ {
    filename <- basename(.)
    crop <-
      str_extract(filename,
                  "(?<=model_)(.*?)(?=\\_no_sequential_interaction.rds)")
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
  rename_with( ~ paste0("coef_", .), -all_of(c("crop", "group_id"))) %>%
  separate(group_id, c("State", "County"), sep = "_") %>%
  filter(crop %in% c("corn", "soy"))
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
  dplyr::select(NAME_0, NAME_1, NAME_2, name, value) %>%
  separate(name, c("model_scenario", "season", "model_ref")) %>%
  mutate(across('season', str_replace, 'su', 'delta_summer_tx')) %>%
  mutate(across('season', str_replace, 'sp', 'delta_spring_tx')) %>%
  filter(!NAME_1  %in% c("Alaska", "Hawaii")) %>%
  mutate_at(vars(NAME_1, NAME_2), str_to_upper) %>%
  pivot_wider(names_from = season, values_from = value) %>%
  rename(State = NAME_1, County = NAME_2) %>%
  mutate(delta_spring_sm = 0,
         delta_summer_sm = 0)
#======================================================================
#calculate yield projections per temperature driver (tx predictor)
#======================================================================
#calculate crop yield projections
calc_yield_predictions_simple <-
  inner_join(cmip6_tx_delta_change_per_model,
             model_simple_coefs) %>%
  dplyr::select(-NAME_0) %>%
  group_by(model_scenario, model_ref, State,County, crop) %>%
  nest() %>%
  summarise(data = map(data, calc_yield_pred_simple)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct()
#----------------------------------------------------------------------
#calculate crop yield projections
calc_yield_predictions_complex <-
  inner_join(cmip6_tx_delta_change_per_model,
             model_complex_coefs) %>%
  dplyr::select(-NAME_0) %>%
  group_by(model_scenario, model_ref, State,County, crop) %>%
  nest() %>%
  summarise(data = map(data, calc_yield_pred_complex)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct()
#===============================================================================
#Join model results to harvested area to transform to production
#===============================================================================
production_model_simple <-calc_yield_predictions_simple %>% 
  inner_join(global_harvested_area %>%
             dplyr::select(-weighted_harvest_area)) %>%
  pivot_longer(predicted_yield:tx_summer_effect) %>%
  mutate(production_change = value * Value) %>%
  inner_join(recent_production_area) %>%
  group_by(model_scenario, name, crop, model_ref) %>%
  summarise_at(vars(production_change, production_total), sum) %>%
  mutate(production_anomaly = production_change / production_total  * 100) %>% 
mutate(name = factor(
  name,
  levels = c(
    "tx_spring_effect",
    "tx_summer_effect",
    "predicted_yield"
  )
)) %>%
  group_by(model_scenario, name, crop) %>%
  mutate(upr = max(production_anomaly)) %>%
  mutate(lwr = min(production_anomaly)) %>%
  ungroup() %>%
  mutate(model_scenario  = factor(
    model_scenario ,
    labels = c("SSP1 1.9",
               "SSP1 2.6",
               "SSP2 4.5",
               "SSP3 7.0")
  )) %>%
  mutate(name  = factor(
    name ,
    labels = c(
      "Spring temperature",
      "Summer temperature",
      "Total"
    )
  ))
  
#----------------------------------------------------------------------
production_model_complex <-calc_yield_predictions_complex %>% 
  inner_join(global_harvested_area %>%
               dplyr::select(-weighted_harvest_area)) %>%
  pivot_longer(predicted_yield:sequential_heat_effect) %>%
  mutate(production_change = value * Value) %>%
  inner_join(recent_production_area) %>%
  group_by(model_scenario, name, crop, model_ref) %>%
  summarise_at(vars(production_change, production_total), sum) %>%
  mutate(production_anomaly = production_change / production_total  * 100) %>% 
  mutate(name = factor(
    name,
    levels = c(
      "tx_spring_effect",
      "tx_summer_effect",
      "sequential_heat_effect",
      "predicted_yield"
    )
  )) %>%
  group_by(model_scenario, name,crop) %>%
  mutate(upr = max(production_anomaly)) %>%
  mutate(lwr = min(production_anomaly)) %>%
  ungroup() %>%
  mutate(model_scenario  = factor(
    model_scenario ,
    labels = c("SSP1 1.9",
               "SSP1 2.6",
               "SSP2 4.5",
               "SSP3 7.0")
  )) %>%
  mutate(name  = factor(
    name ,
    labels = c(
      "Spring temperature",
      "Summer temperature",
      "Spring-Summer temperature interaction",
      "Total"
    )
  ))
#======================================================================
#Plots 
#======================================================================
#plot soybean impact projections
fi4b <-ggplot() +
  geom_bar(
    data = production_model_simple %>%
      filter(crop == "soy") %>% 
      filter(model_ref == "mea"),
    aes(
      x = model_scenario,
      y = production_anomaly,
      fill = name,
      group = name
    ),
    stat = "identity",
    position = "dodge"
  ) +
  geom_errorbar(
    data = production_model_simple %>%
      filter(crop == "soy") %>% 
      filter(model_ref != "mea") %>% 
      group_by(model_scenario,name,crop) %>% 
      mutate(lwr = quantile(production_anomaly,0.05),
             upr = quantile(production_anomaly,0.95)),
    aes(x = model_scenario,
        group = name,
        ymin = lwr,
        ymax = upr),
    width = 0.1,
    position = position_dodge(width = 0.9)
  ) +
  geom_hline(yintercept = 0)+
  theme_bw() +
  ggtitle(label = "b) Total production change for soybean in the US")+
  scale_fill_manual("", values = c("aquamarine4", "tan3", "darkred")) +
  theme(
    panel.grid = element_blank(),
    legend.position="none",
    #legend.position = c(0.25, 0.20),
    legend.background = element_rect(fill = alpha("white", 0)),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 25),
    strip.text.x = element_text(size = 25),
    legend.margin = margin(0),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(angle = 0, vjust = 0.5,size = 25),
    legend.spacing = unit(0, "cm")
  )+
  xlab("")+
  #ylab("Percentage of recent total production (2010-2021)") +
  ylab("") +
  geom_hline(yintercept = 0, size = 1) +
  ylim(c(-20,7))
#----------------------------------------------------------------------
fig4c <-ggplot() +
  geom_bar(
    data = production_model_simple %>%
      filter(crop == "corn") %>%
      filter(model_ref == "mea"),
    aes(
      x = model_scenario,
      y = production_anomaly,
      fill = name,
      group = name
    ),
    stat = "identity",
    position = "dodge"
  ) +
  geom_errorbar(
    data = production_model_simple %>% 
      filter(crop == "corn") %>%
      filter(model_ref != "mea") %>% 
      group_by(model_scenario,name,crop) %>% 
      mutate(lwr = quantile(production_anomaly,0.05),
             upr = quantile(production_anomaly,0.95)),
    aes(x = model_scenario,
        group = name,
        ymin = lwr,
        ymax = upr),
    width = 0.1,
    position = position_dodge(width = 0.9)
  ) +
  geom_hline(yintercept = 0)+
  theme_bw() +
  ggtitle(label = "c) Total production change for maize in the US")+
  scale_fill_manual("", values = c("aquamarine4", "tan3", "darkred")) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.45, 0.20),
    legend.background = element_rect(fill = alpha("white", 0)),
    plot.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 25, hjust = -0.35),
    strip.text.x = element_text(size = 25),
    legend.margin = margin(0),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(angle = 0, vjust = 0.5,size = 25),
    legend.spacing = unit(0, "cm")
  )+
  xlab("")+
  ylab("Percentage of recent total production (2010-2021)") +
  geom_hline(yintercept = 0, size = 1) +
  ylim(c(-20,7))


#----------------------------------------------------------------------
#join the two panels and generate final figure
fig4 <-(fi4b/fig4c)

png(
  file.path(dir_figures, "impact_change_no_interaction.png"),
  width = 21,
  height = 11,
  units = 'in',
  res = 300
) 
print(fig4)
dev.off()
