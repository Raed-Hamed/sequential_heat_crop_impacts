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
dir_shp <- file.path(root_dir, "data/shapefiles")
dir_model <- file.path(root_dir, "data/model_data")
dir_cmip6 <- file.path(root_dir, "data/cmip6/")  
dir_figures <- file.path(root_dir, "figures")
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
#======================================================================
#load spatial info
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
#======================================================================
#load harvested data
#======================================================================
#load USDA harvested soy data
usda_soy_harvested_area <-
  list.files(dir_usda,
             "soybean_harvested_area",
             full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "SOYBEANS - ACRES HARVESTED") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>%
  filter(County != "OTHER COUNTIES") %>%
  mutate(Value = gsub(",", "", Value) %>%  as.numeric) %>% 
  mutate(Value = Value /  2.471) %>%
  arrange(Year) %>%
  dplyr::select(-Data.Item) %>%
  filter(Value != 0) %>%
  group_by(State, County) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(Value), mean) %>%
  ungroup()%>%
  mutate(normalizer  = 1 / sum(Value)) %>%
  mutate(weighted_harvest_area = normalizer * Value) %>%
  dplyr::select(-normalizer) %>% 
  mutate(crop = "soy")
#----------------------------------------------------------------------
#load USDA harvested maize data
usda_maize_harvested_area <-
  list.files(dir_usda,
             "corn_harvested_area",
             full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "CORN, GRAIN - ACRES HARVESTED") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>%
  filter(County != "OTHER COUNTIES") %>%
  mutate(Value = gsub(",", "", Value) %>%  as.numeric) %>%
  mutate(Value = Value /  2.471) %>%
  arrange(Year) %>%
  dplyr::select(-Data.Item) %>%
  filter(Value != 0) %>%
  group_by(State, County) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(Value), mean) %>%
  ungroup() %>% 
  mutate(normalizer  = 1 / sum(Value)) %>%
  mutate(weighted_harvest_area = normalizer * Value) %>%
  dplyr::select(-normalizer) %>% 
  mutate(crop = "corn")
#----------------------------------------------------------------------
#join harvesting areas
harvested_are_2010_2021 <-
  bind_rows(usda_soy_harvested_area, usda_maize_harvested_area)
#----------------------------------------------------------------------
#load average soy recent yield
avg_recent_soy_yield <- list.files(dir_model,
                                   "absolute_recent_soy_yield_us.rds" ,
                                   full.names = TRUE) %>%
  readRDS()
#----------------------------------------------------------------------
#load average corn recent yield
avg_recent_maize_yield <- list.files(dir_model,
                                     "absolute_recent_corn_yield_us.rds" ,
                                     full.names = TRUE) %>%
  readRDS()
#----------------------------------------------------------------------
#calculate recent production area (2010-2021)
crop_recent_production_avg <- avg_recent_soy_yield %>%
  mutate(crop = "soy") %>%
  bind_rows(avg_recent_maize_yield %>%
              mutate(crop = "corn")) %>%
  filter(absolute_yield  != 0) %>%
  group_by(State, County, crop) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(absolute_yield), mean) %>%
  ungroup() %>%
  inner_join(harvested_are_2010_2021) %>%
  mutate(production_total = absolute_yield * Value) %>%
  dplyr::select(-weighted_harvest_area)
#----------------------------------------------------------------------
#clean some memory space
rm(
  avg_recent_maize_yield,
  avg_recent_soy_yield,
  usda_maize_harvested_area,
  usda_soy_harvested_area
)
#======================================================================
#load crop yield projections
#======================================================================
#calculate soy yield projections
calc_yield_predictions_df <-
  list.files(dir_cmip6,
             "crop_yield_projections.rds" ,
             full.names = TRUE) %>%
  readRDS() 
#----------------------------------------------------------------------
#transform yield projections into production anomalies
crop_production_anomaly <- calc_yield_predictions_df %>%
  dplyr::select(
    model_scenario,
    group_id,
    model_ref,
    tx_spring_effect,
    '(Intercept)',
    tx_summer_effect,
    sequential_heat_effect,
    predicted_yield,
    crop
  ) %>%
  separate(group_id, c("State", "County"), sep = "_") %>%
  inner_join(harvested_are_2010_2021 %>%
               dplyr::select(-weighted_harvest_area)) %>%
  pivot_longer(tx_spring_effect:predicted_yield) %>%
  mutate(production_change = value * Value) %>%
  inner_join(crop_recent_production_avg) %>%
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
  filter(name != "(Intercept)")  %>%
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
      "Spring-Summer temperature interaction",
      "Total"
    )
  ))
#======================================================================
#load CMIP6 frequency changes from future minus historic period (to simplify)!
#======================================================================
#load tx frequency changes
cmip6_tx_freq_change <-
  fs::dir_ls(dir_cmip6,
             regexp = "_75.shp$") %>%
  map_dfr(read_sf, .id = "source") %>%
  tibble()  %>%
  filter(!NAME_1  %in% c("Alaska", "Hawaii")) %>% 
  mutate(source = str_trunc(source, 41, side = "left", ellipsis = "")) %>%
  separate(
    source,
    c(
      "model_scenario",
      "to_remove",
      "reference",
      "to_remove1",
      "percentile",
      "to_remove2"
    )
  ) %>%
  dplyr::select(-c(to_remove, to_remove1, to_remove2))
#----------------------------------------------------------------------
#average frequency changes over soybean harvested area using local weights
soy_harvest_area_freq_change <- cmip6_tx_freq_change %>%
  filter(percentile == 75) %>%
  mutate(State = NAME_1 %>%  str_to_upper()) %>%
  mutate(County = NAME_2 %>%  str_to_upper()) %>%
  inner_join(harvested_are_2010_2021 %>%
               filter(crop == "soy")) %>%
  dplyr::select(-geometry, -X.M_incre, -X.M_ampli) %>%
  pivot_longer(M_0:M_23) %>%
  dplyr::select(model_scenario,
                reference,
                State,
                County,
                Value,
                weighted_harvest_area,
                name,
                value) %>%
  drop_na() %>%
  group_by(model_scenario, reference, name) %>%
  nest() %>%
  mutate(weighted_val = map(data, ~ weighted.mean(.$value, .$weighted_harvest_area))) %>%
  dplyr::select(-data) %>%
  unnest(c(weighted_val)) %>%
  ungroup() %>%
  mutate(crop = "Soybean")
#----------------------------------------------------------------------
#average frequency changes over maize harvested area using local weights
maize_harvest_area_freq_change <-cmip6_tx_freq_change %>%
  filter(percentile == 75) %>%
  mutate(State = NAME_1 %>%  str_to_upper()) %>%
  mutate(County = NAME_2 %>%  str_to_upper()) %>%
  inner_join(harvested_are_2010_2021 %>%
               filter(crop == "corn")) %>%
  dplyr::select(-geometry,-X.M_incre,-X.M_ampli) %>%
  pivot_longer(M_0:M_23) %>%
  dplyr::select(model_scenario,
                reference,
                State,
                County,
                Value,
                weighted_harvest_area,
                name,
                value) %>%
  drop_na() %>% 
  group_by(model_scenario, reference, name) %>%
  nest() %>%
  mutate(weighted_val = map(data, ~ weighted.mean(.$value, .$weighted_harvest_area))) %>%
  dplyr::select(-data) %>%
  unnest(c(weighted_val)) %>%
  ungroup() %>%
  mutate(crop = "Maize")
#======================================================================
#Generate summary plot (projections + frequency changes)
#======================================================================
#plot frequency change
fi4a <-bind_rows(soy_harvest_area_freq_change,
                 maize_harvest_area_freq_change) %>%
  filter(reference  == "historic") %>% 
  filter(!name %in% c("mean_incre"))%>% 
  mutate(weighted_val =  weighted_val * 100) %>%
  group_by(model_scenario,reference,crop) %>% 
  mutate(lwr = quantile(weighted_val ,0.05),
         mean = mean(weighted_val ),
         upr = quantile(weighted_val ,0.95)) %>% 
  ungroup() %>% 
  mutate(
    model_scenario = case_when(
      model_scenario == "p119" ~ "ssp119",
      model_scenario == "p126" ~ "ssp126",
      model_scenario == "p245" ~ "ssp245",
      model_scenario == "p370" ~ "ssp370",
      TRUE ~ model_scenario
    )
  ) %>%
  mutate(model_scenario  = factor(
    model_scenario ,
    labels = c("SSP1 1.9",
               "SSP1 2.6",
               "SSP2 4.5",
               "SSP3 7.0")
  )) %>%
  mutate(reference =  paste(reference, "climatology")) %>%
  mutate(reference  = factor(
    reference ,
    levels = c("historic climatology",
               "future climatology")
  )) %>%
  ggplot() +
  geom_errorbar(
    aes(
      x = model_scenario,
      group = crop,
      ymin = lwr,
      ymax = upr
    ),
    width = 0.1,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    aes(
      x = model_scenario,
      y = mean,
      color = crop,
      group = crop
    ),
    size = 6,
    position = position_dodge(width = 0.5)
  ) +
  ggtitle(label = "a) Change in hot-hot event frequency")+
  theme_bw(base_size = 25) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.8, 0.2),
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
  xlab("") +
  ylab("Frequency change in sequential spring-summer hot events (%)")
#----------------------------------------------------------------------
#plot soybean impact projections
fi4b <-ggplot() +
  geom_bar(
    data = crop_production_anomaly %>%
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
    data = crop_production_anomaly %>%
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
  ggtitle(label = "b) Total production of soybean in the US")+
  scale_fill_manual("", values = c("aquamarine4", "tan3", "darkgoldenrod1", "darkred")) +
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
  geom_hline(yintercept = 0, size = 1) 
#----------------------------------------------------------------------
fig4c <-ggplot() +
  geom_bar(
    data = crop_production_anomaly %>%
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
    data = crop_production_anomaly %>% 
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
  ggtitle(label = "c) Total production of maize in the US")+
  scale_fill_manual("", values = c("aquamarine4", "tan3", "darkgoldenrod1", "darkred")) +
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
  geom_hline(yintercept = 0, size = 1) 
#----------------------------------------------------------------------
#join the two panels and generate final figure
fig4 <-fi4a + (fi4b/fig4c)

png(
  file.path(dir_figures, "tx_frequency_and_impact_change.png"),
  width = 20,
  height = 12,
  units = 'in',
  res = 300
)

print(fig4)
dev.off()


