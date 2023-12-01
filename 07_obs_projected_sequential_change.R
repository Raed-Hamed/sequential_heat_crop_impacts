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
source("00_load_functions.R")
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_shp <- file.path(root_dir, "data/shapefiles")
dir_model <- file.path(root_dir, "data/model_data")
dir_cmip6 <- file.path(root_dir, "data/cmip6/")
#======================================================================
#load spatial data 
#======================================================================
#spatial info from crop data for plotting
spatial_info <-
  list.files(dir_model,
             "spatial_county_group_id_ref.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  tibble()
#----------------------------------------------------------------------
#define Americas mask (Adjust)
sf_mask_country <- world %>%
  filter(region_un  == "Americas")
#----------------------------------------------------------------------
#us states for plotting
us_state_shp <- list.files(dir_shp, ".json", full.names = TRUE) %>%
  map_dfr(st_read) %>%
  mutate(NAME_1  = NAME_1 %>%  str_to_upper) %>%
  filter(!NAME_1 %in%  c("ALASKA", "HAWAII"))
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
# #load average soy recent yield
global_absolute_yield <-
  list.files(dir_model,
             "global_absolute_recent_yield.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  group_by(State,County,crop) %>%
  summarise_at(vars(absolute_yield),mean) %>%
  ungroup()
# #----------------------------------------------------------------------
# #load average corn recent yield
# recent_avg_yield_all <-
#   inner_join(global_harvested_area, global_absolute_yield) %>%
#   group_by(crop) %>%
#   nest() %>%
#   mutate(yield_est = map(
#     data,
#     ~ weighted.mean(.$absolute_yield, .$weighted_harvest_area)
#   )) %>%
#   unnest(yield_est) %>%
#   dplyr::select(-data) %>%
#   mutate(crop = case_when(crop == "soy" ~ "soy_us",
#                           crop == "corn" ~ "corn_us",
#                           TRUE ~ crop)) %>% 
#   rename(area = crop)
#----------------------------------------------------------------------
#sum soy and maize crop harvested area (give weights to grid-cells)
sum_crop_recent_harvested_area <- global_harvested_area %>%
  filter(crop %in% c("corn","soy")) %>% 
  pivot_wider(names_from = crop, values_from = Value) %>%
  mutate_at(vars(corn, soy), replace_na, 0) %>%
  mutate(total_harvest_area_pr_county = corn+soy)
#----------------------------------------------------------------------
#combined crop production area
crop_recent_production_avg <- global_absolute_yield %>%
  inner_join(global_harvested_area) %>%
  filter(crop %in% c("corn", "soy")) %>%
  mutate(production_total = absolute_yield * Value) %>%
  group_by(crop) %>%
  summarise_at(vars(production_total), sum)
#======================================================================
#Calculate temperature trends (singular and sequential) in observed data
#======================================================================
#load model df with observed CPC and GLEAM datasets
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
  group_by(zone, crop) %>%
  nest() %>%
  mutate(data = map(
    data,
    ~ group_by(., State, County) %>%
      mutate(group_id = cur_group_id() %>%  as.factor())
  )) %>%
  unnest(c(data)) 
#----------------------------------------------------------------------
#calculate sequential heat frequency in observed data
sequential_heat_freq <-model_df %>%
  filter(zone == "US") %>% 
  filter(crop != "wheat") %>% 
  ungroup() %>% 
  dplyr::select(State,
                County,
                Year,
                mean_tmax_spring,
                mean_tmax_summer) %>%
  distinct() %>% 
  group_by(State, County) %>%
  mutate(tmax_spring_75 = mean_tmax_spring %>%
           quantile(0.75) %>%
           round(1)) %>%
  mutate(tmax_summer_75 = mean_tmax_summer %>%
           quantile(0.75) %>%
           round(1)) %>% 
  mutate(warm_spring = ifelse(mean_tmax_spring >= tmax_spring_75,1,0)) %>% 
  mutate(warm_summer = ifelse(mean_tmax_summer >= tmax_summer_75,1,0)) %>% 
  mutate(sequential_heat = ifelse(warm_spring & warm_summer == 1,1,0)) %>% 
  inner_join(sum_crop_recent_harvested_area) %>% 
  group_by(Year) %>% 
  mutate(total_harvest = sum(total_harvest_area_pr_county, na.rm = TRUE)) 
#----------------------------------------------------------------------
#plot cropland area under sequential heat events
sequential_heat_freq %>%
  mutate(spring_heat_area = sum(ifelse(
    warm_spring == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  
  mutate(summer_heat_area = sum(ifelse(
    warm_summer == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  mutate(sequential_heat_area = sum(ifelse(
    sequential_heat == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  distinct() %>%
  pivot_longer(c(spring_heat_area, summer_heat_area, sequential_heat_area)) %>%
  ggplot(aes(x = Year, y = value, color = name)) + geom_line(size = 1) +
  facet_wrap(~ name)+
  ylab("% cropland area under sequential heat (spring & summer > 75th percentile")
#----------------------------------------------------------------------
#area weighted spring and summer absolute maximum temperature values
weighted_avg_model_data <-model_df %>%
  filter(zone == "US") %>% 
  filter(crop != "wheat") %>% 
  ungroup() %>% 
  dplyr::select(State,
                County,
                Year,
                mean_tmax_spring,
                mean_tmax_summer) %>%
  distinct() %>% 
  inner_join(sum_crop_recent_harvested_area) %>% 
  group_by(Year) %>% 
  mutate(normalizer  = 1 / sum(total_harvest_area_pr_county)) %>%
  mutate(weighted_harvest_area = normalizer * total_harvest_area_pr_county) %>%
  dplyr::select(-normalizer) %>% 
  pivot_longer(mean_tmax_spring:mean_tmax_summer) %>% 
  group_by(Year,name) %>%
  nest() %>%
  mutate(weighted_val = map(data, ~ weighted.mean(.$value, .$weighted_harvest_area))) %>%
  dplyr::select(-data) %>%
  unnest(c(weighted_val)) %>%
  ungroup() %>% 
  group_by(name) %>%
  mutate(percentile_75 = weighted_val %>% quantile(0.75) %>%
           round(1)) 
#----------------------------------------------------------------------
#what are years in observed dataset where we witness sequential heat events ?
weighted_avg_model_data %>%  
  mutate(heat_event = ifelse(weighted_val > percentile_75,1,0)) %>% 
  ungroup() %>% 
  dplyr::select(Year,name,heat_event) %>% 
  pivot_wider(names_from = name, values_from = heat_event) %>% 
  mutate(heat_year = ifelse(mean_tmax_spring & mean_tmax_summer ==1,1,0)) %>% 
  filter(heat_year == 1) %>% 
  .$Year %>% 
  unique()
#----------------------------------------------------------------------
#summarise this in a plot
weighted_avg_model_data %>%
  ggplot(aes(x = Year, y = weighted_val, color = name)) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = percentile_75, color = name), linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(1988, 1991, 2006 , 2010 , 2012)) +
  theme_bw()+
  scale_color_manual("",
                     values = c("orange", "firebrick3")) +
  theme(
    legend.position = c(0.5, 0.5),
    legend.background = element_rect(fill = alpha("white", 0.3),
                                     color = "transparent"),
    legend.title = element_text(size=10),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 10))+
  ylab("°C")

#======================================================================
#load CMIP6 data and do similar calculations over projections
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
  mutate(model_scenario  = factor(
    model_scenario ,
    labels = c("SSP1 1.9",
               "SSP1 2.6",
               "SSP2 4.5",
               "SSP3 7.0")
  )) %>% 
  filter(model_ref == "mea") %>% 
  pivot_longer(delta_summer_tx:delta_spring_tx)
#======================================================================
#weighted temperature change per season and model scenario
cmip6_tx_delta_change_per_model %>% 
  separate(group_id, c("State", "County"), sep = "_") %>%
  inner_join(sum_crop_recent_harvested_area %>% dplyr::select(-soy,-corn)) %>% 
  mutate(normalizer  = 1 / sum(total_harvest_area_pr_county)) %>%
  mutate(weighted_harvest_area = normalizer * total_harvest_area_pr_county) %>%
  dplyr::select(-normalizer) %>% 
  group_by(model_scenario,name) %>%
  nest() %>%
  mutate(weighted_val = map(data, ~ weighted.mean(.$value, .$weighted_harvest_area))) %>%
  dplyr::select(-data) %>%
  unnest(c(weighted_val)) %>%
  ungroup()
#======================================================================
plot_temperature_change <-cmip6_tx_delta_change_per_model %>% 
  mutate(name  = factor(
    name ,
    labels = c("Spring T change",
               "Summer T change"))) %>%
  separate(group_id, c("State", "County"), sep = "_") %>%
  inner_join(spatial_info %>% dplyr::select(-zone,-crop,-group_id) %>% distinct) %>% 
      st_as_sf() %>%
      ggplot(aes(fill = value)) +
      geom_sf(color = "transparent") +
      geom_sf(
        data = us_state_shp,
        color = "black",
        size = 0.00000001,
        fill = "transparent"
      ) +
      theme_void(base_size = 5) +
      scale_fill_continuous_sequential(name = "(°C)", "Reds 3") +
      facet_grid(name ~ model_scenario, switch = "y")+
      theme(
        legend.position = "bottom",
        legend.key.size = unit(1, 'cm'),
        legend.title.align = 0.5,
        legend.title = element_text(size = 15, angle = 0),
        legend.key.width = unit(5, "cm"),
        legend.box = "horizontal",
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15, angle = 90),
        legend.text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm")
      ) +
      guides(fill = guide_colorbar(title.position = "bottom"))+
      coord_sf(xlim = c(-125,-75), ylim = c(25, 50)) +
      coord_sf(expand = FALSE)+
      coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

#----------------------------------------------------------------------
#save wheat eu climate sensitivity plot
png(
  file.path(dir_figures, "delta_temperatures.png"),
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)
print(plot_temperature_change)
dev.off()

