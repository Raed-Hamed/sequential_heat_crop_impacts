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
#----------------------------------------------------------------------
#state spatial info EU
eu_state_shp <-
  list.files(dir_shp, "NUTS_RG_20M_2021_3035.shp", full.names = TRUE) %>%
  st_read() %>%
  filter(LEVL_CODE == 2) %>%
  st_transform(., crs = st_crs("+proj=longlat +datum=WGS84 +no_defs")) %>%
  filter(CNTR_CODE %in% spatial_info$State)
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
#load bootstrap yield predictions
yield_ci_boot <- map_dfr(list.files(dir_model,
                                    "prediction_" ,
                                    full.names = TRUE),
                         ~ {
                           data <- readRDS(.x)
                           data$id <-
                             basename(.x)  %>% file_path_sans_ext
                           data
                         }) %>%
  pivot_longer(-c(group_id,id)) %>%
  drop_na() %>% 
  group_by(group_id, id) %>%
  summarise(
    mean = mean(value),
    median = median(value),
    percentile_5th = quantile(value, 0.05),
    percentile_95th = quantile(value, 0.95)
  )  %>%
  ungroup() %>%
  separate(
    id,
    into = c("type", "crop", "ssp"),
    sep = "_",
    extra = "drop"
  ) %>%
  separate(group_id,
           into = c("State", "County"),
           sep = "_")
#-------------------------------------------------------------------------------
#compute bootstrap production CIs
production_ci_boot <-yield_ci_boot %>%
  inner_join(harvested_are_2010_2021) %>%
  pivot_longer(mean:percentile_95th,
               names_to = "CI_level",
               values_to = "yield_anomaly") %>% 
  mutate(production_anomaly = Value * yield_anomaly) %>%
  distinct() %>% 
  group_by(type,crop,ssp,CI_level) %>%
  summarise_at(
    vars(
      production_anomaly,
      Value
    ),
    sum
  )
#-------------------------------------------------------------------------------
#compute relative CI bootstrap predictions
relative_production_ci <-production_ci_boot %>% 
  mutate(production_prct_change = production_anomaly/Value*100)
#-------------------------------------------------------------------------------
#keep memory storage low
rm(yield_ci_boot)
gc()
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
#calculate soy yield projections
calc_yield_predictions_df <-
  inner_join(cmip6_tx_delta_change_per_model,
             boot_coefs) %>%
  dplyr::select(-NAME_0) %>% 
  group_by(model_scenario, model_ref, group_id,crop) %>%
  nest() %>%
  summarise(data = map(data, calc_yield_predictions)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct() 
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
      "Spring Temperature",
      "Summer Temperature",
      "Spring-Summer interaction",
      "Total"
    )
  )) %>% 
  filter(crop == "soy") 

#======================================================================
#load CMIP6 frequency changes from future minus historic period (to simplify)!
#======================================================================
#load tx frequency changes
cmip6_tx_freq_change <-
  fs::dir_ls(dir_cmip6_frequency_change,
             regexp = ".shp$") %>%
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
soy_harvest_area_freq_change <- cmip6_tx_freq_change %>%
  filter(percentile == 75) %>%
  mutate(State = NAME_1 %>%  str_to_upper()) %>%
  mutate(County = NAME_2 %>%  str_to_upper()) %>%
  inner_join(usda_soy_harvested_area) %>%
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
  mutate(crop = "Soybean")
#----------------------------------------------------------------------
maize_harvest_area_freq_change <-cmip6_tx_freq_change %>%
  filter(percentile == 75) %>%
  mutate(State = NAME_1 %>%  str_to_upper()) %>%
  mutate(County = NAME_2 %>%  str_to_upper()) %>%
  inner_join(usda_maize_harvested_area) %>%
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
    size = 4,
    position = position_dodge(width = 0.5)
  ) +
  ggtitle(label = "a) Change in hot-hot event frequency")+
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.8, 0.2),
    legend.background = element_rect(fill = alpha("white", 0)),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    legend.margin = margin(0),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 0, vjust = 0.5,size = 20),
    legend.spacing = unit(0, "cm")
  )+
  xlab("") +
  ylab("Frequency change in sequential spring-summer hot events (%)")
#----------------------------------------------------------------------
#plot impact projections
fi4b <-ggplot() +
  geom_bar(
    data = crop_production_anomaly %>%
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
    legend.position = c(0.25, 0.20),
    legend.background = element_rect(fill = alpha("white", 0)),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.title.y = element_text(size = 20),
    strip.text.x = element_text(size = 20),
    legend.margin = margin(0),
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(angle = 0, vjust = 0.5,size = 20),
    legend.spacing = unit(0, "cm")
  )+
  xlab("")+
  ylab("Percentage of recent total production (2010-2021)") +
  geom_hline(yintercept = 0, size = 1) 

#----------------------------------------------------------------------
#join the two panels and generate final figure
fig4 <-fi4a +fi4b
pdf("fig4.pdf",height=10,width=20)
print(fig4)
dev.off()


