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
library(colorspace)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_eurostat <- file.path(root_dir, "data/eurostat")
dir_model <- file.path(root_dir, "data/model_data")
#======================================================================
#load helper functions
#======================================================================
source("00_load_functions.R")
#======================================================================
#load shape-files 
#======================================================================
#load EU shape file
shp_eu <- list.files(dir_eurostat, "shp", full.names = TRUE) %>%
  read_sf() %>%
  filter(LEVL_CODE %in% c(0, 1, 2, 3)) %>%
  dplyr::select(NUTS_ID, CNTR_CODE) %>%
  tibble() 
#----------------------------------------------------------------------
#load US shape file
shp_us <- list.files(dir_usda, "sf.rds", full.names = TRUE) %>%
  readRDS() %>%
  dplyr::select(NAME_1, NAME_2) %>%
  rename(State = NAME_1,
         County = NAME_2) %>%
  mutate_at(vars(State, County), ~ toupper(.))
#======================================================================
#load crop yields and de-trend
#======================================================================
#load USDA wheat yield data
usda_wheat_yield <-
  list.files(dir_usda, "usda-nass-wheat", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>%
  mutate(Value = Value / 14.87) %>%
  dplyr::select(-Data.Item) %>%
  arrange(Year) %>%
  filter(Value != 0) 
#----------------------------------------------------------------------
#de-trend USDA wheat yield data
usda_wheat_df <-usda_wheat_yield%>% 
  group_by(State,County) %>% 
  nest() %>% 
  mutate(model_trend = map(data,~lm(Value ~ Year, data = . ))) %>% 
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(State,County,Year,.resid) %>%
  rename(Value = .resid)
#----------------------------------------------------------------------
#load USDA corn yield data
usda_corn_yield <- list.files(dir_usda, "usda-nass-corn", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "CORN, GRAIN - YIELD, MEASURED IN BU / ACRE") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>% 
  filter(County != "OTHER COUNTIES") %>% 
  mutate(Value = Value / 14.87) %>%
  arrange(Year) %>% 
  dplyr::select(-Data.Item) %>% 
  filter(Value != 0)
#----------------------------------------------------------------------
#de-trend USDA corn yield data
usda_corn_df <- usda_corn_yield %>% 
group_by(State, County) %>%
  nest() %>%
  mutate(model_trend = map(data,  ~ lm(Value ~ Year, data = .))) %>%
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(State, County, Year, .resid) %>%
  rename(Value = .resid)
#----------------------------------------------------------------------
#load USDA soy yield data
usda_soy_yield <-
  list.files(dir_usda, "usda-nass-soy", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "SOYBEANS - YIELD, MEASURED IN BU / ACRE") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>%
  filter(County != "OTHER COUNTIES") %>%
  filter(County != "") %>%
  mutate(Value = Value / 14.87) %>%
  arrange(Year) %>%
  dplyr::select(-Data.Item) %>%
  filter(Value != 0) %>%
  distinct()
#----------------------------------------------------------------------
#de-trend USDA soy yield data
usda_soy_df <- usda_soy_yield %>%
  group_by(State, County) %>%
  nest() %>%
  mutate(model_trend = map(data,  ~ lm(Value ~ Year, data = .))) %>%
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(State, County, Year, .resid) %>%
  rename(Value = .resid)
#----------------------------------------------------------------------
#load euro-stat wheat yield data
eurostat_wheat_yield <-
  list.files(dir_eurostat, "cropstats", full.names = TRUE) %>%
  read.csv() %>%
  filter(TYPE == "Yield") %>%
  filter(CROP_NAME == "C1100") %>%
  dplyr::select(IDREGION, YEAR, VALUE) %>%
  distinct() %>%
  rename(NUTS_ID = IDREGION) %>%
  arrange(YEAR) %>%
  filter(VALUE != 0)
#----------------------------------------------------------------------
eurostat_wheat_df <-  eurostat_wheat_yield %>%
  group_by(NUTS_ID) %>%
  nest() %>%
  mutate(model_trend = map(data,  ~ lm(VALUE ~ YEAR, data = .))) %>%
  mutate(augment_trend  = map(model_trend,  ~ broom::augment(.))) %>%
  unnest(cols = c(augment_trend)) %>%
  ungroup() %>%
  dplyr::select(NUTS_ID, YEAR, .resid) %>%
  rename(Value = .resid)
#======================================================================
#store absolute yield estimates for period 2010-2021
#======================================================================
#maize absolute yield average 2010-2021
absolute_recent_corn_yield <-usda_corn_yield %>%
  rename(absolute_yield = Value) %>%
  mutate(crop = "corn") %>% 
  filter(Year %in% 2010:2021) 
#----------------------------------------------------------------------
#soy absolute yield average 2010-2021
absolute_recent_soy_yield <-usda_soy_yield %>%
  rename(absolute_yield = Value) %>%
  mutate(crop = "soy") %>% 
  filter(Year %in% 2010:2021) 
#----------------------------------------------------------------------
#wheat us absolute yield average 2010-2021
absolute_recent_wheat_yield_us <-usda_wheat_yield %>%
  rename(absolute_yield = Value) %>%
  mutate(crop = "wheat_us") %>% 
  filter(Year %in% 2010:2021) 
#----------------------------------------------------------------------
#wheat eu absolute yield average 2010-2021
absolute_recent_wheat_yield_eu <- eurostat_wheat_yield %>%
  inner_join(shp_eu %>% dplyr::select(-geometry)) %>%
  rename(
    County = 1,
    Year = 2,
    State = 4,
    Value = 3
  ) %>% 
  rename(absolute_yield = Value) %>%
  mutate(crop = "wheat_eu") %>% 
  filter(Year %in% 2010:2021) 
#----------------------------------------------------------------------
#join both into one spatial informed data-frame
absolute_recent_yield_all_crops <- absolute_recent_corn_yield %>%
  bind_rows(
    absolute_recent_soy_yield,
    absolute_recent_wheat_yield_us,
    absolute_recent_wheat_yield_eu
  ) 
#======================================================================
#Load recent yield and calculate weighted harvested areas
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
  group_by(State, County) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(Value), mean) %>%
  ungroup()%>%
  mutate(normalizer  = 1 / sum(Value)) %>%
  mutate(weighted_harvest_area = normalizer * Value) %>%
  dplyr::select(-normalizer) %>% 
  mutate(crop = "corn")
#----------------------------------------------------------------------
#load USDA harvested maize data
usda_wheat_harvested_area <-
  list.files(dir_usda,
             "wheat_harvested_area",
             full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "WHEAT, WINTER - ACRES HARVESTED") %>%
  filter(County != "OTHER (COMBINED) COUNTIES") %>%
  filter(County != "OTHER COUNTIES") %>%
  mutate(Value = gsub(",", "", Value) %>%  as.numeric) %>% 
  group_by(State, County) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(Value), mean) %>%
  ungroup()%>%
  mutate(normalizer  = 1 / sum(Value)) %>%
  mutate(weighted_harvest_area = normalizer * Value) %>%
  dplyr::select(-normalizer) %>% 
  mutate(crop = "wheat_us")
#----------------------------------------------------------------------
absolute_recent_wheat_yield_eu <- eurostat_wheat_yield %>%
  inner_join(shp_eu %>% dplyr::select(-geometry)) %>%
  rename(
    County = 1,
    Year = 2,
    State = 4,
    Value = 3
  ) %>% 
  rename(absolute_yield = Value) %>%
  mutate(crop = "wheat_eu") %>% 
  filter(Year %in% 2010:2021) 

#load USDA harvested maize data
eurostat_wheat_harvested_area <-
  list.files(dir_eurostat,
             "cropstats",
             full.names = TRUE) %>%
  read.csv() %>%
  filter(TYPE == "Area") %>%
  filter(CROP_NAME == "C1100") %>%
  dplyr::select(IDREGION, YEAR, VALUE) %>%
  distinct() %>%
  rename(NUTS_ID = IDREGION) %>%
  arrange(YEAR) %>%
  group_by(NUTS_ID) %>%
  drop_na() %>% 
  filter(YEAR %in% 2010:2021) %>%
  summarise_at(vars(VALUE), mean) %>%
  ungroup()%>%
  mutate(normalizer  = 1 / sum(VALUE)) %>%
  mutate(weighted_harvest_area = normalizer * VALUE) %>%
  dplyr::select(-normalizer) %>% 
  mutate(crop = "wheat_eu") %>% 
  inner_join(shp_eu %>% dplyr::select(-geometry)) %>% 
  rename(County = NUTS_ID) %>% 
  rename(Value = VALUE) %>% 
  rename(State = CNTR_CODE)
  #----------------------------------------------------------------------
#join both into one spatial informed data-frame
global_crop_harvested_area <- usda_soy_harvested_area %>%
  bind_rows(
    usda_maize_harvested_area,
    usda_wheat_harvested_area,
    eurostat_wheat_harvested_area
  ) 
#======================================================================
#transform data frame into SF object (i.e. add spatial dimension)
#======================================================================
#add spatial dimension to us wheat data-frame
usda_wheat_sf <- shp_us %>%
  left_join(usda_wheat_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>% 
  mutate(zone = "US") %>% 
  tibble()
#----------------------------------------------------------------------
#add spatial dimension to us corn dataframe
usda_corn_sf <- shp_us %>%
  left_join(usda_corn_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>% 
  mutate(zone = "US")%>% 
  tibble()
#----------------------------------------------------------------------
#add spatial dimension to us soy dataframe
usda_soy_sf <-shp_us %>%
  left_join(usda_soy_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>% 
  mutate(zone = "US")%>% 
  tibble()
#----------------------------------------------------------------------
#add spatial dimension to eu dataframe
eurostat_wheat_sf <- eurostat_wheat_df %>%
  inner_join(shp_eu) %>%
  group_by(NUTS_ID) %>%
  mutate(n = n()) %>%
  rename(
    County = 1,
    Year = 2,
    State = 4,
    Value = 3
  ) %>%
  ungroup() %>% 
  mutate(zone = "EU")%>% 
  tibble()
#----------------------------------------------------------------------
#join both into one spatial informed dataframe
global_wheat_sf <- eurostat_wheat_sf %>%
  bind_rows(usda_wheat_sf) 
#======================================================================
#Diagnostic plot: How many years do we have per county
#======================================================================
usda_corn_sf %>%
  dplyr::select(n, geometry, zone) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  facet_grid(~zone)+
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")
#----------------------------------------------------------------------
usda_soy_sf %>%
  dplyr::select(n, geometry, zone) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  facet_grid(~zone)+
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")
#----------------------------------------------------------------------
global_wheat_sf %>%
  dplyr::select(n, geometry, zone) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  facet_grid(~zone)+
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")
#======================================================================
#save crop yield data
#======================================================================
#save spatial de-trended crop data sets
saveRDS(global_wheat_sf,
        file.path(dir_model, "global_wheat_dtr_sf.rds"))
saveRDS(usda_corn_sf, file.path(dir_model, "us_corn_dtr_sf.rds"))
saveRDS(usda_soy_sf, file.path(dir_model, "us_soy_dtr_sf.rds"))
#----------------------------------------------------------------------
#save absolute yield and harvested area global data sets
saveRDS(
  global_crop_harvested_area,
  file.path(dir_model, "global_crop_harvested_area.rds")
)
saveRDS(
  absolute_recent_yield_all_crops,
  file.path(dir_model, "global_absolute_recent_yield.rds")
)
