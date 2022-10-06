#======================================================================
#clean workspace
#======================================================================
rm(list = ls())
graphics.off()
gc()
#======================================================================
#load libraries
#======================================================================
library(tidyverse)
library(sf)
library(USAboundaries)
library(colorspace)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_crop <- file.path(root_dir, "data/usda-nass")
source("load_functions.R")
#======================================================================
#load crop yields USDA dataframe
#======================================================================
usda_wheat_df <- list.files(dir_crop, "usda-nass-wheat", full.names = TRUE) %>%
  map_dfr(read.csv) %>%
  dplyr::select(Year, State, County, Data.Item, Value) %>%
  filter(Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE") %>%
  mutate(Value = Value / 14.87) %>%
  dplyr::select(-Data.Item)
#======================================================================
#transform USDA dataframe into SF object (i.e. add spatial dimension)
#======================================================================
usda_wheat_sf <- us_counties() %>%  us_spatial_sf %>%
  left_join(usda_wheat_df, by = c("County", "State")) %>%
  drop_na() %>%
  arrange(Year) %>%
  tibble() %>%
  group_by(State, County) %>%
  mutate(n = n()) %>%
  filter(n < 44) %>%
  ungroup() %>%
  st_as_sf()
#======================================================================
#Diagnostic plot: How many years do we have per county
#======================================================================
usda_wheat_sf %>%
  dplyr::select(n, geometry) %>%
  tibble() %>%
  distinct() %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(aes(fill = n)) +
  theme_bw(base_size = 20) +
  scale_fill_continuous_sequential("Viridis")
#======================================================================
#save crop yield data
#======================================================================
saveRDS(usda_wheat_sf, file.path(dir_crop,"usda_wheat_sf.rds"))
