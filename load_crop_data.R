#======================================================================
library(tidyverse)
library(sf)
library(USAboundaries)
library(colorspace)
#======================================================================
#clean workspace
rm(list = ls())
graphics.off()
gc()
#======================================================================
conversation_metric <-14.87
#load functions
us_spatial_sf<-  function(x) {
  #x is either us_counties or us_states from usaboundaries package  
  x %>% 
    # filter(!state_name  %in% 
    #          c("Alaska","Puerto Rico","Hawaii","Washington",
    #            "Oregon","California","Nevada","Idaho","Utah",
    #            "Arizona","New Mexico","Colorado","Wyoming","Montana"))%>% 
    rename(County = name,
           State  = state_name) %>% 
    dplyr::select(County,State,geometry) %>% 
    mutate_at(vars(State,County), ~toupper(.))
  
}
#======================================================================
root_dir <-"C:/Users/rhd630/Desktop/PhD/Academic/paper_4"
root_data <-file.path(root_dir,"data/usda-nass")
#======================================================================
usda_wheat_df <-list.files(root_data, "wheat", full.names = TRUE) %>%
  map_dfr(read.csv) %>% 
  dplyr::select(Year,State,County,Data.Item,Value) %>% 
  filter(Data.Item == "WHEAT, WINTER - YIELD, MEASURED IN BU / ACRE") 
#======================================================================
usda_wheat_sf <- us_counties() %>%  us_spatial_sf %>%
  left_join(usda_wheat_df, by = c("County","State"))%>% 
  drop_na() %>% 
  arrange(Year) 
#======================================================================
usda_shp <-usda_wheat_sf %>% 
  dplyr::select(geometry) %>% 
  tibble() %>% 
  distinct() %>% 
  st_as_sf()
#======================================================================
usda_wheat_sf %>% 
  dplyr::select(Year,Value) %>% 
  filter(Year %in% 2010:2020) %>% 
  ggplot()+
  geom_sf(data = usda_shp, color = "lightgray", alpha = 0.5)+
  geom_sf(aes(fill = Value))+
  facet_wrap(~Year)+
  scale_fill_continuous_sequential("Greens 2")+
  theme_bw()

