#copy pasted to be adjusted
#======================================================================
#load model climate df
climate_df <-#list.files(dir_data,"", full.names = TRUE)%>% 
 #rastertopoints

us_climate <- us_climate_xy%>% 
  #do you want to filter this tibble based on climate zones? (e.g. below)
  #inner_join(cz_raster, by = c("x","y")) %>% 
  #filter(cz %in% c(25,14)) %>% 
  st_as_sf(x = ., coords = c("x", "y"),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") 
#======================================================================
#Aggregate grid based climate information to county scale
df_point_list <- split(dplyr::select(us_climate, -Year), 
                       us_climate$Year)

#Split sf usda_resid country by Years
df_poly_list <- split(usda_resid_county, usda_resid_county$Year)


#join yield and climate sf, na are reproduced as climate_df is filtered
full_Sf_yield_climate <- map2_dfr(df_poly_list, df_point_list,
                                  ~ .x %>% 
                                    st_join(.y, left =FALSE) %>% 
                                    group_by(County,Year,State) %>% 
                                    mutate_at(vars(-geometry,-State,-County,-Year),~(mean(.,na.rm = TRUE))))

#======================================================================
#construct final model data frame
model_us <- full_Sf_yield_climate %>%
  drop_na()%>%
  group_by(County,State) %>% 
  mutate_at(vars(-Year,-County, -geometry,-State), ~ scale(.,center = TRUE, scale = TRUE))%>%
  ungroup() %>% 
  inner_join(us_counties() %>%
               us_spatial_sf %>%
               as_tibble() %>%
               dplyr::select(-geometry),
             by = c("State","County")) %>% 
  drop_na() %>% 
  unique()