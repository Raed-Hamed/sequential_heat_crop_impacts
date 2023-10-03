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
library(broom.mixed)
source("00_load_functions.R")
library(piecewiseSEM)
library(patchwork)
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
dir_data <- file.path(dir_root, "data")
dir_model <- file.path(dir_data, "model_data")
dir_shp <- file.path(dir_data, "shapefiles")
dir_figures <- file.path(dir_root, "figures")
#======================================================================
#load helper function
#======================================================================
#Define the prediction function based on the model equation
calc_yield_predictions <- function(data) {

  coef_intercept <- data$`coef_(Intercept)`
  coef_tmx_spring <- data$coef_mean_tmax_spring
  coef_tmx_summer <- data$coef_mean_tmax_summer
  coef_tmx_interaction <- data$`coef_mean_tmax_spring:mean_tmax_summer`
  coef_sm_spring <- data$coef_mean_sm_spring
  coef_sm_summer <- data$coef_mean_sm_summer
  coef_sm_interaction <- data$`coef_mean_sm_spring:mean_sm_summer`
  

  data$predicted_yield <- coef_intercept +
    coef_tmx_spring * data$mean_tmax_spring  +
    coef_tmx_summer * data$mean_tmax_summer +
    coef_tmx_interaction * data$mean_tmax_spring * data$mean_tmax_summer +
    coef_sm_spring * data$mean_sm_spring  +
    coef_sm_summer * data$mean_sm_summer +
    coef_sm_interaction * data$mean_sm_spring * data$mean_sm_summer 

  return(data)
}
#----------------------------------------------------------------------
#Plot R-square
plot_rsq <- function(crop_select) {
  plot <- county_scale_rsq %>%
    ungroup() %>%
    mutate(crop = factor(
      crop,
      levels = c("corn",
                 "soy",
                 "wheat_us",
                 "wheat_eu")
      ,
      labels = c("Maize US",
                 "Soybean US",
                 "Wheat US",
                 "Wheat EU")
    )) %>%
    filter(crop == crop_select) %>%
    st_as_sf() %>%
    ggplot(aes(fill = rsq)) +
    geom_sf() +
    theme_void(base_size = 5) +
    scale_fill_continuous_sequential(name = "", "reds 3", limit = c(0,1)) +
    facet_grid(~crop)+
    theme(
      legend.position = "bottom",
      legend.key.size = unit(1, 'cm'),
      legend.title.align = 0.5,
      legend.title = element_text(size = 15, angle = 0),
      legend.key.width = unit(5, "cm"),
      legend.box = "horizontal",
      plot.title = element_text(size = 20, hjust = 0.5),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      legend.text = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
      
    ) +
    guides(fill = guide_colorbar(title.position = "bottom"))
  
  return(plot)
}

#======================================================================
#Load spatial data
#======================================================================
#spatial info from crop data
spatial_info <-
  list.files(dir_model,
             "spatial_county_group_id_ref.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  tibble() %>%
  dplyr::select(-group_id) %>%
  distinct() %>% 
  unite(crop, c(crop, zone), sep = "_") %>%
  mutate(
    crop = case_when(
      crop == "corn_US" ~ "corn",
      crop == "soy_US" ~ "soy",
      crop == "wheat_US" ~ "wheat_us",
      crop == "wheat_EU" ~ "wheat_eu",
      TRUE ~ crop
    )
  ) 
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
#Get quantile equivalents for variables of interest
#======================================================================
#load model data 
model_data <-
  list.files(dir_model, "model_data_", full.names = TRUE) %>%
  map_dfr( ~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_data_)(.*?)(?=\\.rds)")
    data <- readRDS(.)
    tibble(crop = crop, data)
  }) %>%
  dplyr::select(-geometry, -group_id) %>%
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
  )  %>% 
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
  unite("group_id",State:County) %>% 
  unite(crop, c(crop, zone), sep = "_") %>%
  mutate(
    crop = case_when(
      crop == "corn_US" ~ "corn",
      crop == "soy_US" ~ "soy",
      crop == "wheat_US" ~ "wheat_us",
      crop == "wheat_EU" ~ "wheat_eu",
      TRUE ~ crop
    )
  )
#----------------------------------------------------------------------
#store spring tx quantile values
spring_tx_quantile_values <-
  model_data %>%
  group_by(crop) %>%
  summarise_at(vars(c(mean_tmax_spring)),
               ~ data.frame(q50 = quantile(., 0.50),
                            q95 = quantile(., 0.95))) %>% 
  mutate(q50 = mean_tmax_spring$q50,
         q95 = mean_tmax_spring$q95) %>%
  dplyr::select(-mean_tmax_spring)
#======================================================================
#load model fit and extract key features
#======================================================================
#load model incl. all predictors
model_fit_full <-
  list.files(dir_model,
             "_full" ,
             full.names = TRUE) %>%
  map(~ {
    filename <- basename(.)
    crop <- str_extract(filename, "(?<=model_)(.*?)(?=\\_full.rds)")
    data <- readRDS(.)
    tibble(crop = crop, model_fit = list(data))
  }) %>%
  bind_rows() 
#----------------------------------------------------------------------
#get model r-squared values for full model
model_fit_full %>% 
  mutate(map_dfr(model_fit, ~ rsquared(.)))
#----------------------------------------------------------------------
#load model incl. all predictors minus sequential interaction terms (to compare)
model_fit_simple <-
  list.files(dir_model,
             "_interaction" ,
             full.names = TRUE) %>%
  map( ~ {
    filename <- basename(.)
    crop <-
      str_extract(filename,
                  "(?<=model_)(.*?)(?=\\_no_sequential_interaction.rds)")
    data <- readRDS(.)
    tibble(crop = crop, model_fit = list(data))
  }) %>%
  bind_rows() 
#----------------------------------------------------------------------
#get model r-squared values for simple model (to compare)
model_fit_simple %>% 
  mutate(map_dfr(model_fit, ~ rsquared(.)))
#----------------------------------------------------------------------
#store tx coefficients of interest
model_tx_coefs <- model_fit_full %>% 
  mutate(summer_coef = map(
    model_fit,
    ~ (.) %>%
      tidy() %>%
      filter(term == "mean_tmax_summer") %>%
      pull(estimate)
  ))%>% 
  mutate(interaction_coef = map(
    model_fit,
    ~ (.) %>%
      tidy() %>%
      filter(term == "mean_tmax_spring:mean_tmax_summer") %>%
      pull(estimate)
  )) %>% 
  dplyr::select(-model_fit) %>% 
  unnest(c(summer_coef,interaction_coef))
#----------------------------------------------------------------------
#store model coefficients 
model_coefs <- model_fit_full %>%
  mutate(model_coefs = map(
    model_fit,
    ~ coef(.) %>%
      .$group_id %>%
      rownames_to_column("group_id")
    
  )) %>%
  dplyr::select(-model_fit) %>%
  unnest(c(model_coefs)) %>%
  rename_with(~ paste0("coef_", .),-all_of(c("crop", "group_id")))
#======================================================================
#join model coefficients and quantile spring TX values
#======================================================================
calc_tx_spring_effect_on_summer_tx_impact <-
  inner_join(spring_tx_quantile_values, model_tx_coefs) %>%
  mutate(
    summer_effect_hot_spring = summer_coef + q95 * interaction_coef,
    summer_effect_mild_spring = summer_coef + q50 * interaction_coef
  ) %>%
  mutate(additional_hot_spring_impacts =
           summer_effect_hot_spring / summer_effect_mild_spring)
#======================================================================
#calculate county scale r-squared
#======================================================================
county_scale_rsq <-inner_join(model_data,
           model_coefs) %>% 
  group_by(crop,group_id) %>% 
  nest() %>%
  summarise(data = map(data, calc_yield_predictions)) %>%
  unnest(data)  %>%
  ungroup() %>%
  distinct() %>% 
  group_by(crop,group_id) %>% 
  summarise(rsq = cor(crop_yield,predicted_yield)^2) %>% 
  ungroup() %>%  
  distinct() %>% 
  separate(group_id, c("State","County"), sep = "_") %>% 
  inner_join(spatial_info)
#----------------------------------------------------------------------
#quick plotting
Maize_us <-plot_rsq(crop_select = "Maize US") +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125,-75), ylim = c(25, 50)) +
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

Soy_us <-plot_rsq(crop_select = "Soybean US") +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125,-75), ylim = c(25, 50)) +
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

Wheat_us <-plot_rsq(crop_select = "Wheat US") +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125,-75), ylim = c(25, 50)) +
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")

Wheat_eu <-plot_rsq(crop_select = "Wheat EU") +
  geom_sf(
    data = eu_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-10, 50), ylim = c(30, 70)) 
#----------------------------------------------------------------------
#Bring together in one plot
rsq_diagnostic_plot <-  Maize_us + Soy_us + Wheat_us + Wheat_eu+
plot_layout(guides = "collect",
            ncol = 2) &
  theme(legend.position = 'bottom')
#----------------------------------------------------------------------
#save R-square diagnostic plot
png(
  file.path(dir_figures, "rsq_diagnostic_plot.png"),
  width = 15,
  height = 15,
  units = 'in',
  res = 300
)
print(rsq_diagnostic_plot)
dev.off()
