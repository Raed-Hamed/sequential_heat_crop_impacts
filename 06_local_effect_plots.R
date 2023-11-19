#======================================================================
#clean work space
#======================================================================
rm(list = ls())
graphics.off()
gc()
#======================================================================
#load libraries
#======================================================================
library(sf)
library(tidyverse)
library(spData)
library(raster)
library(lubridate)
library(colorspace)
library(lme4)
library(piecewiseSEM)
library(lmerTest)
library(broom.mixed)
library(patchwork)
source("00_load_functions.R")
#======================================================================
#load helper functions
#======================================================================
load_summary_boot_coef_data <- function(filename) {
  coef_summary <-
    list.files(dir_model,
               filename ,
               full.names = TRUE) %>%
    readRDS()  %>%
    filter(coef != "(Intercept)") %>%
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
#----------------------------------------------------------------------
local_coef_plot <- function(input_data) {
  coef_plot <-
    input_data %>%
    # mutate(coef = factor(
    #   coef,
    #   levels = c(
    #     "mean_sm_spring",
    #     "mean_sm_summer",
    #     "mean_sm_spring:mean_tmax_spring",
    #     "mean_sm_spring:mean_sm_summer",
    #     "mean_tmax_spring",
    #     "mean_tmax_summer",
    #     "mean_sm_summer:mean_tmax_summer",
    #     "mean_tmax_spring:mean_tmax_summer"
    #   )
    #   ,
    #   labels = c(
    #     expression("a) SM"[spring]),
    #     expression("b) SM"[summer]),
    #     expression("c) TX:SM"[spring]),
    #     expression("d) SM"[spring:summer]),
    #     expression("e) TX"[spring]),
    #     expression("f) TX"[summer]),
    #     expression("g) TX:SM"[summer]),
    #     expression("h) TX"[spring:summer])
    #   )
    # )) %>%
    mutate(coef = factor(
      coef,
      levels = c(
        "mean_tmax_spring",
        "mean_tmax_summer",
        "mean_sm_spring",
        "mean_sm_summer",
        "mean_sm_spring:mean_tmax_spring",
        "mean_sm_summer:mean_tmax_summer",
        "mean_tmax_spring:mean_tmax_summer",
        "mean_sm_spring:mean_sm_summer"
      )
      ,
      labels = c(
        expression("a) T"[spring]),
        expression("b) T"[summer]),
        expression("c) M"[spring]),
        expression("d) M"[summer]),
        expression("e) T:M"[spring]),
        expression("f) T:M"[summer]),
        expression("g) T"[spring:summer]),
        expression("h) M"[spring:summer])
      )
    )) %>%
    inner_join(spatial_info) %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = mean), color = "transparent") +
    #comment in our out if you want to plot significance
    geom_sf(
      data = . %>%  filter(significant == 0),
      fill = "gray",
      color = "transparent",
      size = 0.00000001,
      alpha = 1
    ) +
    scale_fill_continuous_diverging(
      name = "Standardized ß coefficient",
      "Green-Brown",
      rev = TRUE,
      limit = c(-0.50, 0.50)
    ) +
    theme_void(base_size = 5) +
    facet_wrap( ~ coef, labeller = label_parsed, nrow = 2)  +
    #US coordinates
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
      axis.title.y = element_text(size = 20, angle = 90),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      #axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title.position = "bottom"))
}
#----------------------------------------------------------------------
local_coef_maize <- function(input_data) {
  coef_plot <-
    input_data %>%
    mutate(coef = factor(
      coef,
      levels = c(
        "mean_tmax_spring",
        "mean_tmax_summer",
        "mean_sm_spring",
        "mean_sm_summer",
        "mean_sm_spring:mean_tmax_spring",
        "mean_sm_summer:mean_tmax_summer",
        "mean_tmax_spring:mean_tmax_summer",
        "mean_sm_spring:mean_sm_summer"
      )
      ,
      labels = c(
        expression("i) T"[spring]),
        expression("j) T"[summer]),
        expression("k) M"[spring]),
        expression("l) M"[summer]),
        expression("m) T:M"[spring]),
        expression("n) T:M"[summer]),
        expression("o) T"[spring:summer]),
        expression("p) M"[spring:summer])
      )
    )) %>%
    inner_join(spatial_info) %>%
    st_as_sf() %>%
    ggplot() +
    geom_sf(aes(fill = mean), color = "transparent") +
    #comment in our out if you want to plot significance
    geom_sf(
      data = . %>%  filter(significant == 0),
      fill = "gray",
      color = "transparent",
      size = 0.00000001,
      alpha = 1
    ) +
    scale_fill_continuous_diverging(
      name = "Standardized ß coefficient",
      "Green-Brown",
      rev = TRUE,
      limit = c(-0.50, 0.50)
    ) +
    theme_void(base_size = 5) +
    facet_wrap( ~ coef, labeller = label_parsed, nrow = 2)  +
    #US coordinates
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
      axis.title.y = element_text(size = 20, angle = 90),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      #axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    guides(fill = guide_colorbar(title.position = "bottom"))
}
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
dir_data <- file.path(dir_root, "data")
dir_echam <- file.path(dir_data, "echam_climatology")
dir_storyline <- file.path(dir_data, "echam_storylines")
dir_obs_climate <- file.path(dir_data, "obs_climate_data")
dir_crop <- file.path(dir_data, "crop_data")
dir_shp <- file.path(dir_data, "shapefiles")
dir_model <- file.path(dir_data, "model_data")
dir_figures <- file.path(dir_root, "figures")
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
#load bootstrap coefficient summary per crop incl. CIs
#======================================================================
#load coefficient bootstrap for corn model
boot_corn_model <-
  load_summary_boot_coef_data("coef_ci_corn_std.rds")
#----------------------------------------------------------------------
#load coefficient bootstrap for soy model
boot_soy_model <- load_summary_boot_coef_data("coef_ci_soy_std.rds")
#----------------------------------------------------------------------
#load coefficient bootstrap for wheat us model
boot_wheat_eu_model <-
  load_summary_boot_coef_data("coef_ci_wheat_eu_std.rds")
#----------------------------------------------------------------------
#load coefficient bootstrap for wheat eu model
boot_wheat_us_model <-
  load_summary_boot_coef_data("coef_ci_wheat_us_std.rds")
#======================================================================
#create plots (plot spatially varying coefficients)
#======================================================================
#soy US plot
plot_soy_coef <- local_coef_plot(boot_soy_model) +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125, -75), ylim = c(25, 50)) +
  ggtitle(label = "",
          subtitle = "") +
  ylab("Soybean-US")+
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
#----------------------------------------------------------------------
#maize US plot
plot_corn_coef <- local_coef_maize(boot_corn_model) +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125, -75), ylim = c(25, 50)) +
  ggtitle(label = "",
          subtitle = "") +
  ylab("Maize-US")+
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
#----------------------------------------------------------------------
#Wheat US plot
plot_wheat_us_coef <- local_coef_plot(boot_wheat_us_model) +
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-125,-75), ylim = c(25, 50)) +
  ggtitle(label = "",
          subtitle = "") +
  coord_sf(crs = "+proj=lcc +lon_0=-90 +lat_1=33 +lat_2=45")
#----------------------------------------------------------------------
#Wheat EU plot
plot_wheat_eu_coef <- local_coef_plot(boot_wheat_eu_model) +
  geom_sf(
    data = eu_state_shp,
    color = "black",
    size = 0.00000001,
    fill = "transparent"
  ) +
  coord_sf(xlim = c(-10, 50), ylim = c(30, 70)) +
  ggtitle(label = "",
          subtitle = "") 
#======================================================================
#Patch figures and save them 
#======================================================================
#save soy + maize climate sensitivity plot (version contained plot spacer)
# maize_soy_coef_map <-
#   plot_soy_coef + plot_spacer() + plot_corn_coef +
#   plot_layout(guides = "collect",
#               ncol = 1,
#               heights = c(5, 0.1, 5)) &
#   theme(legend.position = 'bottom')
#----------------------------------------------------------------------
#save soy + maize climate sensitivity plot
maize_soy_coef_map <-
  plot_soy_coef  + plot_corn_coef +
  plot_layout(guides = "collect",
              ncol = 1) &
  theme(legend.position = 'bottom',
        plot.margin = grid::unit(c(0, 0, 0, 0), "mm"))

#----------------------------------------------------------------------
#create png file and store
png(
  file.path(dir_figures, "maize_soy_coef_map.png"),
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)
print(maize_soy_coef_map)
dev.off()
#----------------------------------------------------------------------
#save wheat us climate sensitivity plot
png(
  file.path(dir_figures, "wheat_us_coef_map.png"),
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)
print(plot_wheat_us_coef)
dev.off()
#----------------------------------------------------------------------
#save wheat eu climate sensitivity plot
png(
  file.path(dir_figures, "wheat_eu_coef_map.png"),
  width = 15,
  height = 10,
  units = 'in',
  res = 300
)
print(plot_wheat_eu_coef)
dev.off()

