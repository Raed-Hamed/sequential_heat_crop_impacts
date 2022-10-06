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
library(raster)
library(USAboundaries)
library(colorspace)
library(lme4)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_crop <- file.path(root_dir, "data/usda-nass")
dir_model <-file.path(root_dir, "data/model_df")
source("load_functions.R")
#======================================================================
#load model df
#======================================================================
model_df <-list.files(dir_model,"model_data", full.names = TRUE) %>%
  readRDS() %>% 
  #optional to scale data
  group_by(group_id,County, State) %>% 
  mutate_at(vars(kdd_spring, kdd_summer,gdd_spring, gdd_summer,crop_yield),
            ~ scale(., center = TRUE, scale = TRUE)) %>% 
  ungroup()
#======================================================================
#panel model fit
#======================================================================
panel_model <-lmer(
    crop_yield ~ kdd_summer * kdd_spring + gdd_spring  + (1 |group_id) + (1 |Year),
    data = model_df
  )

panel_model_coefs <-summary_panel$coefficients %>% 
  data.frame()

model_space_predictions <-
  expand.grid(kdd_spring = seq(-4,4,0.05),
              kdd_summer = seq(-4,4,0.05)) %>%
  mutate(
    crop_yield =
      panel_model_coefs[3, 1] * kdd_spring +
      panel_model_coefs[2, 1] * kdd_summer + 
      panel_model_coefs[5, 1] * kdd_spring * kdd_summer 
  ) %>%
  unique()

dependence_plot <-
  #define x,y,z vars for contour lines
  ggplot(data  = model_space_predictions,
         aes(x = kdd_spring,
             y = kdd_summer,
             z = crop_yield)) +
  #add contour lines
  geom_contour(aes(colour = ..level..),
               size = 2,
               show.legend = FALSE) +
  #add points
  # geom_point(
  #   data = obs_df,
  #   aes(x = kdd_spring, y = kdd_summer, fill = crop_yield),
  #   color = "black",
  #   size = 4,
  #   shape = 21,
  #   inherit.aes = FALSE,
  #   alpha = 1
  # ) +
  #add title
#  ggtitle(label = state_name)+
  #add axes
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  #adjust general plot features
  xlab(expression("Standardized spring KDD")) +
  ylab(expression("Standardized summer KDD")) +
  scale_color_continuous_diverging(palette = "Red-Green") +
  scale_fill_continuous_diverging(name = "", palette = "Red-Green") +
  theme_classic(base_size = 20)
#======================================================================
#County scale crop data
#======================================================================
model_pipe <-
  model_df %>%
  group_by(group_id) %>%
  nest() %>%
  mutate(
    model_fit = map(data,  ~ lm(
      crop_yield ~ Year + kdd_summer * kdd_spring + gdd_spring, data = .
    )),
    model_coefs = map(model_fit, ~ broom::tidy(., conf.int = TRUE))
  )
#======================================================================
#extract model details
#======================================================================
#prepare sf layer to spatially plot model outputs
model_sf_to_join <-
  model_df %>%
  dplyr::select(State, County, geometry,group_id) %>%
  distinct()

#extract model coefficients at county scale and join sf layer
model_coefficients <- model_pipe %>%
  dplyr::select(model_coefs, group_id) %>%
  unnest(cols = c(model_coefs)) %>%
  ungroup() %>%
  inner_join(model_sf_to_join) %>%
  filter(!term %in% c("(Intercept)", "Year")) %>%
  drop_na() %>% 
  st_as_sf()

#prepare mask for counties used in the county scale modeling
us_States_Sf <- us_states() %>%
  us_spatial_sf %>%
  tibble() %>%
  inner_join(model_coefficients %>%
               tibble() %>%
               distinct() %>%
               dplyr::select(State)) %>%
  distinct() %>%
  st_as_sf()

#extract and highlight significant cells
model_pvals <- model_pipe %>%
  dplyr::select(model_coefs, group_id) %>%
  unnest(cols = c(model_coefs)) %>%
  ungroup() %>%
  inner_join(model_sf_to_join) %>% 
  st_as_sf() %>% 
  filter(!term %in% c("(Intercept)","Year")) %>% 
  filter(p.value <= 0.05) %>% 
  st_centroid()
#======================================================================
#quick plotting
#======================================================================
#county scale coefs
county_scale_coefs_plot <-county_scale_coefs(model_coefs_df = model_coefficients,
                   model_pvals_df = model_pvals,
                   sf_mask = us_States_Sf)


#distribution of sig model coefficents
model_coefficients %>% 
  tibble() %>% 
  filter(p.value <= 0.1) %>% 
  ggplot()+
  geom_density(aes(estimate), fill = "red", alpha = 0.4)+
  facet_wrap(~term)+
  theme_bw()

#plot dependence plot per state
plot_dependence_per_state(
  observation_df = model_df,
  model_coefs = model_coefficients,
  state_name =  "MONTANA",
  pval_thrsh = 0.1
)

