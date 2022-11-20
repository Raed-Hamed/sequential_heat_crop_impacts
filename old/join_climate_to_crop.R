##TO REMOVE
##

#======================================================================
#load libraries
#======================================================================
library(tidyverse)
library(sf)
library(raster)
library(USAboundaries)
library(colorspace)
#======================================================================
#clean workspace
#======================================================================
# rm(list = ls())
# graphics.off()
# gc()
#======================================================================
#load directories
#======================================================================
root_dir <- "C:/Users/rhd630/Desktop/PhD/Academic/paper_4"
dir_data <- file.path(root_dir, "data/cpc")
#======================================================================
#load cpc kdd spring dataset
kdd_spring <- list.files(dir_data, "kdd", full.names = TRUE) %>%
  brick(varname = "kdd_spring") %>%
  crop(usda_shp) %>%
  mask(usda_shp) %>%
  rasterToPoints() %>%
  data.frame() %>%
  tibble() %>%
  pivot_longer(X1979.1:X1979.43) %>%
  group_by(x, y) %>%
  mutate(Year = 1979:2021) %>%
  dplyr::select(-name) %>%
  ungroup()
#======================================================================
#load cpc kdd summer dataset
kdd_summer <- list.files(dir_data, "kdd", full.names = TRUE) %>%
  brick(varname = "kdd_summer") %>%
  crop(usda_shp) %>%
  mask(usda_shp) %>%
  rasterToPoints() %>%
  data.frame() %>%
  tibble() %>%
  pivot_longer(X1979.1:X1979.43) %>%
  group_by(x, y) %>%
  mutate(Year = 1979:2021) %>%
  dplyr::select(-name) %>%
  ungroup()
#======================================================================
kdd_data_full <- kdd_spring %>%
  rename(kdd_spring = value) %>%
  inner_join(kdd_summer) %>%
  rename(kdd_summer = value)
#======================================================================
#kdd as sf
kdd_full_sf <- kdd_data_full %>%
  st_as_sf(x = .,
           coords = c("x", "y"),
           crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
#======================================================================
#take crop data and filter for years at least 33
usda_wheat_sf$Year %>%  unique() %>%  length() * .70

clean_crop_data <- usda_wheat_sf %>%
  filter(n > 30) %>%
  dplyr::select(-n) %>%
  rename(crop_yield = Value)
#======================================================================
#Aggregate grid based climate information to county scale
df_point_list <- split(dplyr::select(kdd_full_sf, -Year), 
                       kdd_full_sf$Year)

#Split sf usda_resid country by Years
df_poly_list <- split(clean_crop_data, clean_crop_data$Year)

#join yield and climate sf, na are reproduced as climate_df is filtered
full_Sf_yield_climate <- map2_dfr(
  df_poly_list,
  df_point_list,
  ~ .x %>%
    st_join(.y, left = FALSE) %>%
    group_by(County, Year, State) %>%
    mutate_at(vars(-geometry,-State,-County,-Year),  ~
                (mean(., na.rm = TRUE)))
)
#======================================================================
#create ID per county,state combination
full_data_id <- full_Sf_yield_climate %>%
  group_by(County, State) %>%
  mutate(group_id = cur_group_id()) %>%
  group_by(group_id,County, State) %>% 
  mutate_at(vars(kdd_spring, kdd_summer, crop_yield),
            ~ scale(., center = TRUE, scale = TRUE)) %>% 
  ungroup()

#======================================================================
#detrended_full_data_id <-full_data_id

saveRDS(full_data_id, "model_data.rds")
#======================================================================
#panel model fit
#======================================================================
library(lme4)
est <-
  lmer(
    crop_yield ~ kdd_summer * kdd_spring + (1 |group_id) + (1 |Year),
    data = full_data_id
  )

summary(est)
#======================================================================
#County scale crop data
#======================================================================
model_pipe <-
  full_data_id %>%
  group_by(group_id) %>%
  nest() %>%
  mutate(
    model_fit = map(data,  ~ lm(
      crop_yield ~ Year + kdd_summer * kdd_spring, data = .
    )),
    model_coefs = map(model_fit, ~ broom::tidy(., conf.int = TRUE))
  )
#======================================================================
sf_join_data <-
  full_data_id %>%
  dplyr::select(State, County, geometry,group_id) %>%
  distinct()
#======================================================================
model_coefficients <- model_pipe %>%
  dplyr::select(model_coefs, group_id) %>%
  unnest(cols = c(model_coefs)) %>%
  ungroup() %>%
  inner_join(sf_join_data) %>%
  filter(!term %in% c("(Intercept)", "Year")) %>%
  drop_na() %>% 
  st_as_sf()
#======================================================================
#signigicant cells
model_pvals <- model_pipe %>%
  dplyr::select(model_coefs, group_id) %>%
  unnest(cols = c(model_coefs)) %>%
  ungroup() %>%
  inner_join(sf_join_data) %>% 
  st_as_sf() %>% 
  filter(!term %in% c("(Intercept)","Year")) %>% 
  filter(p.value <= 0.05) %>% 
  st_centroid()
#======================================================================
#load yield crop mask
us_States_Sf <- us_states() %>%
  us_spatial_sf %>%
  tibble() %>%
  inner_join(model_coefficients %>%
               tibble() %>%
               distinct() %>%
               dplyr::select(State)) %>%
  distinct() %>%
  st_as_sf()
#======================================================================
#plotting coefficient destribution
model_coefficients %>% 
  ggplot() + 
  geom_histogram(aes(x=estimate), fill = "red") +
  facet_grid(~term)+
  theme_bw()
#======================================================================
library(RColorBrewer)
breaks_plot <-seq(-0.5,0.5,0.1)
n <- breaks_plot %>%  length()
getPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
colors_plot <- getPalette(n) #%>%  rev()
limit <-c(-0.5,0.5)
#======================================================================
ggplot() +
  geom_sf(data = us_States_Sf, fill = "darkgrey") +
  geom_sf(data = model_coefficients,
          aes(fill = estimate),
          color = "transparent") +
  geom_sf(data = model_pvals, size = 0.1) +
  facet_wrap( ~ term) +
  geom_sf(
    data = model_coefficients %>%
      dplyr::select(geometry) %>%
      tibble() %>%
      distinct() %>%
      st_as_sf(),
    fill = "transparent",
    size = 0.0001,
    color = "transparent",
    alpha = 0.5
  ) +
  scale_fill_gradient2(
    name = "β coefficient",
    midpoint = 0,
    low = "firebrick4",
    mid = "white",
    high = "forestgreen",
    breaks = seq(-0.5, 0.5, 0.1),
    limits = c(-0.5,  0.5)
  ) +
  # scale_fill_stepsn(
  #   name = "β coefficient",
  #   colours = colors_plot,
  #   breaks = breaks_plot,
  #   limit = limit,
  #   labels = function(x) {
  #     x <- as.character(x %>%  round(1))
  #     x[!1:0] <- " "
  #     x
  #   }
  # ) +
  geom_sf(data = us_States_Sf, fill = "transparent") +
  theme_bw(base_size = 15) +
  theme(
    legend.text = element_text(size = 13),
    legend.key.width = unit(1, "cm"),
    legend.position = "bottom",
    legend.spacing.x = unit(0, 'cm')
  ) +
  guides(fill = guide_legend(
    label.position = "bottom",
    nrow = 1,
    override.aes = list(size = 1.5)
  ))
#======================================================================
model_coefficients %>% 
  tibble() %>% 
  filter(p.value <= 0.1) %>% 
  ggplot()+
  geom_density(aes(estimate), fill = "red", alpha = 0.4)+
  facet_wrap(~term)+
  theme_bw()
#north dakota / Montana
#======================================================================
plot_dependence_per_state <-function(state_name,pval_thrsh){

   # state_name <- "MONTANA"
   # pval_thrsh <- 0.1
  
   sign_counties <-model_coefficients %>% 
   filter(p.value <= pval_thrsh) %>% 
     tibble() %>% 
     .$County
   
  #Get coefficients for model with interactions for a specific state
  coefs_df <- model_coefficients %>% 
    filter(p.value <= pval_thrsh) %>% 
    dplyr::select(County, State, term, estimate) %>%
    group_by(term) %>%
    filter(State == state_name) %>% 
    tibble() %>% 
    group_by(term) %>% 
    summarise_at(vars(estimate),~mean(., na.rm = TRUE))
  
  #Get coefficients with confidence bounds
  coefs_df_confidence_interval <- model_coefficients%>%
    tibble() %>% 
    filter(p.value <= pval_thrsh) %>% 
    dplyr::select(County, State, term, estimate, conf.low, conf.high) %>%
    group_by(term) %>%
    filter(State == state_name) %>% 
    summarize_at(vars(-County, -State),  ~ (mean(., na.rm = TRUE)))
  
  #get points montana
  obs_df <-full_data_id %>% 
    ungroup() %>% 
    filter(County %in% sign_counties) %>% 
    dplyr::select(-geometry,-group_id)%>% 
    na.exclude() %>% 
    filter(State == state_name) %>%
    unique()
  
  #model space predictions
  model_space_predictions <-
    expand.grid(kdd_spring = seq(-4,4,0.05),
                kdd_summer = seq(-4,4,0.05)) %>%
    mutate(
      crop_yield =
        pull(coefs_df[1, 2]) * kdd_spring +
        pull(coefs_df[2, 2]) * kdd_summer + 
        pull(coefs_df[3, 2]) * kdd_spring * kdd_summer 
    ) %>%
    unique()
  
  #interaction plot
  interaction_df <-
    expand.grid(
      kdd_summer = obs_df$kdd_summer,
      kdd_spring = c(
        obs_df$kdd_spring %>%  quantile(0.05) %>%  round(1),
        obs_df$kdd_spring %>%  quantile(0.5) %>%  round(1),
        obs_df$kdd_spring %>%  quantile(0.95) %>%  round(1)
      )
    ) %>%
    mutate(
      crop_yield =
        pull(coefs_df_confidence_interval[1, 2]) * kdd_spring +
        pull(coefs_df_confidence_interval[2, 2]) * kdd_summer + 
        pull(coefs_df_confidence_interval[3, 2]) * kdd_spring * kdd_summer ,
      
      crop_yield.low =
        pull(coefs_df_confidence_interval[1, 3]) * kdd_spring +
        pull(coefs_df_confidence_interval[2, 3]) * kdd_summer + 
        pull(coefs_df_confidence_interval[3, 3]) * kdd_spring * kdd_summer ,
      
      crop_yield.high =
        pull(coefs_df_confidence_interval[1, 4]) * kdd_spring +
        pull(coefs_df_confidence_interval[2, 4]) * kdd_summer + 
        pull(coefs_df_confidence_interval[3, 4]) * kdd_spring * kdd_summer ,
    ) %>% 
    unique() %>% 
    as.data.frame() %>%   
    tibble() %>% 
    mutate(kdd_spring = factor(paste(
      kdd_spring)))

  #Change factor names
  levels(interaction_df$kdd_spring) <- c("percentile_5th",
                                         "percentile_50th",
                                         "percentile_95th")
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
    geom_point(
      data = obs_df,
      aes(x = kdd_spring, y = kdd_summer, fill = crop_yield),
      color = "black",
      size = 4,
      shape = 21,
      inherit.aes = FALSE,
      alpha = 1
    ) +
    #add title
    ggtitle(label = state_name)+
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
  interaction_plot <-
    ggplot(interaction_df ,
           aes(
             y = crop_yield,
             x = kdd_summer,
             color = kdd_spring,
             linetype = kdd_spring
           )) +
    #add line
    geom_line(size = 2) +
    #add CI
    geom_ribbon(
      data = interaction_df,
      aes(
        ymin = crop_yield.low,
        ymax = crop_yield.high,
        fill = kdd_spring
      ),
      linetype = 0,
      alpha = 0.2,
      show.legend = FALSE
    ) +
    #Adjust legend
    scale_color_manual("Spring KDD",
                       values = c("#bc7124", "#999999", "#209281"),
                       labels = c("Extreme cool (5%)", "Normal (50%)", "Extreme warm (95%)")
    ) +
    #Adjust linetype
    scale_linetype_manual(values = c("dashed", "solid" , "dashed")) +
    guides(linetype = FALSE) +
    #adjust CI
    scale_fill_manual(values = c("#bc7124", "#999999", "#209281"),
                      labels = c("")) +
    #adjust general plot features
    ylab("Standardized Yield Anomaly") +
    xlab("Standardized summer KDD") +
    theme_classic(base_size = 20) +
    theme(
      legend.position = c(0.40, 0.20),
      legend.background = element_rect(fill = alpha("white", 0))
    )+
    xlim(-2, 2)
  
  library(patchwork)
  conmbined_plot <- dependence_plot +
    interaction_plot +
    plot_annotation(tag_levels = list(c("(a)", "(b)")))
  
  return(conmbined_plot)
}
#======================================================================
plot_dependence_per_state("KANSAS", 0.1)
#======================================================================
#incl. gdd
#detrend at county scale