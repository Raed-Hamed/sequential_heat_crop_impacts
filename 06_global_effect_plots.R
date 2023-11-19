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
#======================================================================
#load directories
#======================================================================
dir_root <- getwd()
dir_data <- file.path(dir_root, "data")
dir_model <- file.path(dir_data, "model_data")
dir_figures <- file.path(dir_root, "figures")
dir_usda <- file.path(dir_root, "data/usda-nass")
dir_eurostat <- file.path(dir_root, "data/eurostat")
#======================================================================
#helper functions
#======================================================================
bivariate_interaction_plot <- function(model, observation, plot_title) {
  ggplot(data  = model,
         aes(x = tx_summer,
             y = tx_spring,
             z = yield_anomaly_predict)) +
    stat_summary_2d(
      data =  observation,
      aes(x = mean_tmax_summer, y = mean_tmax_spring , z = crop_yield),
      color = "darkgray",
      size = 0.005,
      binwidth = c(0.25, 0.25),
      fun = mean,
      inherit.aes = FALSE,
      alpha = 1
    ) +
    #add contour lines
    geom_contour(
      aes(colour = ..level..),
      size = 1.5,
      linejoin = "round",
      show.legend = FALSE,
      bins = 15
    ) +
    geom_contour(
      data = model %>% 
        mutate(yield_anomaly_predict_rnd = round(yield_anomaly_predict)) %>% 
        filter(yield_anomaly_predict_rnd  == 0),
      aes(colour = ..level..),
      color = "black",
      size = 1.5,
      linejoin = "round",
      linetype = "twodash",
      show.legend = FALSE,
      bins = 2
    ) +
    #add title
    ggtitle(label = plot_title)+
    #add axes
    geom_hline(yintercept = 0, size = 1, color = "darkgray") +
    geom_vline(xintercept = 0, size = 1, color = "darkgray") +
    #adjust general plot features
    xlab("Summer maximum temperature anomaly (°C)") +
    ylab("Spring maximum temperature anomaly (°C)") +
    scale_color_continuous_diverging(palette = "Red-Green"#,
                                     #limit = c(-50, 50)
                                     ) +
    scale_fill_continuous_diverging(name = "Detrended yield anomaly (t/ha)",
                                    palette = "Red-Green"#,
                                    #limit = c(-50, 50)
                                    ) +
    ylim(c(-4.5, 4.5)) +
    xlim(c(-4.5, 4.5)) +
    theme_bw(base_size = 20) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.size = unit(1, 'cm'),
      legend.title.align = 0.5,
      legend.title = element_text(size = 15, angle = 0),
      legend.key.width = unit(2, "cm"),
      legend.box = "horizontal",
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      legend.text = element_text(size = 15)
      
    ) +
    guides(fill = guide_colorbar(title.position = "bottom"))
  
}
#----------------------------------------------------------------------
delta_slope_plot <- function(crop,plot_title) {
  
  ggplot(
    delta_slope_effect %>%
      filter(area == crop),
    aes(
      y = estimate ,
      x = summer_tx,
      color = as.factor(spring_tx)
    )
  ) +
    #add v-line
    # geom_vline(xintercept = 0, linetype = "dashed") +
    # #add h-line
    # geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0,color = "black") +
    geom_vline(xintercept = 0, color = "black") +
    #add line
    geom_line(size = 1, linetype = "twodash") +
    geom_ribbon(
      aes(ymin = lwr ,
          ymax = upr,
          fill = spring_tx),
      linetype = 1,
      size = 0.1,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    #Adjust legend
    scale_color_manual(
      name = "Spring maximum temperature anomaly (°C)",
      values = c("darkcyan", "#999999", "firebrick3"),
      labels =
        c("Cold (5%)",
          "Normal (50%)",
          "Hot (95%)")
    ) +
    guides(linetype = "none") +
    #adjust CI
    scale_fill_manual(values = c("darkcyan", "#999999", "firebrick3"),
                      labels = c("")) +
    #adjust general plot features
    #ylab("Detrended yield anomaly (t/ha)") +
    ylab("Percentage of recent yield (%)") +
    xlab("Summer maximum temperature anomaly (°C)") +
    #add title
    ggtitle(label = plot_title)+
    theme_bw(base_size = 20) +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = alpha("white", 1),
                                       color = "black"),
      legend.title = element_text(size = 15),
      strip.text.x = element_text(size = 15),
      strip.text.y = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 15)
    )  +
    xlim(c(-4.5, 4.5)) +
    ylim(c(-25,25))
}

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
#load average soy recent yield
global_absolute_yield <-
  list.files(dir_model,
             "global_absolute_recent_yield.rds" ,
             full.names = TRUE) %>%
  readRDS() %>% 
  group_by(State,County,crop) %>% 
  summarise_at(vars(absolute_yield),mean) %>% 
  ungroup()
#----------------------------------------------------------------------
#load average corn recent yield
recent_avg_yield_all <-
  inner_join(global_harvested_area, global_absolute_yield) %>%
  group_by(crop) %>%
  nest() %>%
  mutate(yield_est = map(
    data,
    ~ weighted.mean(.$absolute_yield, .$weighted_harvest_area)
  )) %>%
  unnest(yield_est) %>%
  dplyr::select(-data) %>%
  mutate(crop = case_when(crop == "soy" ~ "soy_us",
                          crop == "corn" ~ "corn_us",
                          TRUE ~ crop)) %>% 
rename(area = crop)
#======================================================================
#Load model data and output
#======================================================================
#load model data
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
  filter(County != "AT34") %>% 
  group_by(zone, crop) %>%
  nest() %>%
  mutate(data = map(
    data,
    ~ group_by(., State, County) %>%
      mutate(group_id = cur_group_id() %>%  as.factor())
  )) %>%
  unnest(c(data)) %>%
  group_by(County, State, zone, crop, group_id) %>%
  mutate_at(
    vars(
      crop_yield,
      mean_sm_spring,
      mean_sm_summer,
      mean_tmax_spring,
      mean_tmax_summer
    ),
    ~ scale(., center = TRUE, scale = FALSE)
  )
#----------------------------------------------------------------------
#store spring tx quantile values
spring_tx_quantile_values <-
  model_df %>%
  group_by(crop,zone) %>%
  summarise_at(vars(c(mean_tmax_spring)),
               ~ data.frame(
                 q05 = quantile(., 0.05),
                 q50 = quantile(., 0.50),
                 q95 = quantile(., 0.95)
               )) %>%
  mutate(q05 = mean_tmax_spring$q05,
         q50 = mean_tmax_spring$q50,
         q95 = mean_tmax_spring$q95) %>%
  dplyr::select(-mean_tmax_spring)
#----------------------------------------------------------------------
#load model outputs
model_output <-  list.files(dir_model,
                              "model_.+full.rds" ,
                              full.names = TRUE) %>%
  map(readRDS)
#======================================================================
#store fixed effects and create model space predictions
#======================================================================
fixed_effects <- tibble(
  area = c("corn_us", "soy_us", "wheat_eu", "wheat_us"),
  model_fit = c(model_output)
) %>% 
  group_by(area) %>% 
  mutate(model_fit = map(model_fit, ~tidy(.))) %>% 
  unnest(c(model_fit))%>% 
  filter(effect == "fixed") %>% 
  mutate(upr = estimate + std.error*1.96) %>% 
  mutate(lwr = estimate - std.error*1.96) %>% 
  dplyr::select(-group,-statistic)  
#----------------------------------------------------------------------
#create model space predictions (Bivariate plot)
model_space_predictions <-
  expand.grid(tx_summer = seq(-4.5, 4.5, 0.05),
              tx_spring = seq(-4.5, 4.5, 0.05)) %>%
  crossing(fixed_effects) %>%
  dplyr::select(tx_summer, tx_spring, area, term, estimate) %>%
  pivot_wider(names_from = "term", values_from = "estimate") %>%
  group_by(area) %>%
  nest() %>%
  mutate(yield_anomaly_predict = map(
    data,
    ~ with(
      .,
      mean_tmax_spring * tx_spring +
        mean_tmax_summer * tx_summer +
        `mean_tmax_spring:mean_tmax_summer` *
        tx_spring * tx_summer
    )
  )) %>%
  unnest(c(data, yield_anomaly_predict)) %>%
  ungroup() 
#======================================================================
#plot bivariate prediction space
#======================================================================
#WHEAT-EU plot
bivariate_wheat_eu_plot <- bivariate_interaction_plot(
  model = model_space_predictions %>%
    filter(area == "wheat_eu") ,
    # mutate(yield_anomaly_predict = 
    #          yield_anomaly_predict / 
    #          avg_recent_wheat_yield_eu %>%
    #          pull *100)
  observation = model_df %>% 
    filter(zone == "EU") %>%
    filter(County != "AT34") %>%
    # mutate(crop_yield = 
    #          crop_yield / 
    #          avg_recent_wheat_yield_eu %>%
    #          pull *100) %>% 
    filter(crop == "wheat"),
  plot_title = "b) Wheat-EU"
) + geom_hline(
  data = spring_tx_quantile_values %>%
    filter(zone == "EU"),
  aes(yintercept = q05),
  size = 1.5,
  color = "darkcyan",
  linetype = "twodash"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(zone == "EU"),
  aes(yintercept = q95),
  size = 1.5,
  color = "firebrick3",
  linetype = "twodash"
)

#----------------------------------------------------------------------
#WHEAT-EU plot
bivariate_wheat_us_plot <-bivariate_interaction_plot(
  model = model_space_predictions %>%
    filter(area == "wheat_us"),
  observation = model_df %>%  filter(zone == "US") %>%
    filter(crop == "wheat"),
  plot_title = "a) Wheat-US"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(zone == "US") %>% 
    filter(crop == "wheat"),
  aes(yintercept = q05),
  size = 1.5,
  color = "darkcyan",
  linetype = "twodash"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(zone == "US") %>% 
    filter(crop == "wheat"),
  aes(yintercept = q95),
  size = 1.5,
  color = "firebrick3",
  linetype = "twodash"
)
#----------------------------------------------------------------------
#MAIZE-US plot
bivariate_corn_us_plot <-
  bivariate_interaction_plot(
  model = model_space_predictions %>%
    filter(area == "corn_us"),
    # mutate(yield_anomaly_predict = 
    #          yield_anomaly_predict / 
    #          avg_recent_maize_yield %>%
    #          pull *100)
  observation = model_df  %>%
    # mutate(crop_yield = 
    #          crop_yield / 
    #          avg_recent_maize_yield %>%
    #          pull *100) %>% 
    filter(crop == "corn"),
  plot_title = "b) Maize-US"
) + geom_hline(
  data = spring_tx_quantile_values %>%
    filter(crop  == "corn"),
  aes(yintercept = q05),
  size = 1.5,
  color = "darkcyan",
  linetype = "twodash"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(crop == "corn"),
  aes(yintercept = q95),
  size = 1.5,
  color = "firebrick3",
  linetype = "twodash"
)

#----------------------------------------------------------------------
#Soy-US plot
bivariate_soy_us_plot <-
bivariate_interaction_plot(
  model = model_space_predictions %>%
    filter(area == "soy_us") , 
    # mutate(yield_anomaly_predict = 
    #          yield_anomaly_predict / 
    #          avg_recent_soy_yield %>%
    #          pull *100)
  observation = 
    model_df %>%  
    # mutate(crop_yield = 
    #          crop_yield / 
    #          avg_recent_soy_yield %>%
    #          pull *100) %>% 
    filter(crop == "soy"),
  plot_title = "a) Soybean-US"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(crop  == "corn"),
  aes(yintercept = q05),
  size = 1.5,
  color = "darkcyan",
  linetype = "twodash"
)+ geom_hline(
  data = spring_tx_quantile_values %>%
    filter(crop == "corn"),
  aes(yintercept = q95),
  size = 1.5,
  color = "firebrick3",
  linetype = "twodash"
)
#======================================================================
#calculate sequential heat delta slope effect
#======================================================================
delta_slope_effect <- model_df %>%
  ungroup() %>%
  mutate(zone = str_to_lower(zone)) %>%
  unite("area", c(crop, zone)) %>%
  inner_join(fixed_effects) %>%
  dplyr::select(area,
                State,
                County,
                term,
                mean_tmax_spring,
                mean_tmax_summer,
                estimate,
                upr,
                lwr)  %>%
  mutate(term = paste0("coef_", term)) %>%
  pivot_longer(estimate:lwr) %>%
  distinct() %>%
  pivot_wider(names_from = "term", values_from = "value") %>%
  group_by(area, name) %>%
  nest() %>%
  mutate(predictor_space = map(
    data,
    ~ expand.grid(
      summer_tx = .$mean_tmax_summer,
      spring_tx = c(
        .$mean_tmax_spring %>%
          quantile(0.05) %>%
          round(1),
        .$mean_tmax_spring %>%
          quantile(0.5) %>%
          round(1),
        .$mean_tmax_spring %>%
          quantile(0.95) %>%
          round(1)
      )
    )
  )) %>%
  mutate(
    predict_yield_anomaly = map2(
      predictor_space,
      data,
      ~ .x$spring_tx * .y$coef_mean_tmax_spring +
        .x$summer_tx * .y$coef_mean_tmax_summer +
        .x$spring_tx * .x$summer_tx * .y$`coef_mean_tmax_spring:mean_tmax_summer`
    )
  ) %>%
  dplyr::select(-data) %>%
  unnest(c(predict_yield_anomaly, predictor_space)) %>%
  ungroup() %>%
  mutate(spring_tx = spring_tx %>%  names() %>%
           as.factor()) %>%
  distinct() %>%
  pivot_wider(names_from = name,
              values_from = predict_yield_anomaly) %>% 
  inner_join(recent_avg_yield_all) %>% 
  mutate(estimate = estimate / yield_est *100,
         upr = upr/yield_est*100,
         lwr = lwr/yield_est*100)
#======================================================================
#plot sequential heat delta slope effect
#======================================================================
#delta slope plot maize US
maize_us_ds <-delta_slope_plot(crop = "corn_us", plot_title = "d) Maize-US")+
  guides(color = guide_legend(override.aes = list(size = 3))) 

#------------------------------------------------------------------------
#delta slope plot maize US
soy_us_ds <-delta_slope_plot(crop = "soy_us", plot_title = "c) Soybean-US")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
#------------------------------------------------------------------------
#delta slope plot Wheat US
wheat_us_ds <-delta_slope_plot(crop = "wheat_us", plot_title = "c) Wheat-US")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
#------------------------------------------------------------------------
#delta slope plot maize US
wheat_eu_ds <-delta_slope_plot(crop = "wheat_eu", plot_title = "d) Wheat-EU")+
  guides(color = guide_legend(override.aes = list(size = 3))) 
#======================================================================
#Bring plots together and save
#======================================================================
maize_soy_plot <- (bivariate_soy_us_plot + bivariate_corn_us_plot) /
  ((soy_us_ds  + maize_us_ds) + plot_layout(guides = "collect") &
     theme(legend.position = 'bottom')
  )

png(
  file.path(dir_figures, "bivariate_maize_soy.png"),
  width = 15,
  height = 15,
  units = 'in',
  res = 300
)
print(maize_soy_plot)
dev.off()
#------------------------------------------------------------------------
wheat_eu_us_plot <- (bivariate_wheat_us_plot + bivariate_wheat_eu_plot) /
  ((wheat_us_ds  + wheat_eu_ds) + plot_layout(guides = "collect") &
     theme(legend.position = 'bottom')
  )

png(
  file.path(dir_figures, "bivariate_wheat_us_eu.png"),
  width = 14,
  height = 15,
  units = 'in',
  res = 300
)
print(wheat_eu_us_plot)
dev.off()

