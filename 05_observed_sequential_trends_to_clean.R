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
source("00_load_functions.R")
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_usda <- file.path(root_dir, "data/usda-nass")
dir_shp <- file.path(root_dir, "data/shapefiles")
dir_model <- file.path(root_dir, "data/model_data")
dir_cmpi6 <- file.path(root_dir, "data/cmip6/")
dir_cmip6_delta_change <-
  file.path(dir_cmpi6, "output/2023_06_13_T_diff/")
dir_cmip6_frequency_change <-
  file.path(dir_cmpi6, "output/2023_06_13_f_change/")
#======================================================================
#load spatial data 
#======================================================================
#spatial info from crop data for plotting
spatial_info <-
  list.files(dir_model,
             "spatial_county_group_id_ref.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  tibble()
#----------------------------------------------------------------------
#define Americas mask (Adjust)
sf_mask_country <- world %>%
  filter(region_un  == "Americas")
#----------------------------------------------------------------------
#us states for plotting
us_state_shp <- list.files(dir_shp, ".json", full.names = TRUE) %>%
  map_dfr(st_read)
#======================================================================
#load USDA datasets and compute mean US harvest estimates for (2010-2021)
#======================================================================
#load USDA harvested yield data
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
 ungroup()
#----------------------------------------------------------------------
#load USDA harvested yield data
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
 ungroup()
#----------------------------------------------------------------------
#load USDA harvested yield data
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
  mutate(Value = Value /  2.471) %>%
  arrange(Year) %>%
  dplyr::select(-Data.Item) %>%
  filter(Value != 0) %>%
  group_by(State, County) %>%
  filter(Year %in% 2010:2021) %>%
  summarise_at(vars(Value), mean) %>%
  ungroup()
#----------------------------------------------------------------------
#load average soy recent yield
avg_recent_soy_yield <- list.files(dir_model,
                                   "avg_recent_soy_yield_us.rds" ,
                                   full.names = TRUE) %>% 
  readRDS()
#----------------------------------------------------------------------
#load average corn recent yield
avg_recent_maize_yield <- list.files(dir_model,
                                     "avg_recent_corn_yield_us.rds" ,
                                     full.names = TRUE) %>% 
  readRDS()
#----------------------------------------------------------------------
#combined crop harvested area
crop_recent_harvested_area <- usda_soy_harvested_area %>%
  mutate(crop = "soy") %>%
  bind_rows(usda_maize_harvested_area %>%
              mutate(crop = "corn"))

#summed crop harvested area
sum_crop_recent_harvested_area <- crop_recent_harvested_area %>%
  pivot_wider(names_from = crop, values_from = Value) %>%
  mutate_at(vars(corn, soy), replace_na, 0) %>%
  mutate(total_harvest_area_pr_county = corn+soy)
#----------------------------------------------------------------------
#combined crop production area
crop_recent_production_avg <-avg_recent_soy_yield%>%
  mutate(crop = "soy") %>%
  bind_rows(avg_recent_maize_yield %>%
              mutate(crop = "corn")) %>% 
  inner_join(crop_recent_harvested_area) %>% 
  mutate(production_total = absolute_yield*Value) %>% 
  group_by(crop) %>% 
  summarise_at(vars(production_total),sum)
#======================================================================
#load model data and standard deviation values per variable
#======================================================================
#load model df with data demeaned at county level
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
  group_by(zone, crop) %>%
  nest() %>%
  mutate(data = map(
    data,
    ~ group_by(., State, County) %>%
      mutate(group_id = cur_group_id() %>%  as.factor())
  )) %>%
  unnest(c(data)) 
#----------------------------------------------------------------------
#scaled model data
model_df_scaled <-model_df %>%
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
#calculate sequential heat frequency
sequential_heat_freq <-model_df %>%
  filter(zone == "US") %>% 
  filter(crop != "wheat") %>% 
  ungroup() %>% 
  dplyr::select(State,
                County,
                Year,
                mean_tmax_spring,
                mean_tmax_summer) %>%
  distinct() %>% 
  group_by(State, County) %>%
  mutate(tmax_spring_75 = mean_tmax_spring %>%
           quantile(0.75) %>%
           round(1)) %>%
  mutate(tmax_summer_75 = mean_tmax_summer %>%
           quantile(0.75) %>%
           round(1)) %>% 
  mutate(warm_spring = ifelse(mean_tmax_spring >= tmax_spring_75,1,0)) %>% 
  mutate(warm_summer = ifelse(mean_tmax_summer >= tmax_summer_75,1,0)) %>% 
  mutate(sequential_heat = ifelse(warm_spring & warm_summer == 1,1,0)) %>% 
  inner_join(sum_crop_recent_harvested_area) %>% 
  group_by(Year) %>% 
  mutate(total_harvest = sum(total_harvest_area_pr_county, na.rm = TRUE)) 
#----------------------------------------------------------------------
sequential_heat_freq %>%
  mutate(spring_heat_area = sum(ifelse(
    warm_spring == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  
  mutate(summer_heat_area = sum(ifelse(
    warm_summer == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  mutate(sequential_heat_area = sum(ifelse(
    sequential_heat == 1, total_harvest_area_pr_county  , 0
  )) /
    total_harvest  * 100) %>%
  distinct() %>%
  pivot_longer(c(spring_heat_area, summer_heat_area, sequential_heat_area)) %>%
  ggplot(aes(x = Year, y = value, color = name)) + geom_line(size = 1) +
  facet_wrap(~ name)+
  ylab("% cropland area under sequential heat (spring & summer > 75th percentile")
#----------------------------------------------------------------------
weighted_avg_model_data <-model_df %>%
  filter(zone == "US") %>% 
  filter(crop != "wheat") %>% 
  ungroup() %>% 
  dplyr::select(State,
                County,
                Year,
                mean_tmax_spring,
                mean_tmax_summer) %>%
  distinct() %>% 
  inner_join(sum_crop_recent_harvested_area) %>% 
  group_by(Year) %>% 
  mutate(normalizer  = 1 / sum(total_harvest_area_pr_county)) %>%
  mutate(weighted_harvest_area = normalizer * total_harvest_area_pr_county) %>%
  dplyr::select(-normalizer) %>% 
  pivot_longer(mean_tmax_spring:mean_tmax_summer) %>% 
  group_by(Year,name) %>%
  nest() %>%
  mutate(weighted_val = map(data, ~ weighted.mean(.$value, .$weighted_harvest_area))) %>%
  dplyr::select(-data) %>%
  unnest(c(weighted_val)) %>%
  ungroup() %>% 
  group_by(name) %>%
  mutate(percentile_75 = weighted_val %>% quantile(0.75) %>%
           round(1)) 
  
                


weighted_avg_model_data %>%  
  mutate(heat_event = ifelse(weighted_val > percentile_75,1,0)) %>% 
  ungroup() %>% 
  dplyr::select(Year,name,heat_event) %>% 
  pivot_wider(names_from = name, values_from = heat_event) %>% 
  mutate(heat_year = ifelse(mean_tmax_spring & mean_tmax_summer ==1,1,0)) %>% 
  filter(heat_year == 1) %>% 
  .$Year %>% 
  unique()
  



weighted_avg_model_data %>%
  ggplot(aes(x = Year, y = weighted_val, color = name)) +
  geom_line(size = 1) +
  geom_hline(aes(yintercept = percentile_75, color = name), linetype = "dashed", size = 1) +
  geom_vline(xintercept = c(1988, 1991, 2006 , 2010 , 2012)) +
  theme_bw()+
  scale_color_manual("",
    values = c("orange", "firebrick3")) +
  theme(
    legend.position = c(0.5, 0.5),
    legend.background = element_rect(fill = alpha("white", 0.3),
                                     color = "transparent"),
    legend.title = element_text(size=10),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 10))+
  ylab("°C")
  
  
  


  
#======================================================================
#load CMIP6 data
#======================================================================
#load TX delta changes
cmip6_tx_delta_change_per_model <-
  list.files(dir_cmip6_delta_change,
             "counties_difference_per_model.shp",
             full.names = TRUE) %>%
  read_sf() %>%
  tibble() %>%
  pivot_longer("245_su_M0":"119_sp_mea") %>%
  dplyr::select(NAME_0,NAME_1,NAME_2,name,value,geometry) %>% 
  separate(name, c("model_scenario", "season","model_ref")) %>%
  mutate(across('season', str_replace, 'su', 'Summer')) %>%
  mutate(across('season', str_replace, 'sp', 'Spring')) %>% 
  filter(!NAME_1  %in% c("Alaska", "Hawaii"))
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
#quick plot of frequency changes
cmip6_tx_freq_change %>% 
  mutate(model_scenario = str_remove_all(model_scenario, "[[:alpha:]]")) %>% 
  filter(reference == "historic") %>% 
  filter(percentile == "75") %>% 
  st_as_sf() %>% 
  ggplot()+
  geom_sf(aes(fill = mean_incre), color = "transparent")+
  geom_sf(data = . %>%  filter(`%M_agree` > 0.5),color = "black",
          size = 0.00000001, fill = "transparent")+
  facet_wrap(~model_scenario)+
  scale_fill_continuous_diverging("Green-Brown")+
  theme_bw()+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, 'cm'),
    legend.title.align = 0.5,
    legend.title = element_text(size = 10, angle = 0),
    legend.key.width = unit(3, "cm"),
    legend.box = "horizontal",
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 10)
    
  )+
  guides(fill = guide_colorbar(title.position = "bottom"))
#======================================================================
#load model coefficients + Model df
#======================================================================
#load local coefficients and adjust units to reflect actual data units
local_coef <-
  list.files(dir_model,
             "local_coef_all_interactions.rds" ,
             full.names = TRUE) %>%
  readRDS() %>%
  tibble() %>%
  filter(zone == "US") %>%
  filter(crop %in% c("corn", "soy")) %>%
  filter(
    coef %in% c(
      "mean_tmax_spring",
      "mean_tmax_summer",
      "mean_tmax_spring:mean_tmax_summer"
    )
  ) %>% 
  dplyr::select(-sd_val) %>% 
  inner_join(model_df_sd_vals)%>% 
  mutate(est_orig_units = std_est*sd_val) %>% 
  dplyr::select(-std_est,-sd_val) %>% 
  pivot_wider(names_from = coef, values_from = est_orig_units)
#----------------------------------------------------------------------
#join coefficients to CMIP6 data (Take only average)
join_coef_climate <-cmip6_tx_delta_change_per_model %>% 
  filter(model_ref == "mea") %>% 
  mutate(State = str_to_upper(NAME_1),
         County = str_to_upper(NAME_2)) %>% 
  dplyr::select(-NAME_1,-NAME_2,-NAME_0) %>% 
  dplyr::select(-geometry) %>% 
  pivot_wider(names_from = season, values_from = value) %>% 
  inner_join(local_coef)
#----------------------------------------------------------------------
#Get future yield projections
yield_prediction_future <-join_coef_climate %>% 
  group_by(State,County,crop,model_ref,CI_level,model_scenario) %>% 
  nest() %>%
  mutate(total_yield_impact = map(
    data,
    ~ with(
      data = .,
      mean_tmax_spring   * Spring   +
        mean_tmax_summer   * Summer  +
        Spring * Summer * `mean_tmax_spring:mean_tmax_summer`
    )
  )) %>%
  mutate(spring_tx_contribution = map(data,
                                      ~ with(data = .,
                                             mean_tmax_spring  * Spring))) %>%
  mutate(summer_tx_contribution = map(data,
                                      ~ with(data = .,
                                             mean_tmax_summer  * Summer))) %>%
  mutate(interaction_contribution = map(
    data,
    ~ with(data = .,
           Spring * Summer * `mean_tmax_spring:mean_tmax_summer`)
  ))
#----------------------------------------------------------------------
#Quick plot
yield_prediction_future %>%
  filter(crop == "corn") %>% 
  filter(model_scenario  == "245") %>% 
  ungroup() %>%
  dplyr::select(-data) %>%
  unnest(c(
    total_yield_impact,
    interaction_contribution,
    summer_tx_contribution,
    spring_tx_contribution
  )) %>% 
  pivot_longer(total_yield_impact:interaction_contribution) %>% 
  mutate(name = factor(name, levels = c("spring_tx_contribution",
                                        "summer_tx_contribution",
                                        "interaction_contribution",
                                        "total_yield_impact"
  ))) %>% 
  inner_join(avg_recent_maize_yield) %>% 
  #filter(name == "interaction_contribution") %>% 
  mutate(prct_of_recent_yield = value/absolute_yield*100) %>% 
  inner_join(spatial_info) %>% 
  filter(CI_level == "point_est") %>% 
  st_as_sf() %>% 
  ggplot()+
  geom_sf(
    data = sf_mask_country,
    color = "black",
    size = 0.00001,
    fill = "lightgray"
  ) +
  geom_sf(aes(fill = value), color = "transparent")+
  geom_sf(
    data = us_state_shp,
    color = "black",
    size = 0.00001,
    fill = "transparent"
  ) +
  facet_wrap(~name)+
  scale_fill_continuous_diverging(name = "Yield anomaly (t/ha)",
                                  "Red-Green")+
  theme_bw(base_size = 5) +
  coord_sf(xlim = c(-125, -75), ylim = c(25, 48)) +
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(1, 'cm'),
    legend.title.align = 0.5,
    legend.title = element_text(size = 10, angle = 0),
    legend.key.width = unit(3, "cm"),
    legend.box = "horizontal",
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 10)
    
  )+
  guides(fill = guide_colorbar(title.position = "bottom"))
# theme(
#   legend.position = "right",
#   legend.key.size = unit(2.6, 'cm'),
#   legend.title.align = 0.5,
#   legend.title = element_text(size = 8, angle = -90),
#   legend.key.width = unit(1, "cm"),
#   legend.box = "vertical",
#   strip.text.x = element_text(size = 10),
#   strip.text.y = element_text(size = 10),
#   legend.text = element_text(size = 10)
#   
# ) +
# guides(fill = guide_colorbar(title.position = "right"))
#=========================================================





production_prediction <-yield_prediction_future%>%
  ungroup() %>%
  dplyr::select(-data) %>%
  unnest(c(
    total_yield_impact,
    interaction_contribution,
    summer_tx_contribution,
    spring_tx_contribution
  )) %>% 
  pivot_longer(total_yield_impact:interaction_contribution ) %>% 
  mutate(name = factor(name, levels = c("spring_tx_contribution",
                                        "summer_tx_contribution",
                                        "interaction_contribution",
                                        "total_yield_impact"
  ))) %>% 
  inner_join(crop_recent_harvested_area) %>% 
  mutate(production = value*Value) %>% 
  group_by(model_scenario,crop,CI_level,name) %>% 
  summarise_at(vars(production), sum) %>% 
  ungroup() %>%
  inner_join(crop_recent_production_avg) %>% 
  mutate(production = production/production_total*100) %>% 
  pivot_wider(names_from =CI_level, values_from =production) 


ggplot() +
  # geom_rect(
  #   data = production_prediction,
  #   aes(
  #     xmin = "Central Brazil",
  #     xmax = "United States",
  #     ymin = 0,
  #     ymax = value
  #   ),
  #   fill = "gray",
  #   alpha = 0.05,
  #   color = "transparent"
# ) +
geom_bar(data = production_prediction,
         aes(
           x = model_scenario,
           y = point_est,
           fill = name,
           group = name
         ),
         stat = "identity",
         position = "dodge") +
  geom_errorbar(
    data = production_prediction,
    aes(x = model_scenario,
        group = name,
        ymin = lwr,
        ymax = upr),
    width = 0.1,
    position = position_dodge(width = 0.9)
  ) +
  geom_point(
    data = production_prediction,
    aes(x = model_scenario,y = point_est,group =name
    ),color = "black",
    position = position_dodge(width = 0.9)
  ) +
  # geom_point(
  #   data = delta_change_production %>%
  #     group_by(experiment) %>%
  #     mutate_at(vars(value), sum),
  #   aes(
  #     x = 2.1,
  #     y = value
  #   ),color = "darkgray"
  # ) +
  geom_vline(xintercept = 3.5, size = 1, linetype = "dashed")+
  facet_grid(~crop, scale = "free",
             labeller = labeller(crop = c("corn" = "Maize",
                                          "soy" = "Soybean"))) +
  xlab("SSP") +
  ylab("Percentage of recent total production (2010-2021)") +
  geom_hline(yintercept = 0, linewidth = 1) +
  #  scale_fill_manual("", values = c("firebrick", "tan", "#FFBF00")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.15, 0.20),
    legend.background = element_rect(fill = alpha("white", 1),
                                     color = "black"),
    legend.text=element_text(size=12),
    axis.title.y = element_text(size = 12),
    strip.text.x = element_text(size = 15)#,
    # axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

