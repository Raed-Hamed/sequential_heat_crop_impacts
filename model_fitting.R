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
library(spData)
#======================================================================
#load directories
#======================================================================
root_dir <- getwd()
dir_data <- file.path(root_dir, "data/cpc")
dir_crop <- file.path(root_dir, "data/usda-nass")
dir_eurostat <- file.path(root_dir, "data/eurostat")
dir_model <-file.path(root_dir, "data/model_data")
source("load_functions.R")
#======================================================================
#load model df
#======================================================================
model_df <-list.files(dir_model,"model_data_maize.rds", full.names = TRUE) %>%
  readRDS() %>% 
  drop_na() %>% 
  #optional to scale data
  group_by(group_id,County,State) %>%
  mutate_at(
    vars(
      kdd_spring,
      kdd_summer,
      gdd_spring,
      gdd_summer,
      crop_yield,
      mean_sm_spring,
      mean_tmax_spring,
      mean_sm_summer,
      mean_tmax_summer
    ),
    ~ scale(., center = TRUE, scale = TRUE)) %>%
  ungroup() %>%
  filter(n >= 30) %>% 
  distinct()
#======================================================================
#County scale model fit
#======================================================================
model_pipe <-
  model_df %>%
  drop_na() %>%
  group_by(County, State, group_id) %>%
  nest() %>%
  mutate(
    model_fit = map(data,  ~ lm(
      crop_yield ~ mean_sm_spring * mean_sm_summer, data = .
    )),
    model_coefs = map(model_fit, ~ broom::tidy(., conf.int = TRUE)),
    fitted_vals = map(model_fit,  ~ .$fitted.values),
    rsq_full_model = map2(data, fitted_vals, ~ rsq(.x$crop_yield, .y)),
    residuals = map(model_fit, ~ residuals(.)),
    #loo_model = map2(data, model_fit, ~ fun_loo_per_grid(.x, .y)),
    #rsq_loo_model = map2(data, loo_model, ~ rsq(.x,$crop_yield, .y)),
    model_pvalue = map(model_fit, ~ lmp(.)),
    #diagnostics significance tests
    dwtest = map(
      model_fit,
      ~ lmtest::dwtest(., alternative = "two.sided")$p.value
    ),
    bptest = map(model_fit, ~ lmtest::bptest(.)$p.value),
    vif_val = map(model_fit, ~ car::vif(.) %>%
                    tibble(var = names(.), vif_val = .)),
    shapirotest = map(residuals, ~ shapiro.test(.)$p.value),
    resettest = map2(
      model_fit,
      data,
      ~ lmtest::resettest(
        formula(.x),
        data = .y,
        type = "regressor",
        power = 2
      )$p.value
    )
    
  )
#======================================================================
#Load files for spatial plotting 
#======================================================================
#prepare sf layer to be joined to model output for spatial plotting
model_df_spatial <-
  model_df %>%
  dplyr::select(State, County, geometry,group_id,zone) %>%
  distinct()

#prepare mask for us states for plotting
us_States_Sf <- us_states() %>%
  us_spatial_sf %>%
  tibble() %>%
  inner_join(model_pipe %>%
               tibble() %>%
               distinct() %>%
               dplyr::select(State)) %>%
  distinct() %>%
  st_as_sf()

#prepare mask for eu states for plotting
shp_eu <- list.files(dir_eurostat, "shp", full.names = TRUE) %>%
  read_sf() %>%
  filter(LEVL_CODE %in% c(0)) %>%
  dplyr::select(NUTS_ID, CNTR_CODE) %>% 
  rename(County =1,
         State = 2) %>% 
  tibble() %>% 
  inner_join(model_pipe %>%
               tibble() %>%
               distinct() %>%
               dplyr::select(State)) %>%
  distinct() %>%
  st_as_sf()
#======================================================================
#extract model details
#======================================================================
#Durbin-Watson test analyzes the null hypothesis that residuals 
#from the regression are not autocorrelated 
autocorrelation_resid_test <-model_pipe%>%
  dplyr::select(dwtest, group_id) %>%
  unnest(cols = c(dwtest)) %>%
  ungroup() %>%
  mutate(dwtest = ifelse(dwtest < 0.05,1,0)) %>% 
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#the Breusch-Pagan tests for heteroscedasticity in the data. 
#i.e. scatter of errors is different depending on covariate values
#If so the standard errors may be unreliable.
heteroscedasticity_resid_test <-model_pipe%>%
  dplyr::select(bptest, group_id) %>%
  unnest(cols = c(bptest)) %>%
  ungroup() %>%
  mutate(bptest = ifelse(bptest < 0.05,1,0)) %>% 
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#the shapirotest tests for nromality of residuals 
#predictive ability of predictors not the same
#across the full range of the dependent variable
normality_resid_test <-model_pipe%>%
  dplyr::select(shapirotest, group_id) %>%
  unnest(cols = c(shapirotest)) %>%
  ungroup() %>%
  mutate(shapirotest = ifelse(shapirotest < 0.05,1,0)) %>% 
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#the shapirotest tests for nromality of residuals 
#predictive ability of predictors not the same
#across the full range of the dependent variable
power_term_test <-model_pipe%>%
  dplyr::select(resettest, group_id) %>%
  unnest(cols = c(resettest)) %>%
  ungroup() %>%
  mutate(resettest = ifelse(resettest < 0.05,1,0)) %>% 
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#VIF test assesses multicolinearity among variables
vif_test <- model_pipe %>%
  dplyr::select(vif_val, group_id) %>%
  unnest(cols = c(vif_val)) %>%
  ungroup() %>%
  pivot_wider(names_from = var, values_from = vif_val) %>%
  mutate(max_vif = pmax(
    mean_sm_spring,
    mean_sm_summer,
    `mean_sm_spring:mean_sm_summer`
  )) %>% 
  mutate(vif_test = ifelse(max_vif > 5,1,0)) %>%
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#extract model coefficients at county scale and join sf layer
model_coefficients <- model_pipe %>%
  dplyr::select(model_coefs, group_id) %>%
  unnest(cols = c(model_coefs)) %>%
  ungroup() %>%
  inner_join(model_df_spatial) %>%
 # filter(!term %in% c("(Intercept)", "Year")) %>%
  drop_na() %>% 
  st_as_sf()

#extract model rsquare at county scale and join sf layer
model_rsquare <- model_pipe %>%
  dplyr::select(rsq_full_model, group_id) %>%
  unnest(cols = c(rsq_full_model)) %>%
  ungroup() %>%
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()

#extract model rsquare at county scale and join sf layer
model_sig <- model_pipe %>%
  dplyr::select(model_pvalue, group_id) %>%
  unnest(cols = c(model_pvalue)) %>%
  ungroup() %>%
  inner_join(model_df_spatial) %>%
  drop_na() %>% 
  st_as_sf()
#======================================================================
#plot (Check load functions for base script)
#======================================================================
#plot EU coefs
plot_model_coefs(area = "EU", pval_trsh = 0.1)
#plot US coefs
plot_model_coefs(area = "US", pval_trsh = 0.05)
#plot EU rsq  
plot_model_rsq(area = "EU",minimum_rsq_val = 0.2)
#plot US rsq
plot_model_rsq(area = "US",minimum_rsq_val = 0.6)
#======================================================================
#midwest subset
midwest_us <- c("ILLINOIS",
                "INDIANA",
                "IOWA",
                "MICHIGAN",
                "MINNESOTA",
                "OHIO",
                "WISCONSIN",
                "KANSAS",
                "MISSOURI",
                "NEBRASKA",
                "NORTH DAKOTA",
                "SOUTH DAKOTA")

#Get coefficients for model with interactions for a specific state
group_select <- model_coefficients %>%
  filter(State %in% midwest_us) %>%
  tibble() %>%
  dplyr::select(term,estimate,group_id) %>% 
  pivot_wider(names_from = "term", values_from = "estimate") %>%
  filter(mean_sm_spring < 0 &
           mean_sm_summer > 0 &
           `mean_sm_spring:mean_sm_summer` < 0)


#Get coefficients for model with interactions for a specific state
coefs_df <- model_coefficients %>% 
  filter(group_id %in% c(group_select$group_id)) %>% 
 # filter(p.value <= pval_thrsh) %>% 
  dplyr::select(County, State, term, estimate) %>%
  group_by(term) %>%
  filter(State %in% midwest_us) %>% 
  tibble() %>% 
  group_by(term) %>% 
  summarise_at(vars(estimate),~median(., na.rm = TRUE))


#get points montana
obs_df <-model_df %>% 
  ungroup() %>% 
  filter(group_id %in% c(group_select$group_id)) %>% 
  dplyr::select(-geometry,-group_id)%>% 
  na.exclude() %>% 
  unique()

obs_df %>%
  ggplot(aes(x = mean_sm_spring, y = mean_sm_summer, color = crop_yield)) +
  geom_text(aes(label= Year)) +
  scale_color_continuous_diverging("Red-Green")+
  theme_bw()+
  geom_vline(xintercept = 0)+
  geom_hline(yintercept = 0)
  










#quick plotting
#======================================================================
#county scale coefs
county_scale_coefs_plot <-
  county_scale_coefs(
    model_coefs_df = model_coefficients,
    model_pvals_df = model_pvals,
    sf_mask = us_States_Sf
  )


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

#======================================================================
#panel model fit
#======================================================================
panel_model <-lmer(
  crop_yield ~  kdd_summer * kdd_spring + gdd_spring  + gdd_summer + (1 |group_id),
  data = model_df %>%  filter(zone == "EU")
)

panel_model_coefs <-panel_model$coefficients %>% 
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
