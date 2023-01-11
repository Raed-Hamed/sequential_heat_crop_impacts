#load functions
us_spatial_sf<-  function(x) {
  #x is either us_counties or us_states from usaboundaries package  
  x %>% 
    rename(County = name,
           State  = state_name) %>% 
    dplyr::select(County,State,geometry) %>% 
    mutate_at(vars(State,County), ~toupper(.))
  
}
#======================================================================
#load temperature data
load_cpc_t <-function(filename, varname, shapefile){
  
  list.files(dir_data, filename, full.names = TRUE) %>%
    brick(varname = varname) %>%
    crop(shapefile) %>%
    mask(shapefile) %>%
    rasterToPoints() %>%
    data.frame() %>%
    tibble() %>%
    pivot_longer(X1980:X2021) %>%
    group_by(x, y) %>%
    mutate(Year = 1980:2021) %>%
    dplyr::select(-name) %>%
    ungroup()
  
}
#======================================================================
#print lm model pvalue
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}
#======================================================================
#leave-one-out function
fun_loo_per_grid <- function(x,model){
  predictions <-rep(NA, length(x$crop_yield))
  List_year <-unique(x$Year)
  
  for (i in 1:length(List_year)){
    
    Training_i <- x[x$Year!=List_year[i],]
    Test_i <- x[x$Year==List_year[i],]
    Mod_i <- lm(formula(model),data =Training_i)
    Y_lm_i <- predict(Mod_i, newdata = Test_i,type="response")
    predictions[x$Year==List_year[i]] <-Y_lm_i
    
  }
  return(predictions) 
}
#======================================================================
#r-squared function
rsq <- function (x, y) cor(x, y) ^ 2
#======================================================================
#plot coefficients
plot_model_coefs <- function(area, pval_trsh) {
  p1 <- model_coefficients %>%
    dplyr::select(estimate, term, zone) %>%
    filter(zone == area) %>%
    ggplot() +
    {
      if (area == "US")
        geom_sf(data = us_States_Sf,
                fill = "lightgray",
                color = "transparent")
    } +
    {
      if (area == "EU")
        geom_sf(data = shp_eu,
                fill = "lightgray",
                color = "transparent")
    } +
    geom_sf(aes(fill  = estimate), color = "transparent") +
    geom_sf(
      data = model_coefficients %>% filter(p.value <= pval_trsh) %>%
        filter(zone == area),
      fill = "transparent"
    ) +
    
    facet_wrap( ~ term) +
   # scale_fill_continuous_diverging("Blue-Red 3", limit = c(-1, 1)) +
    scale_fill_continuous_diverging("Blue-Red 3", limit = c(-2000, 2000)) +
    

    {
      if (area == "EU")
        geom_sf(data = shp_eu,
                fill = "transparent",
                color = "black")
    } +
    {
      if (area == "EU")
        coord_sf(c(-20, 40), c(30, 60))
    } +
    {
      if (area == "US")
        geom_sf(data = us_States_Sf,
                fill = "transparent",
                color = "black")
    } +
    theme_bw(base_size = 15)
  
  return(p1)
}
#======================================================================
#plot R-square
plot_model_rsq <- function(area,minimum_rsq_val){
  
  p1 <-model_rsquare %>%
    dplyr::select(rsq_full_model,zone) %>%
    filter(zone == area) %>%
    filter(rsq_full_model >= minimum_rsq_val) %>% 
    ggplot() +
    {
      if (area == "US")
        geom_sf(data = us_States_Sf,
                fill = "lightgray",
                color = "transparent")
    } +
    {
      if (area == "EU")
        geom_sf(data = shp_eu,
                fill = "lightgray",
                color = "transparent")
    } +
    geom_sf(aes(fill = rsq_full_model), color = "transparent") +
    geom_sf(data = model_sig %>% filter(model_pvalue <= 0.05) %>%
              filter(zone == area), fill = "transparent") +
    scale_fill_continuous_sequential("PuRd", limit = c(0, 1)) +
    {
      if (area == "EU")
        geom_sf(data = shp_eu,
                fill = "transparent",
                color = "black")
    } +
    {
      if (area == "EU")
        coord_sf(c(-20, 40), c(30, 60))
    } +
    {
      if (area == "US")
        geom_sf(data = us_States_Sf,
                fill = "transparent",
                color = "black")
    } +
    theme_bw(base_size = 15)
  
  return(p1)
}
#======================================================================
library(RColorBrewer)
breaks_plot <-seq(-0.5,0.5,0.1)
n <- breaks_plot %>%  length()
getPalette <- colorRampPalette(brewer.pal(11, "PiYG"))
colors_plot <- getPalette(n) #%>%  rev()
limit <-c(-0.5,0.5)

#plot coefficients at county scale
county_scale_coefs <- function(model_coefs_df, model_pvals_df, sf_mask){
  
  ggplot() +
    geom_sf(data = sf_mask, fill = "darkgrey") +
    geom_sf(data = model_coefs_df,
            aes(fill = estimate),
            color = "transparent") +
    geom_sf(data = model_pvals_df, size = 0.1) +
    facet_wrap( ~ term) +
    geom_sf(
      data = model_coefs_df %>%
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
      breaks = seq(-0.8, 0.8, 0.1),
      limits = c(-0.8,  0.8)
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
  geom_sf(data = sf_mask, fill = "transparent") +
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
  
  
}
  
  
  

#plot dependence plot based on averaging county scale coefficients per state
plot_dependence_per_state <-function(model_coefs,observation_df,state_name,pval_thrsh){
  
  # state_name <- "MONTANA"
  # pval_thrsh <- 0.1
  
  sign_counties <-model_coefs %>% 
    filter(p.value <= pval_thrsh) %>% 
    tibble() %>% 
    .$County
  
  #Get coefficients for model with interactions for a specific state
  coefs_df <- model_coefs %>% 
    filter(p.value <= pval_thrsh) %>% 
    dplyr::select(County, State, term, estimate) %>%
    group_by(term) %>%
    filter(State == state_name) %>% 
    tibble() %>% 
    group_by(term) %>% 
    summarise_at(vars(estimate),~mean(., na.rm = TRUE))
  
  #Get coefficients with confidence bounds
  coefs_df_confidence_interval <- model_coefs%>%
    tibble() %>% 
    filter(p.value <= pval_thrsh) %>% 
    dplyr::select(County, State, term, estimate, conf.low, conf.high) %>%
    group_by(term) %>%
    filter(State == state_name) %>% 
    summarize_at(vars(-County, -State),  ~ (mean(., na.rm = TRUE)))
  
  #get points montana
  obs_df <-observation_df %>% 
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
