#load climate df
climate_sf <-
  list.files(dir_model, "climate_sf", full.names = TRUE)  %>%
  readRDS() %>% 
  drop_na()


test <- climate_sf %>%
  group_by(x, y) %>%
  nest() %>%
  mutate(cor_spring = map(data, ~ cor.test(.$mean_sm_spring, .$kdd_spring)$estimate)) %>%
  mutate(pvl_spring = map(data, ~ cor.test(.$mean_sm_spring, .$kdd_spring)$p.value)) %>%
  mutate(cor_summer = map(data, ~ cor.test(.$mean_sm_summer, .$kdd_summer)$estimate)) %>%
  mutate(pvl_summer = map(data, ~ cor.test(.$mean_sm_summer, .$kdd_summer)$p.value)) 
  
cor_out <-test %>% 
  dplyr::select(-data) %>% 
  unnest(c(cor_spring,pvl_spring,cor_summer,pvl_summer)) %>% 
  ungroup() %>% 
  pivot_longer(c(-x,-y)) %>% 
  separate(name,c("type","season")) %>% 
  pivot_wider(names_from = type, values_from = value)

ggplot() +
  geom_tile(data = cor_out, aes(x = x, y = y, fill = cor)) +
  geom_point(data = cor_out %>% filter(pvl <= 0.05), (aes(x = x, y = y)), size = 0.00000001) +
  scale_fill_continuous_sequential("Blues 3", rev = TRUE, c(-1,0))+
  facet_grid(~ season) +
  theme_bw()

test2 <- climate_sf %>%
  group_by(x, y) %>%
  nest() %>%
  mutate(cor_kdd = map(data, ~ cor.test(.$kdd_summer, .$kdd_spring)$estimate)) %>%
  mutate(pvl_kdd = map(data, ~ cor.test(.$kdd_summer, .$kdd_spring)$p.value)) %>%
  mutate(cor_gdd = map(data, ~ cor.test(.$gdd_summer, .$gdd_spring)$estimate)) %>%
  mutate(pvl_gdd = map(data, ~ cor.test(.$gdd_summer, .$gdd_spring)$p.value)) 

cor_out2 <-test2 %>% 
  dplyr::select(-data) %>% 
  unnest(c(cor_kdd,pvl_kdd,cor_gdd,pvl_gdd)) %>% 
  ungroup() %>% 
  pivot_longer(c(-x,-y)) %>% 
  separate(name,c("type","var")) %>% 
  pivot_wider(names_from = type, values_from = value)

ggplot() +
  geom_tile(data = cor_out2, aes(x = x, y = y, fill = cor)) +
  geom_point(data = cor_out2 %>% filter(pvl <= 0.05), (aes(x = x, y = y)), size = 0.00000001) +
  scale_fill_continuous_sequential("Reds 3", limit = c(0,1))+
  facet_grid(~ var) +
  theme_bw()
# +
#   geom_sf(
#     data = model_df_spatial %>%  filter(zone == "EU") %>% dplyr::select(geometry) %>%  st_as_sf(),
#     fill = "transparent",
#     size = 0.00000001
#   )
