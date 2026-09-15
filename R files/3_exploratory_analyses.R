############################
### Exploratory analyses ###
############################

# Map the LFAs
lfas <- readRDS("LFAPolysSF.rds")%>% 
  filter(LFA %in% 34:38) %>% 
  mutate(area = as.numeric(st_area(geometry))/1000)
islands <- lfas %>% 
  filter(area <= 200000) %>% 
  st_union()
lfas <- st_difference(lfas, islands) %>% 
  filter(area > 200000)
locations <- data.frame(
  place = c("Yarmouth", "St. John", "Lobster Bay", "St. Mary's Bay", 
            "Grand Manan Island", "Cape Sable Island"),
  X = c(249598, 260222, 265972, 256562, 191394, 288945),
  Y = c(4858144, 5017486, 4839875, 4926698, 4947340, 4814107)
) 
p_lfa <- ggplot(map_fundy)+geom_sf()+
  geom_sf(data = lfas, aes(fill = LFA), alpha = 0.3)+
  scale_fill_discrete(guide = "none")+
  geom_sf_text(data = lfas, aes(label = LFA))+
  geom_point(data = locations, aes(X, Y), colour = "blue")+
  ggrepel::geom_label_repel(data = locations, aes(X, Y, label = place),
                            min.segment.length = 0,
                            size = 3,
                            point.padding = 1,
                            segment.color = 'grey50')+
  labs(x = "Longitude", y = "Latitude")+
  theme_light()+
  coord_sf(xlim = c(130218, 461283),  ylim = c(4705414, 5094169), expand = F)
# Map the defined subareas
areas_sf_map <- mutate(areas_sf, 
                       Area = factor(name, levels = c("SMB", "OSMB", "JR", "GLB", "OUT"),
                                     labels = c("SMB", "Outer SMB", "JR", "GLB", "34 Outside")))
p_areas <- ggplot(map_fundy)+geom_sf()+
  geom_sf(data = areas_sf_map, aes(fill = Area), alpha = 0.3)+
  geom_sf_label(data = areas_sf_map, aes(label = Area))+
  scale_fill_discrete(guide = "none")+
  labs(x = "Longitude", y = "Latitude")+
  theme_light()+
  coord_sf(xlim = c(130218, 341283),  ylim = c(4705414, 4994169), expand = F)
cowplot::plot_grid(p_lfa, p_areas, labels = c('(a)', '(b)'), label_size = 12, nrow = 1,
                   align = "v", 
                   hjust = c(-2, -2))

# Map of tows by survey for selected representative years
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year %in% rep_years), 
             aes(x = X * 1000, y = Y * 1000, col = fct_rev(survey)), 
             alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  facet_wrap(~year)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "", y = "", col = "Survey")+
  theme(axis.text = element_text(size = rel(0.5)))+
  coord_sf(xlim = c(130218, 461283),  ylim = c(4705414, 5094169), expand = F)+
  scale_x_continuous(breaks = c(-67, -65))+
  scale_y_continuous(breaks = c(43, 44, 45))

# Map of tows by survey 2005-10
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year %in% 2005:2010), 
             aes(x = X * 1000, y = Y * 1000, col = fct_rev(survey)), 
             alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  facet_wrap(~year)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "Longitude", y = "Latitude", col = "Survey")+
  ggtitle("Map of survey tows")
# Map of tows by survey 2011-16
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year %in% 2011:2016), 
             aes(x = X * 1000, y = Y * 1000, col = fct_rev(survey)), 
             alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  facet_wrap(~year)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "Longitude", y = "Latitude", col = "Survey")+
  ggtitle("Map of survey tows")
# Map of tows by survey 2017-22
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year %in% 2017:2022), 
             aes(x = X * 1000, y = Y * 1000, col = fct_rev(survey)), 
             alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  facet_wrap(~year)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "Longitude", y = "Latitude", col = "Survey")+
  ggtitle("Map of survey tows")
# Map of tows by survey 2023
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year == 2023), 
             aes(x = X * 1000, y = Y * 1000, col = fct_rev(survey)), 
             alpha = 0.4, size = 0.01) +
  scale_colour_viridis_d(option = "D", direction = -1)+
  theme_light() +
  guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(x = "Longitude", y = "Latitude", col = "Survey")+
  ggtitle("Map of survey tows (2023)")

# Histogram of tows by year and  month
ggplot(lob_dat_covars, aes(year, fill = month))+
  geom_histogram(alpha = 0.4, bins = 19) +
  labs(x = "Year", y = "Tows", fill = "Month") +
  scale_fill_discrete(labels = c("May", "June", "July", "August", "September", "October"))+
  facet_wrap(~survey)+
  theme_light()+
  ggtitle("Survey tows by year and month")

# Histogram of counts, by survey
ggplot(lob_dat_covars, aes(count, fill = fct_rev(survey)))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  scale_fill_viridis_d(option = "D", direction = -1)+
  scale_y_continuous(transform = "log1p", breaks = 10^(0:4))+
  labs(x = "Count of lobster recruits", y = "Frequency (log scale)", fill = "Gear type") +
  theme_light()

# Histogram of swept areas, by gear type
ggplot(lob_dat_covars, aes(area_km2, fill = fct_rev(survey)))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  scale_fill_viridis_d(option = "D", direction = -1)+
  scale_y_continuous(transform = "log1p", breaks = 10^(0:4))+
  labs(x = "Swept area (km2)", y = "Frequency (log scale)", fill = "Gear type") +
  theme_light()

# Histogram of density, by survey
ggplot(lob_dat_covars, aes(count_per_km2, fill = survey))+
  geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
  scale_fill_viridis_d(option = "D", direction = 1)+
  scale_y_continuous(transform = "log1p", breaks = 10^(0:4))+
  labs(x = "Lobsters per km2", y = "Frequency (log scale)", 
       fill = "Survey", title = "Lobster recruit density") +
  theme_light()

# Map of counts
ggplot(map_fundy) + geom_sf() +
  geom_point(data = lob_dat_covars, aes(x = X * 1000, y = Y * 1000, col = count), alpha = 0.4, size = 0.01) +
  scale_colour_viridis_c(option = "B", direction = -1, transform = "log1p", breaks = 10^(0:3))+
  theme_light()+
  labs(x = "Longitude", y = "Latitude", col = "Lobster recruit count")

# Map of density
ggplot(map_fundy) + geom_sf() +
  geom_point(data = lob_dat_covars, aes(x = X * 1000, y = Y * 1000, col = count_per_km2), alpha = 0.4, size = 0.01) +
  scale_colour_viridis_c(option = "B", direction = -1, transform = "log1p", breaks = 10^(0:4))+
  theme_light()+
  labs(x = "Longitude", y = "Latitude", col = "Lobster recruit density")

# Map of density in representative years
ggplot(map_fundy) + geom_sf() +
  geom_point(data = filter(lob_dat_covars, year %in% rep_years), 
             aes(x = X * 1000, y = Y * 1000, col = count_per_km2), 
             alpha = 0.4, size = 0.01) +
  facet_wrap(~year)+
  scale_colour_viridis_c(option = "B", direction = -1, transform = "log1p", breaks = 10^(0:4))+
  theme_light()+
  labs(x = "Longitude", y = "Latitude", col = "Lobsters per km2",
       title = "Lobster recruit density")+
  theme(axis.text = element_text(size = rel(0.5)))

############# Initial variable selection #######
# Make boxplots of temperatures in each month, and their average
dat_long <- lob_dat_covars %>% 
  dplyr::select(month5, month6, month7, month8, month9, month10, temperature) %>% 
  pivot_longer(everything(), 
               names_to = "month", 
               values_to = "temp") %>% 
  mutate(month = factor(month, levels = c("month5", "month6", "month7", "month8", "month9", "month10", "temperature")))
ggplot(dat_long, aes(x = month, y = temp))+
  geom_boxplot()+
  scale_x_discrete(labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Mean"))+
  theme_light()+ 
  labs(x = "Month", y = "Temperature")+
  ggtitle("Boxplots of temperatures each month\nand averaged over the months")

# Make scatterplots, get pairwise Spearman correlations, and calculate VIFs for each covariate group
# Curvature
ggpairs(lob_dat_covars[, 34:39],
        upper = list(continuous = wrap("cor", method = "spearman")))+ggtitle("Curvature")
vif(lob_dat_covars[, 34:39])
# Conclusion: remove mean, max, and min
# Relative position
ggpairs(lob_dat_covars[, 28:30],
        upper = list(continuous = wrap("cor", method = "spearman")))+ggtitle("Relative position")
vif(lob_dat_covars[, 28:30])
# Conclusion: remove DMV
# Roughness
ggpairs(lob_dat_covars[, 40:43],
        upper = list(continuous = wrap("cor", method = "spearman")))+ggtitle("Roughness")
vif(lob_dat_covars[, 40:43])
# Conclusion: These are all highly correlated; keep only adj_SD
# Oceanography
ggpairs(lob_dat_covars[, 44:49],
        upper = list(continuous = wrap("cor", method = "spearman")))+ggtitle("Oceanography")
vif(lob_dat_covars[, 44:49])
# Conclusion: keep seabed disturbance
# Make correlation matrix and calculate VIFs for the chosen covariates
cor_mat <- cor(lob_dat_covars[, c(27:29, 31:36, 42, 48, 51)], method = "spearman", use = "pairwise.complete")
psych::corPlot(cor_mat, main = "All selected covariates")
vif(lob_dat_covars[, c(27:29, 31:36, 42, 48, 51)])

# Make maps of each covariate across the study area
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = bathymetry))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Depth (m)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = broad_BPI))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Broad BPI (m)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = fine_BPI))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Fine BPI (m)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = slope))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Slope (degrees)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = eastness))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Eastness", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = northness))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Northness", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = curvature_profile))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Curvature - profile", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = curvature_planform))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Curvature - planform", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = curvature_twisting))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Curvature - twisting", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = adj_SD))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Adjusted SD (m)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = seabed_dist))+
  scale_fill_continuous()+
  theme_light()+
  labs(fill = "Seabed disturbance index", x = "Longitude", y = "Latitude")
# Mape the temperature for every year
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year %in% 2005:2010), aes(xmetres, ymetres, fill = month_temps))+
  scale_fill_continuous()+
  facet_wrap(~year)+
  theme_light()+
  labs(fill = "Temperature (C)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year %in% 2011:2016), aes(xmetres, ymetres, fill = month_temps))+
  scale_fill_continuous()+
  facet_wrap(~year)+
  theme_light()+
  labs(fill = "Temperature (C)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year %in% 2017:2022), aes(xmetres, ymetres, fill = month_temps))+
  scale_fill_continuous()+
  facet_wrap(~year)+
  theme_light()+
  labs(fill = "Temperature (C)", x = "Longitude", y = "Latitude")
ggplot(map_fundy)+
  geom_sf()+
  geom_tile(data = filter(grid_yrs_all_covars, year == 2023), aes(xmetres, ymetres, fill = month_temps))+
  scale_fill_continuous()+
  facet_wrap(~year)+
  theme_light()+
  labs(fill = "Temperature (C)", x = "Longitude", y = "Latitude")
# Standardize covariates to prepare for modelling
lob_mgcv_dat <- lob_dat_covars %>% 
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>% 
  filter(!if_any(bathymetry:temperature, is.na) & bathymetry < 0) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(logdepth))/sd(logdepth),
         broad_BPI_std = (broad_BPI-mean(broad_BPI))/sd(broad_BPI),
         fine_BPI_std = (fine_BPI-mean(fine_BPI))/sd(fine_BPI),
         sqrtslope_std = (sqrtslope-mean(sqrtslope))/sd(sqrtslope),
         curv_prof_std = curvature_profile/sd(curvature_profile),
         curv_plan_std = curvature_planform/sd(curvature_planform),
         curv_twist_std = curvature_twisting/sd(curvature_twisting),
         logadjSD_std = (logadjSD-mean(logadjSD))/sd(logadjSD),
         temp_std = (month_temps-mean(month_temps))/sd(month_temps))
# Fit GAM in mgcv
mgcvmod <- mgcv::gam(count ~ year_fac+survey+moult+s(logdepth_std, k = 3)+s(broad_BPI_std, k = 3)+
                       s(fine_BPI_std, k = 3)+s(sqrtslope_std, k = 3)+
                       s(eastness, k = 3)+s(northness, k = 3)+s(curv_prof_std, k = 3)+s(curv_plan_std, k = 3)+
                       s(curv_twist_std, k = 3)+s(logadjSD_std, k = 3)+s(seabed_dist, k = 3)+s(temp_std, k = 3), family = mgcv::nb(), 
                     data = lob_mgcv_dat, offset = lob_mgcv_dat$logarea, select = TRUE)
summary(mgcvmod)
plot(mgcvmod, rug = TRUE, pages = 1)
# Fit second gam with fewer covariates
mgcvmod2 <- mgcv::gam(count ~ year_fac+survey+moult+s(logdepth_std, k = 3)+
                        s(broad_BPI_std, k=3)+
                        s(logadjSD_std)+s(seabed_dist, k = 3)+s(temp_std, k=3), family = mgcv::nb(), 
                      data = lob_mgcv_dat, offset = lob_mgcv_dat$logarea)
summary(mgcvmod2)
plot(mgcvmod2, rug = TRUE, pages = 1)