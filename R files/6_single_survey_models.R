############################
### Single-survey models ###
############################

# Recreate mesh for scallop survey
scallop_mesh <- make_mesh(filter(lob_mod_dat, survey == "scallop"), c("X", "Y"), mesh = barrier_mesh$mesh)
scallop_barrier_mesh <- add_barrier_mesh(scallop_mesh, st_as_sf(map_fundy), plot = TRUE, proj_scaling = 1000)
# Fit scallop survey model
scallopmod <- sdmTMB(count ~ 0+year_fac+moult+s(logdepth_std)+s(broad_BPI_std), 
                     family = nbinom2(),
                     data = filter(lob_mod_dat, survey == "scallop"), mesh = scallop_barrier_mesh,
                     offset = filter(lob_mod_dat, survey == "scallop")$log_area,
                     spatial = "on", spatiotemporal = "ar1", time = "year", 
                     extra_time = 2020, silent = FALSE)
summary(scallopmod)
# Residual analysis
set.seed(11)
scallop_resids_dat <- mutate(scallopmod$data, spresids = residuals(scallopmod))
plot_resids_spatial(scallop_resids_dat, plot_type = "histogram")
plot_resids_spatial(scallop_resids_dat, plot_type = "qq-plot")
plot_resids_spatial(scallop_resids_dat, plot_type = "map")
scalloppreds <- predict(scallopmod, newdata = filter(grid_yrs, year != 2020), type = "link", 
                        offset = filter(grid_yrs, year != 2020)$log_area, return_tmb_object = TRUE)
# Make the SD and CV for the predictions using simulations
sims_scallop <- predict(scallopmod, newdata = filter(grid_yrs, year != 2020), type = "link", offset = filter(grid_yrs, year != 2020)$log_area, nsim = 200)
scalloppreds$data$sd <- round(apply(exp(sims_scallop), 1, function(x) sd(x)), 2)
scalloppreds$data$cv <- round(apply(exp(sims_scallop), 1, function(x) sd(x) / mean(x)), 2)
# Make plots for each year
plot_map(scalloppreds$data, years = 2005:2010, return = "response")
plot_map(scalloppreds$data, years = 2011:2016, return = "response")
plot_map(scalloppreds$data, years = 2017:2023, return = "response")
plot_map(scalloppreds$data, years = 2005:2010, return = "cv")
plot_map(scalloppreds$data, years = 2011:2016, return = "cv")
plot_map(scalloppreds$data, years = 2017:2022, return = "cv")
plot_map(scalloppreds$data, years = 2023, return = "cv")
plot_map(scalloppreds$data, years = 2005:2010, return = "spatiotemporal")
plot_map(scalloppreds$data, years = 2011:2016, return = "spatiotemporal")
plot_map(scalloppreds$data, years = 2017:2022, return = "spatiotemporal")
plot_map(scalloppreds$data, years = 2023, return = "spatiotemporal")
# Calculate scallop survey index
scallop_index <- get_index(scalloppreds, bias_correct = TRUE,
                           area = 4)
scallop_index_scaled_single <- scallop_index %>% 
  mutate(est = est/scallop_index$est[1],
         lwr = lwr/scallop_index$est[1],
         upr = upr/scallop_index$est[1],
         se_natural = se_natural/scallop_index$est[1])
scallop_index_scaled_single <- rbind(scallop_index_scaled_single[1:15,],
                                     c(2020, NA, NA, NA, NA, NA, NA, "index"),
                                     scallop_index_scaled_single[16:18,]) %>% 
  mutate(year = as.numeric(year), est = as.numeric(est), 
         lwr = as.numeric(lwr), upr = as.numeric(upr), 
         log_est = as.numeric(log_est), se = as.numeric(se),
         se_natural = as.numeric(se_natural))
plot_index(scallop_index_scaled_single)
# Calculate scallop survey centre of gravity
scallop_cog <- get_cog(scalloppreds, format = "wide", bias_correct = FALSE)
scallop_cog1000 <- mutate(scallop_cog, est_x = 1000*est_x, est_y = 1000*est_y,
                          lwr_x = 1000*lwr_x, lwr_y = 1000*lwr_y,
                          upr_x = 1000*upr_x, upr_y = 1000*upr_y)
ggplot(map_fundy)+
  geom_sf()+
  geom_point(data = scallop_cog1000, mapping = aes(est_x, est_y, colour = year)) +
  geom_linerange(data = scallop_cog1000, mapping = aes(est_x, est_y, xmin = lwr_x, xmax = upr_x, colour = year)) +
  geom_linerange(data = scallop_cog1000, mapping = aes(est_x, est_y, ymin = lwr_y, ymax = upr_y, colour = year)) +
  scale_colour_viridis_c(option = "magma", direction = -1)+theme_light()+labs(x = "Longitude", y = "Latitude")+
  labs(colour = "Year")

# Recreate mesh for lobster survey
lobster_mesh <- make_mesh(filter(lob_mod_dat, survey == "lobster"), c("X", "Y"), mesh = barrier_mesh$mesh)
lobster_barrier_mesh <- add_barrier_mesh(lobster_mesh, st_as_sf(map_fundy), plot = TRUE, proj_scaling = 1000)
# Fit lobster survey model
lobstermod <- sdmTMB(count ~ 0+year_fac+moult+s(logdepth_std)+broad_BPI_std, 
                     family = nbinom2(), 
                     data = filter(lob_mod_dat, survey == "lobster"), mesh = lobster_barrier_mesh,
                     offset = filter(lob_mod_dat, survey == "lobster")$log_area,
                     spatial = "on", spatiotemporal = "ar1", time = "year", silent = FALSE)
summary(lobstermod)

# Residual analysis
set.seed(11)
lobster_resids_dat <- mutate(lobstermod$data, spresids = residuals(lobstermod))
plot_resids_spatial(lobster_resids_dat, plot_type = "histogram")
plot_resids_spatial(lobster_resids_dat, plot_type = "qq-plot")
plot_resids_spatial(lobster_resids_dat, plot_type = "map")                                                                                                                                                                                 

lobsterpreds <- predict(lobstermod, newdata = grid_yrs, type = "link", 
                        offset = grid_yrs$log_area, return_tmb_object = TRUE)
# Make the SD and CV for the predictions using simulations
sims_lobster <- predict(lobstermod, newdata = grid_yrs, type = "link", offset = grid_yrs$log_area, nsim = 200)
lobsterpreds$data$sd <- round(apply(exp(sims_lobster), 1, function(x) sd(x)), 2)
lobsterpreds$data$cv <- round(apply(exp(sims_lobster), 1, function(x) sd(x) / mean(x)), 2)
# Make plots for every year
plot_map(lobsterpreds$data, years = 2005:2010, return = "response")
plot_map(lobsterpreds$data, years = 2011:2016, return = "response")
plot_map(lobsterpreds$data, years = 2017:2022, return = "response")
plot_map(lobsterpreds$data, years = 2023, return = "response")
plot_map(lobsterpreds$data, years = 2005:2010, return = "cv")
plot_map(lobsterpreds$data, years = 2011:2016, return = "cv")
plot_map(lobsterpreds$data, years = 2017:2022, return = "cv")
plot_map(lobsterpreds$data, years = 2023, return = "cv")
plot_map(lobsterpreds$data, years = 2005:2010, return = "spatiotemporal")
plot_map(lobsterpreds$data, years = 2011:2016, return = "spatiotemporal")
plot_map(lobsterpreds$data, years = 2017:2022, return = "spatiotemporal")
plot_map(lobsterpreds$data, years = 2023, return = "spatiotemporal")
# Calculate lobster survey index
lobster_index <- get_index(lobsterpreds, bias_correct = TRUE,
                           area = 4)
lobster_index_scaled_single <- lobster_index %>% 
  mutate(est = est/lobster_index$est[1],
         lwr = lwr/lobster_index$est[1],
         upr = upr/lobster_index$est[1],
         se_natural = se_natural/lobster_index$est[1])
plot_index(lobster_index_scaled_single)
# Calculate lobster survey centre of gravity
lobster_cog <- get_cog(lobsterpreds, format = "wide", bias_correct = FALSE)
lobster_cog1000 <- mutate(lobster_cog, est_x = 1000*est_x, est_y = 1000*est_y,
                          lwr_x = 1000*lwr_x, lwr_y = 1000*lwr_y,
                          upr_x = 1000*upr_x, upr_y = 1000*upr_y)
ggplot(map_fundy)+
  geom_sf()+
  geom_point(data = lobster_cog1000, mapping = aes(est_x, est_y, colour = year)) +
  geom_linerange(data = lobster_cog1000, mapping = aes(est_x, est_y, xmin = lwr_x, xmax = upr_x, colour = year)) +
  geom_linerange(data = lobster_cog1000, mapping = aes(est_x, est_y, ymin = lwr_y, ymax = upr_y, colour = year)) +
  scale_colour_viridis_c(option = "magma", direction = -1)+theme_light()+labs(x = "Longitude", y = "Latitude")+
  labs(colour = "Year")
