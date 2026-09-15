############################################################
### Fit an integrated sdmTMB model and conduct analyses  ###
############################################################

# Make the mesh
mesh <- make_mesh(lob_mod_dat, xy_cols = c("X", "Y"), cutoff = 5,
                  max.edge = c(40, 80))
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, st_as_sf(map_fundy), proj_scaling = 1000, plot = TRUE)
# Fit the model (this is the one discussed in the paper)
stmod <- sdmTMB(count ~ 0+year_fac+survey+moult+s(logdepth_std)+s(broad_BPI_std),
                family = nbinom2(),
                spatial_varying = ~survey,
                data = lob_mod_dat, mesh = barrier_mesh,
                offset = lob_mod_dat$log_area,
                spatial = "on", spatiotemporal = "ar1",
                share_range = TRUE, time = "year")
summary(stmod)
# Verify that the broad BPI covariate is important
stmodnobpi <- sdmTMB(count ~ 0+year_fac+survey+moult+s(logdepth_std),
                     family = nbinom2(),
                     spatial_varying = ~survey,
                     data = lob_mod_dat, mesh = barrier_mesh,
                     offset = lob_mod_dat$log_area,
                     spatial = "on", spatiotemporal = "ar1",
                     share_range = TRUE, time = "year")
summary(stmodbpi)
# Null model for calculating deviance explained
nullmod <- sdmTMB(count ~ 1, family = nbinom2(), 
                  data = lob_mod_dat, mesh = barrier_mesh,
                  offset = lob_mod_dat$log_area,
                  spatial = "off", spatiotemporal = "off",
                  time = "year",
                  control = sdmTMBcontrol(
                    start = list(ln_phi = stmod$parlist$ln_phi), # start it at desired value
                    map = list(ln_phi = factor(NA)) # don't estimate it
                  )) 
1-deviance(stmod)/deviance(nullmod) # 0.8517155
1-deviance(stmodnobpi)/deviance(nullmod) # 0.8511648

AIC(nullmod, stmod, stmodbpi)
cAIC(stmod) # 27957.97
cAIC(stmodbpi) # 27929.73

################### Results and residuals ###########
# Print model parameter estimates and standard errors
summary(stmod)
# Display conditional effects of covariates (on the link scale)
visreg::visreg(stmod, xvar = "logdepth_std", gg = TRUE)+theme_light()+
  labs(title = "Conditional effect of standardized log depth (on link scale)",
       x = "Standardized log depth", y = "Link (log) scale")
visreg::visreg(stmod, xvar = "broad_BPI_std", gg = TRUE)+theme_light()+
  labs(title = "Conditional effect of standardized broad BPI (on link scale)",
       x = "Standardized broad BPI", y = "Link (log) scale")

# Get randomized quantile residuals
set.seed(11)
nb_resids_dat <- mutate(stmod$data, spresids = residuals(stmod))
# Plot the residuals in different ways
plot_resids_spatial(nb_resids_dat, plot_type = "histogram")
plot_resids_spatial(nb_resids_dat, plot_type = "qq-plot")
plot_resids_spatial(nb_resids_dat, plot_type = "map")
# Make maps of residuals for each year
plot_resids_spatial(filter(nb_resids_dat, year %in% rep_years), plot_type = "map")+
  facet_wrap(~year)+
  theme(axis.text = element_text(size = rel(0.5)))
plot_resids_spatial(filter(nb_resids_dat, year %in% 2005:2010), plot_type = "map")+
  facet_wrap(~year)+
  theme(axis.text = element_text(size = rel(0.5)))
plot_resids_spatial(filter(nb_resids_dat, year %in% 2011:2016), plot_type = "map")+
  facet_wrap(~year)+
  theme(axis.text = element_text(size = rel(0.5)))
plot_resids_spatial(filter(nb_resids_dat, year %in% 2017:2022), plot_type = "map")+
  facet_wrap(~year)+
  theme(axis.text = element_text(size = rel(0.5)))
plot_resids_spatial(filter(nb_resids_dat, year == 2023), plot_type = "map")+
  facet_wrap(~year)+
  theme(axis.text = element_text(size = rel(0.5)))

################### Predictions and derived quantities ############
# Make density predictions on a grid
stpreds <- predict(stmod, newdata = grid_yrs, type = "link", 
                   offset = grid_yrs$log_area, return_tmb_object = TRUE)
# Make the SD and CV for the predictions using simulations
sims <- predict(stmod, newdata = grid_yrs, type = "link", offset = grid_yrs$log_area, nsim = 200)
stpreds$data$sd <- round(apply(exp(sims), 1, function(x) sd(x)), 2)
stpreds$data$cv <- round(apply(exp(sims), 1, function(x) sd(x) / mean(x)), 2)

# Map the predictions (response scale) and spatial random field (link scale)
plot_map(stpreds$data, years = rep_years, return = "response")
plot_map(stpreds$data, years = rep_years, return = "responselog")
plot_map(stpreds$data, years = 2023, return = "spatial")
plot_map(stpreds$data, years = rep_years, return = "spatiotemporal")
plot_map(stpreds$data, years = 2023, return = "spatial-varying")
# Map the SD and CV
plot_map(stpreds$data, years = rep_years, return = "sd")
plot_map(stpreds$data, years = rep_years, return = "cv")

# Plots for all years
plot_map(stpreds$data, years = 2005:2010, return = "response")
plot_map(stpreds$data, years = 2011:2016, return = "response")
plot_map(stpreds$data, years = 2017:2022, return = "response")
plot_map(stpreds$data, years = 2023, return = "response")

plot_map(stpreds$data, years = 2005:2010, return = "cv")
plot_map(stpreds$data, years = 2011:2016, return = "cv")
plot_map(stpreds$data, years = 2017:2022, return = "cv")
plot_map(stpreds$data, years = 2023, return = "cv")

plot_map(stpreds$data, years = 2005:2010, return = "spatiotemporal")
plot_map(stpreds$data, years = 2011:2016, return = "spatiotemporal")
plot_map(stpreds$data, years = 2017:2022, return = "spatiotemporal")
plot_map(stpreds$data, years = 2023, return = "spatiotemporal")

# Create abundance index
st_index <- get_index(stpreds, bias_correct = TRUE,
                      area = 4)
st_index_scaled <- st_index %>% 
  mutate(est = est/st_index$est[1],
         lwr = lwr/st_index$est[1],
         upr = upr/st_index$est[1])
# Plot the index
p_index <- plot_index(st_index_scaled)+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous(limits = c(0, 4))
p_index
# Create centre of gravity
cog <- get_cog(stpreds, format = "wide", bias_correct = FALSE)
cog1000 <- mutate(cog, est_x = 1000*est_x, est_y = 1000*est_y,
                  lwr_x = 1000*lwr_x, lwr_y = 1000*lwr_y,
                  upr_x = 1000*upr_x, upr_y = 1000*upr_y)
# Plot the centre of gravity (on different scales)
cog_out <- ggplot(map_fundy)+
  geom_sf()+
  geom_point(data = cog1000, mapping = aes(est_x, est_y, colour = year)) +
  geom_linerange(data = cog1000, mapping = aes(est_x, est_y, xmin = lwr_x, xmax = upr_x, colour = year)) +
  geom_linerange(data = cog1000, mapping = aes(est_x, est_y, ymin = lwr_y, ymax = upr_y, colour = year)) +
  scale_colour_viridis_c(option = "magma", direction = -1)+
  theme_light()+labs(x = "", y = "")+
  labs(colour = "Year")+
  theme(legend.position = "none")+
  coord_sf(xlim = c(130218, 461283),  ylim = c(4705414, 5094169), expand = F)
cog_in <- ggplot(map_fundy)+
  geom_sf()+
  geom_point(data = cog1000, mapping = aes(est_x, est_y, colour = year)) +
  geom_linerange(data = cog1000, mapping = aes(est_x, est_y, xmin = lwr_x, xmax = upr_x, colour = year)) +
  geom_linerange(data = cog1000, mapping = aes(est_x, est_y, ymin = lwr_y, ymax = upr_y, colour = year)) +
  scale_colour_viridis_c(option = "magma", direction = -1)+theme_light()+labs(x = "", y = "")+
  labs(colour = "Year")+
  coord_sf(xlim = c(217035, 301254),  ylim = c(4822294, 4930449), expand = F)
cowplot::plot_grid(cog_out, cog_in, labels = c('(a)', '(b)'), label_size = 12, nrow = 1,
                   align = "hv")

############### Index in each LFA ##########
# Calculate recruit index for each LFA
preds34 <- predict(stmod, newdata = grid34_yrs, type = "link", 
                   offset = grid34_yrs$log_area, return_tmb_object = TRUE)
index34 <- get_index(preds34, bias_correct = TRUE,
                     area = 4)

preds35 <- predict(stmod, newdata = grid35_yrs, type = "link", 
                   offset = grid35_yrs$log_area, return_tmb_object = TRUE)
index35 <- get_index(preds35, bias_correct = TRUE,
                     area = 4)

preds36 <- predict(stmod, newdata = grid36_yrs, type = "link", 
                   offset = grid36_yrs$log_area, return_tmb_object = TRUE)
index36 <- get_index(preds36, bias_correct = TRUE,
                     area = 4)

preds38 <- predict(stmod, newdata = grid38_yrs, type = "link", 
                   offset = grid38_yrs$log_area, return_tmb_object = TRUE)
index38 <- get_index(preds38, bias_correct = TRUE,
                     area = 4)
# Plot all indices together
index_all_LFAs <- cbind(rbind(index34, index35, index36, index38),
                        est_total = rep(st_index$est, 4),
                        LFA = factor(rep(c("34", "35", "36", "38"), each = 19))) %>% 
  mutate(est_scaled = est/st_index$est[1],
         lwr_scaled = lwr/st_index$est[1],
         upr_scaled = upr/st_index$est[1],
         est_pct_total = est/est_total)
p_index_lfa <- ggplot(index_all_LFAs, aes(year, est_scaled)) + 
  geom_line(aes(colour = LFA)) +
  geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = LFA), alpha = 0.1) +
  labs(x = "Year", y = "Index value")+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_colour_manual(values = c("blue","black","orange","firebrick2"), guide = "none")+
  scale_fill_manual(values = c("blue","black","orange","firebrick2"), guide = "none")+
  theme_light()
# Plot LFA indices as percent of total study area index
p_pct_lfa <- ggplot(index_all_LFAs, aes(year, est_pct_total, fill = LFA)) + 
  geom_col(position = "fill")+
  labs(x = "Year", y = "Proportion of total estimate", fill = "Area")+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous()+
  scale_fill_manual(values = c("blue","black","orange","firebrick2"))+
  theme_light()
cowplot::plot_grid(p_index_lfa, p_pct_lfa,labels = c('(a)', '(b)'), label_size = 12, nrow = 1,
                   rel_widths = c(1.75, 2))

############### Index in each subarea ##########
# Calculate indices for each subarea
predsSMB <- predict(stmod, newdata = gridSMB_yrs, type = "link", 
                    offset = gridSMB_yrs$log_area, return_tmb_object = TRUE)
indexSMB <- get_index(predsSMB, bias_correct = TRUE,
                      area = 4)

predsOSMB <- predict(stmod, newdata = gridOSMB_yrs, type = "link", 
                     offset = gridOSMB_yrs$log_area, return_tmb_object = TRUE)
indexOSMB <- get_index(predsOSMB, bias_correct = TRUE,
                       area = 4)

predsJR <- predict(stmod, newdata = gridJR_yrs, type = "link", 
                   offset = gridJR_yrs$log_area, return_tmb_object = TRUE)
indexJR <- get_index(predsJR, bias_correct = TRUE,
                     area = 4)

predsGLB <- predict(stmod, newdata = gridGLB_yrs, type = "link", 
                    offset = gridGLB_yrs$log_area, return_tmb_object = TRUE)
indexGLB <- get_index(predsGLB, bias_correct = TRUE,
                      area = 4)

predsOUT <- predict(stmod, newdata = gridOUT_yrs, type = "link", 
                    offset = gridOUT_yrs$log_area, return_tmb_object = TRUE)
indexOUT <- get_index(predsOUT, bias_correct = TRUE,
                      area = 4)
# Plot all of these indices together
index_all_areas <- cbind(rbind(indexSMB, indexOSMB, indexJR, indexGLB, indexOUT),
                         lfa34 = rep(index34$est, 5),
                         region = factor(rep(c("SMB", "OSMB", "JR", "GLB", "OUT"), each = 19),
                                         levels = c("SMB", "OSMB", "JR", "GLB", "OUT"),
                                         labels = c("SMB", "Outer SMB", "JR", "GLB", "34 Outside"))) %>% 
  mutate(est_scaled = est/st_index$est[1],
         est_scaled34 = est/index34$est[1],
         est_pct_lfa34 = est/lfa34)
p_index_areas <- ggplot(index_all_areas, aes(year, est_scaled, colour = region, fill = region)) + 
  geom_line(linewidth = 1.5) +
  labs(x = "Year", y = "Index value")+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous()+
  scale_colour_manual(values = c("blue","black","grey","orange","firebrick2"), guide = "none")+
  theme_light()
# Plot each subarea index as a percent of the LFA 34 index
p_pct_areas <- ggplot(index_all_areas, aes(year, est_pct_lfa34, fill = region)) + 
  geom_col(position = "fill")+
  labs(x = "Year", y = "Proportion of LFA 34 estimate", fill = "Area")+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous()+
  scale_fill_manual(values = c("blue","black","grey","orange","firebrick2"))+
  theme_light()

cowplot::plot_grid(p_index_areas, p_pct_areas,labels = c('(a)', '(b)'), label_size = 12, nrow = 1,
                   rel_widths = c(1.5, 2))