################################################
### Libraries and functions for the analyses ###
################################################

# Load libraries
library(terra) # For processing spatial data files
library(sf) # For processing spatial data files
library(concaveman) # For processing spatial data files
library(ggplot2) # For plotting
library(dplyr) # For data processing
library(tidyr) # For data processing
library(readr) # For reading in data
library(forcats) # For data processing
library(lubridate) # For data processing
library(tibble) # For data processing
library(knitr) # For making tables
library(GGally) # For plotting
library(cowplot) # For plotting
library(mgcv) # For modelling
library(sdmTMB) # For modelling
library(fmesher) # For mesh construction
library(sdmTMBextra) # For mesh construction
library(future) # For parallel processing
library(foreach) # For parallel processing
library(doFuture) # For parallel processing
library(doRNG) # For parallel processing
library(future.apply) # For parallel processing
library(Hmisc) # For %nin% operator

# Choose a set of years to plot
rep_years <- c(2005, 2010, 2013, 2017, 2020, 2023)

######## Data processing functions##########

# Reads in and processes LFA shapefile
make_shape <- function(LFA = 34:38){
  lfas <- readRDS("LFAPolysSF.rds")%>% 
    filter(LFA %in% !!LFA) %>% 
    mutate(area = as.numeric(st_area(geometry))/1000) 
  shape <- lfas %>% 
    filter(area > 200000) %>% 
    st_union() %>% 
    st_simplify(dTolerance = 500) %>% 
    st_transform(crs_UTM20N)
  islands <-  lfas %>% 
    filter(area <= 200000 & area > 2000) %>% 
    st_union() %>% 
    st_simplify(dTolerance = 500) %>% 
    st_transform(crs_UTM20N)
  shape <- st_difference(shape, islands)
  return(shape)
}

# Reads in the trawl survey data sets
read_count_data <- function(shape){
  # Lobster survey
  ILTS <- readRDS("LobsterILTS_70-82.rds")[[1]] %>% 
    as_tibble() %>% 
    rename(trip = TRIP_ID, tow = SET_NO, date = SET_DATE, year = YEAR, long = SET_LONG, lat = SET_LAT,
           area_km2 = sweptArea, gear = GEAR, count = N_Gear_Corrected, count_not_corrected = N_Not_Gear_Corrected) %>% 
    mutate(trip = as.character(trip), tow = as.character(tow)) %>% 
    filter(year <= 2023 & year >= 2005) %>% 
    select(-count_not_corrected)
  # Scallop survey
  ISAS <- read_csv("scallop_dredge_lobster_abund_data_July-17-2025.csv", col_types = "cccTifiddcdd") %>%
    select(CRUISE, TOW_NO, year, TOW_DATE, SLONG, SLAT, sample_area_m2, GEAR, abun_raw) %>%
    rename(trip = CRUISE, tow = TOW_NO, date = TOW_DATE, long = SLONG, lat = SLAT, gear = GEAR, count = abun_raw) %>%
    mutate(area_km2 = sample_area_m2/1000000, .keep = "unused", .after = lat) %>% 
    filter(year <= 2023 & year >= 2005)
  dat <- bind_rows(ILTS, ISAS) %>%
    arrange(year) %>%
    mutate(year_fac = as_factor(year), .after = year) %>%
    mutate(month = as_factor(as.numeric(month(date))), .after = year) %>% 
    mutate(gear = fct_relevel(fct(gear), "NEST", "280 BALLOON", "DREDGE"),
           survey = fct_recode(gear, lobster = "NEST", lobster = "280 BALLOON", scallop = "DREDGE"),
           moult = fct_rev(fct_collapse(month, Post = c("8", "9", "10"), Pre = c("5", "6", "7"))),
           day_of_year = yday(date),
           .before = area_km2) %>%
    mutate(log_area = log(area_km2), .after = area_km2) %>%
    mutate(count_per_km2 = count/area_km2)
  # Add UTM columns to data
  crs_UTM20N <- 32620
  dat <- add_utm_columns(dat, ll_names = c("long", "lat"), utm_crs = crs_UTM20N) %>%
    mutate(xmetres = X*1000, ymetres = Y*1000)
  dat_sf <- st_as_sf(dat, coords = c("xmetres", "ymetres"))
  st_crs(dat_sf) <- crs_UTM20N
  dat_sf <- st_intersection(dat_sf, shape)
  dat <- as.data.frame(dat_sf) %>%
    mutate(xmetres = X*1000, ymetres = Y*1000) %>%
    dplyr::select(-geometry)
  return(dat)
}

# Adds topographic and environmental covariates to data frame
add_covars <- function(dat, covars_list){
  n <- ncol(dat)
  sp_dat <- terra::vect(dat, geom = c("xmetres", "ymetres"))
  for(i in 1:length(covars_list)){
    covar <- terra::extract(covars_list[[i]], sp_dat)
    name <- names(covars_list)[i]
    dat <- dat %>% mutate(temporary = covar[,2])
    colnames(dat)[n+i] <- name
  }
  return(dat)
}

# Adds temperatures (averaged over months) to data frame
add_temps <- function(dat, temp_list){
  sp_dat <- terra::vect(dat, geom = c("xmetres", "ymetres")) 
  temperature <- NULL
  for(i in 1:length(temp_list)){
    temperature <- c(temperature, terra::extract(temp_list[[i]], sp_dat[which(sp_dat$year == 2004+i), ])[,2])
  }
  res <- cbind(dat, temperature)
  return(res)
}
# Adds temperatures for each month to data frame
add_month_temps <- function(dat, mega_list){
  sp_dat <- terra::vect(dat, geom = c("xmetres", "ymetres"))
  for(i in 1:length(mega_list)){
    rast_list <- mega_list[[i]]
    temperature <- NULL
    for(j in 1:length(rast_list)){
      temperature <- c(temperature, terra::extract(rast_list[[j]], sp_dat[which(sp_dat$year == 2004+j), ])[,2])
    }
    dat <- cbind(temperature, dat)
    names(dat)[1] <- paste0("month", i+4)
  }
  dat <- dat %>% 
    mutate(month_temps = case_when(
      month == 5 ~ month5,
      month == 6 ~ month6,
      month == 7 ~ month7,
      month == 8 ~ month8,
      month == 9 ~ month9,
      month == 10 ~ month10,
    ))
  return(dat)
}

# Makes a prediction grid for the study area (or any LFA)
make_prediction_grid <- function(covars_list, density_km = 2, sweptarea = 1,
                                 lfas = 34:38){
  grid_shape <- readRDS("bathy_by_lfa_37split_noMidas_noPass.rds") %>%
    filter(PID %in% as.character(lfas)) %>% 
    concaveman() %>% 
    st_transform(crs_UTM20N)
  grid <- grid_shape %>% 
    st_make_grid(cellsize = c(density_km, density_km), what = "centers") %>% # 4km^2 size cells default
    st_intersection(grid_shape) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    mutate(xmetres = X*1000, ymetres = Y*1000, survey = as_factor("lobster"), 
           moult = as_factor("Post"), month = 8, day_of_year = 213) # August 1
  sp_grid <- terra::vect(grid, geom = c("xmetres", "ymetres"))
  for(i in 1:length(covars_list)){
    covar <- terra::extract(covars_list[[i]], sp_grid)
    name <- names(covars_list)[i]
    grid <- grid %>% mutate(temporary = covar[,2])
    colnames(grid)[8+i] <- name
  }
  grid <- mutate(grid, "area_km2" = sweptarea)
  grid <- mutate(grid, "log_area" = log(sweptarea))
  return(grid)
}

# Makes a prediction grid for the subareas
make_prediction_grid_areas <- function(covars_list, density_m = 2000, sweptarea = 1,
                                       area){
  grid_shape <- readRDS("areas.rds") %>%
    dplyr::filter(name %in% area) %>% 
    dplyr::select(geometry) %>% 
    st_transform(crs_UTM20N)
  grid <- grid_shape %>% 
    st_make_grid(cellsize = c(density_m, density_m), what = "centers") %>% # 4km^2 size cells default
    st_intersection(grid_shape) %>% 
    st_coordinates() %>% 
    as.data.frame() %>% 
    rename(xmetres = X, ymetres = Y) %>%
    mutate(X = xmetres/1000, Y = ymetres/1000, survey = as_factor("lobster"), 
           moult = as_factor("Post"), month = 8, day_of_year = 213) # August 1
  sp_grid <- terra::vect(grid, geom = c("xmetres", "ymetres"))
  for(i in 1:length(covars_list)){
    covar <- terra::extract(covars_list[[i]], sp_grid)
    name <- names(covars_list)[i]
    grid <- grid %>% mutate(temporary = covar[,2])
    colnames(grid)[8+i] <- name
  }
  grid <- mutate(grid, "area_km2" = sweptarea)
  grid <- mutate(grid, "log_area" = log(sweptarea))
  return(grid)
}

# Makes a barrier mesh 
make_barrier_mesh <- function(dat, map){
  mesh_dat <- as.matrix(dat[dat$year == 2023, c("X", "Y")])
  bnd <- fm_nonconvex_hull(mesh_dat, convex = -0.075)
  nonconvex_mesh <- fm_mesh_2d(
    boundary = bnd,
    cutoff = 4,
    max.edge = c(20, 60),
    offset = c(20, 60)
  )
  mesh <- make_mesh(dat, c("X", "Y"), mesh = nonconvex_mesh)
  plot(mesh)
  barrier_mesh <- add_barrier_mesh(mesh, map, plot = TRUE, proj_scaling = 1000)
  return(barrier_mesh)
}

# Calculates variance inflation factor
vif <- function(x) {
  v <- vapply(seq_along(x), function(i) {
    rsq <- summary(lm(x[[i]] ~ . , data = x[,-i, drop = FALSE]))$r.squared
    1 / (1 - rsq)
  }, FUN.VALUE = numeric(1L))
  setNames(v, colnames(x))
}

########### Plotting functions #########

# Plots model residuals as a histogram, QQ-plot, or map
plot_resids_spatial <- function(data, plot_type = c("histogram", "qq-plot", "map")){
  if(plot_type == "histogram"){
    plot <- ggplot(data, aes(spresids))+
      geom_histogram()+
      labs(x = "Randomized quantile residuals")
  } else if(plot_type == "qq-plot"){
    plot <- ggplot(data, aes(sample = spresids))+
      geom_qq()+
      geom_qq_line()+
      labs(x = "Theoretical", y = "Sample")
  } else {
    plot <- ggplot(map_fundy) + geom_sf() + 
      geom_point(data = data, mapping = aes(x = xmetres, y = ymetres, col = spresids), alpha = 0.5, size = 0.5) + 
      scale_colour_gradient2() +
      labs(x = "Longitude", y = "Latitude", col = "Residual")
  }
  plot <- plot+theme_light()
  return(plot)
}

# Plots the abundance index
plot_index <- function(dat){
  p <- ggplot(dat, aes(year, est)) + geom_line() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
    labs(x = "Year", y = "Index value")+
    scale_x_continuous(breaks = seq(min(dat$year), max(dat$year), by = 1))+
    scale_y_continuous()+
    theme_light()
  return(p)
}

# Maps the predicted abundance, standard deviation, coefficient of variation, or random field values
plot_map <- function(dat, years = "all", 
                     return = c("link", "response", "responselog", "spatial", "spatiotemporal", 
                                "spatial-varying", "sd", "cv"),
                     axes = TRUE) {
  if(return == "link"){
    dat$response <- dat$est
    label = "Log counts"
    trans = "identity"
  } else if (return == "spatial") {
    dat$response <- dat$omega_s
    label = "Spatial RF"
    trans = "identity"
  } else if (return == "spatiotemporal"){
    dat$response <- dat$epsilon_st
    label = "Spatiotemporal RF"
    trans = "identity"
  } else if (return == "spatial-varying"){
    dat$response <- dat$zeta_s_surveyscallop
    label = "Spatial-varying survey RF"
    trans = "identity"
  } else if (return == "sd"){
    dat$response <- dat$sd
    label = "Standard Deviation"
    trans = "identity"
  } else if (return == "cv"){
    dat$response <- dat$cv
    label = "Coefficient of Variation"
    trans = "identity"
  } else if (return == "response") {
    dat$response <- exp(dat$est)
    label = "Counts/km2"
    trans = "identity"
  } else {
    dat$response <- exp(dat$est)
    label = "Response (counts/km2)"
    trans = "log10"
  }
  if(is.numeric(years)){
    dat <- filter(dat, year %in% years)
  }
  p <- ggplot(map_fundy) +
    geom_sf() +
    geom_tile(data = dat, aes(xmetres, ymetres, fill = response)) +
    coord_sf()+
    scale_fill_viridis_c(option = "B", direction = -1, transform = trans)+
    labs(x = "", y = "", fill = label)+
    theme_light()+
    coord_sf(xlim = c(130218, 461283),  ylim = c(4705414, 5094169), expand = F)+
    scale_x_continuous(breaks = c(-67, -65))+
    scale_y_continuous(breaks = c(43, 44, 45))
  if(length(years) != 1){
    p <- p+facet_wrap(~year)
  }
  if(axes == FALSE){
    p <- p+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())
  } else {
    p <- p+theme(axis.text = element_text(size = rel(0.5)))
  }
  return(p)
}

# Plots the simulated indices against the operating model index
plot_sim_indices <- function(true, list){
  index_list <- lapply(list, function(x) x[[6]])
  est_list <- lapply(index_list, function(x) dplyr::select(x, est))
  dat <- data.frame(year = 2005:2023, true = true)
  for(i in seq(length(est_list))){
    est <- est_list[[i]]
    names(est) <- paste("est", as.character(i))
    dat <- dplyr::mutate(dat, est)
  }
  gc()
  dat <- pivot_longer(dat, !year, names_to = "series", values_to = "index")
  p <- ggplot()+
    geom_line(data = dat[dat[ , 2] != "true", ], mapping = aes(year, index, col = series, alpha = 0.2)) +
    geom_line(data = dat[dat[ , 2] == "true", ], mapping = aes(year, index), colour = "black", linewidth = 2)+
    scale_x_continuous(breaks = seq(min(dat$year), max(dat$year), by = 1))+
    scale_y_continuous()+
    theme_light()+
    theme(legend.position="none")+
    labs(x = "Year", y = "Index value")
  return(p)
}

########### Simulation functions ##########
# Simulates negative binomial-distributed count data at all sampling locations
get_nb_sims <- function(mod, template, type = c("mle-eb", "mle-mvn")){
  sim_dat <- simulate(mod, type = type, newdata = template,
                      offset = template$log_area)
  sim_dat <- cbind(template, sim = as.numeric(sim_dat))
  return(sim_dat)
}
# Estimates parameters and abundance index from simulated data
get_pars_and_index <- function(dat, mesh, grid, area = 4, scale, model = "combined", fix, bc = FALSE){
  if(model == "scallop"){
    sim_mod <- sdmTMB(sim ~ 0+year_fac+moult+s(logdepth_std)+s(broad_BPI_std), 
                      family = nbinom2(), data = dat, mesh = mesh, 
                      offset = dat$log_area+fix, extra_time = 2020,
                      spatial = "on", spatiotemporal = "ar1",
                      share_range = TRUE, time = "year", do_index = TRUE, 
                      predict_args = list(newdata = grid, offset = grid$log_area), 
                      index_args = list(area = area))
  } else if (model == "lobster"){
    sim_mod <- sdmTMB(sim ~ 0+year_fac+moult+s(logdepth_std)+broad_BPI_std, 
                      family = nbinom2(), data = dat, mesh = mesh,
                      offset = dat$log_area,
                      spatial = "on", spatiotemporal = "ar1",
                      share_range = TRUE, time = "year", do_index = TRUE, 
                      predict_args = list(newdata = grid, offset = grid$log_area), 
                      index_args = list(area = area))
  } else {
    sim_mod <- sdmTMB(sim ~ 0+year_fac+survey+moult+s(logdepth_std)+s(broad_BPI_std), 
                      family = nbinom2(), data = dat, mesh = mesh,
                      offset = dat$log_area, 
                      spatial_varying = ~survey,
                      spatial = "on", spatiotemporal = "ar1",
                      share_range = TRUE, time = "year", do_index = TRUE, 
                      predict_args = list(newdata = grid, offset = grid$log_area), 
                      index_args = list(area = area))
  }
  gc()
  ln_phi <- sim_mod$parlist$ln_phi
  ln_tau_O <- sim_mod$parlist$ln_tau_O
  ln_tau_E <- sim_mod$parlist$ln_tau_E
  ln_kappa <- sim_mod$parlist$ln_kappa
  if(model == "combined"){
    betascallop <- sim_mod$parlist$b_j[20]
  } else {
    betascallop <- NA
  }
  sim_index <- get_index(sim_mod, bias_correct = bc)
  sim_index_scaled <- sim_index %>%
    select(year, est, lwr, upr, se_natural) %>% 
    mutate(est = est/scale,
           lwr = lwr/scale,
           upr = upr/scale,
           se_natural = se_natural/scale)
  rm(sim_mod, sim_index)
  if(model == "scallop"){
    sim_index_scaled <- rbind(sim_index_scaled[1:15, ], c(2020, NA, NA, NA, NA), sim_index_scaled[16:18,])
    dplyr::mutate(sim_index_scaled, est = as.numeric(est),
                  lwr = as.numeric(lwr),
                  upr = as.numeric(upr),
                  se_natural = as.numeric(se_natural))
    rownames(sim_index_scaled) <- 1:nrow(sim_index_scaled)
  }
  gc()
  list <- list(ln_phi, ln_tau_O, ln_tau_E, ln_kappa, betascallop,
               sim_index_scaled)
  return(list)
}

########### Performance metric functions ##########
# Calculates the relative errors of the simulated indices vs. the true index
get_rel_errors <- function(true, list, n = 100){
  index_list <- lapply(list, function(x) x[[6]])
  est_list <- lapply(index_list, function(x) dplyr::select(x, est))
  rel_errors <- lapply(est_list, function(x) (x-true)/true)
  rel_errors <- matrix(unlist(rel_errors), ncol = n, byrow = FALSE)
}

# Calculates the RMSEs
get_RMSEs <- function(true, list, n = 100){
  index_list <- lapply(list, function(x) x[[6]])
  est_list <- lapply(index_list, function(x) dplyr::select(x, est))
  sq_errors <- lapply(est_list, function(x) (x-true)^2)
  sq_errors <- matrix(unlist(sq_errors), ncol = n, byrow = FALSE)
  RMSEs <- rowMeans(sq_errors)
}

# Calculates the 95% confidence interval coverage rates
get_CovRates <- function(true, list, n = 100){
  index_list <- lapply(list, function(x) x[[6]])
  lower_bounds <- lapply(index_list, function(x) dplyr::select(x, lwr))
  upper_bounds <- lapply(index_list, function(x) dplyr::select(x, upr))
  contains_true <- mapply(function(x, y) x <= true & y >= true, lower_bounds, upper_bounds, SIMPLIFY = FALSE)
  contains_true <- matrix(unlist(contains_true), ncol = n, byrow = FALSE)
  coverage_rates <- 100*rowMeans(contains_true)
}

# Calculates the 95% confidence interval average widths
get_AveWidths <- function(list, n = 100){
  index_list <- lapply(list, function(x) x[[6]])
  lower_bounds <- lapply(index_list, function(x) dplyr::select(x, lwr))
  upper_bounds <- lapply(index_list, function(x) dplyr::select(x, upr))
  widths <- mapply(function(x, y) y-x, lower_bounds, upper_bounds, SIMPLIFY = FALSE)
  widths <- matrix(unlist(widths), ncol = n, byrow = FALSE)
  ave_widths <- rowMeans(widths)
}