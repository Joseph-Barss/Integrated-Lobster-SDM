#######################
### Data processing ###
#######################

# Read in coastline
map_data <- read_sf("AtlanicUSCdnCoast_v2/AtlanicUSCdnCoast_v2.shp")$geometry
map_data <- st_make_valid(map_data)
# Focus on Bay of Fundy area
box <- c(xmin = -68, ymin = 42, xmax = -63, ymax = 46.5)
map_fundy <- st_crop(map_data, box)
# Transform the CRS (needed for sdmTMB)
crs_WG84 <- st_crs(map_data)
crs_UTM20N <- 32620
map_data <- st_transform(map_data, crs_UTM20N)
map_fundy <- st_transform(map_fundy, crs_UTM20N)
# Make polygon of the study area
shape <- make_shape()

# Read in covariate rasters
bathymetry <- rast("OFI_lobster_covariates/ofi_DEM_UTMZ20.tif")
broad_BPI <- rast("OFI_lobster_covariates/sbbpi_25x50.tif")
fine_BPI <- rast("OFI_lobster_covariates/sfbpi_1x20.tif")
DMV <- rast("OFI_lobster_covariates/sdmv_20circlenw.tif")
slope <- rast("OFI_lobster_covariates/slope_11x11nw.tif")
eastness <- rast("OFI_lobster_covariates/eastness_11x11nw.tif")
northness <- rast("OFI_lobster_covariates/northness_11x11nw.tif")
curvature_profile <- rast("OFI_lobster_covariates/profc_11x11nw.tif")
curvature_planform <- rast("OFI_lobster_covariates/planc_11x11nw.tif")
curvature_twisting <- rast("OFI_lobster_covariates/twistc_11x11nw.tif")
curvature_mean <- rast("OFI_lobster_covariates/meanc_11x11nw.tif")
curvature_max <- rast("OFI_lobster_covariates/maxc_11x11nw.tif")
curvature_min <- rast("OFI_lobster_covariates/minc_11x11nw.tif")
VRM <- rast("OFI_lobster_covariates/vrm_11x11nw.tif")
SAPA <- rast("OFI_lobster_covariates/sapa_11x11nw.tif")
adj_SD <- rast("OFI_lobster_covariates/adj_SD_11x11nw.tif")
RIE <- rast("OFI_lobster_covariates/rie_11x11nw.tif")
sed_mob_freq <- rast("Li_et_al_2024_covariates/smf_utmz20.tif")
comb_shear_vel <- rast("Li_et_al_2024_covariates/ust_utmz20.tif")
total_curr_spd <- rast("Li_et_al_2024_covariates/Utotl_utmz20.tif")
grainsize <- rast("Li_et_al_2024_covariates/gs_phi_utmz20.tif")
seabed_dist <- rast("Li_et_al_2024_covariates/SDI_utmz20.tif")
sed_mob_ind <- rast("Li_et_al_2024_covariates/SMI_utmz20.tif")
# Put covariates in list
covars_list <- list(bathymetry = bathymetry, broad_BPI = broad_BPI, 
                    fine_BPI = fine_BPI, DMV = DMV, slope = slope,
                    eastness = eastness, northness = northness, 
                    curvature_profile = curvature_profile, 
                    curvature_planform = curvature_planform,
                    curvature_twisting = curvature_twisting, 
                    curvature_mean = curvature_mean, 
                    curvature_max = curvature_max,
                    curvature_min = curvature_min, VRM = VRM, SAPA = SAPA,
                    adj_SD = adj_SD, RIE = RIE, sed_mob_freq = sed_mob_freq,
                    comb_shear_vel = comb_shear_vel,
                    total_curr_spd = total_curr_spd,
                    grainsize = grainsize,
                    seabed_dist = seabed_dist,
                    sed_mob_ind = sed_mob_ind)
# Read in temperatures
avg_temp_list <- list()
may_list <- list()
jun_list <- list()
jul_list <- list()
aug_list <- list()
sep_list <- list()
oct_list <- list()
path <- "OFI_lobster_covariates/Monthly_BtmTemp_cropped_to_BOF_and_GOM_domain/UTMz20/"
for(i in 1:19){
  yr <- 2004+i
  temp_may <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_5.tif"))
  temp_jun <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_6.tif"))
  temp_jul <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_7.tif"))
  temp_aug <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_8.tif"))
  temp_sep <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_9.tif"))
  temp_oct <- rast(paste0(path, yr, "/BOF_GOM_BottomT_", yr, "_10.tif"))
  avg_temp_list[[i]] <- terra::mean(temp_may, temp_jun, temp_jul, temp_aug,
                                    temp_sep, temp_oct, na.rm = TRUE)
  may_list[[i]] <- temp_may 
  jun_list[[i]] <- temp_jun 
  jul_list[[i]] <- temp_jul 
  aug_list[[i]] <- temp_aug 
  sep_list[[i]] <- temp_sep 
  oct_list[[i]] <- temp_oct 
}
# Put monthly temperatures in list
list_of_lists <- list(may_list, jun_list, jul_list, aug_list,
                      sep_list, oct_list)

# Read in and process the trawl survey data
lob_data <- read_count_data(shape)
lob_dat_covars <- add_covars(lob_data, covars_list)
lob_dat_covars <- add_temps(lob_dat_covars, avg_temp_list)
lob_dat_covars <- add_month_temps(lob_dat_covars, list_of_lists)

lob_mod_dat <- lob_dat_covars %>% 
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>%
  filter(!if_any(bathymetry:temperature, is.na) & bathymetry < 3.5) %>% 
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

################# Make prediction grid ###########
#Make prediction grids, removing depths under 3.5m
grid <- make_prediction_grid(covars_list,
                             density_km = 2, sweptarea = 1)
grid_sf <- terra::vect(grid, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], grid_sf)[,2]
grid <- cbind(grid, temp23)
grid <- filter(grid, bathymetry < -3.5 & !is.na(bathymetry))
grid_yrs_all_covars <- replicate_df(grid, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
grid_yrs <- grid_yrs_all_covars %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(X:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)

################# Make grids for each LFA (37 split between 36 and 38) #########
# LFA 34
grid34 <- make_prediction_grid(covars_list,
                               density_km = 2, sweptarea = 1, lfas = 34)
grid34_sf <- terra::vect(grid34, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], grid34_sf)[,2]
grid34 <- cbind(grid34, temp23)
grid34 <- filter(grid34, bathymetry < -3.5 & !is.na(bathymetry))
grid34_yrs <- replicate_df(grid34, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists) 
grid34_yrs <- grid34_yrs %>% 
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(X:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# LFA 35
grid35 <- make_prediction_grid(covars_list,
                               density_km = 2, sweptarea = 1, lfas = 35)
grid35_sf <- terra::vect(grid35, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], grid35_sf)[,2]
grid35 <- cbind(grid35, temp23)
grid35 <- filter(grid35, bathymetry < -3.5 & !is.na(bathymetry))
grid35_yrs <- replicate_df(grid35, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists) 
grid35_yrs <- grid35_yrs %>% 
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(X:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# LFA 36
grid36 <- make_prediction_grid(covars_list,
                               density_km = 2, sweptarea = 1, lfas = 36)
grid36_sf <- terra::vect(grid36, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], grid36_sf)[,2]
grid36 <- cbind(grid36, temp23)
grid36 <- filter(grid36, bathymetry < -3.5 & !is.na(bathymetry))
grid36_yrs <- replicate_df(grid36, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists) 
grid36_yrs <- grid36_yrs %>% 
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(X:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# LFA 38
grid38 <- make_prediction_grid(covars_list,
                               density_km = 2, sweptarea = 1, lfas = 38)
grid38_sf <- terra::vect(grid38, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], grid38_sf)[,2]
grid38 <- cbind(grid38, temp23)
grid38 <- filter(grid38, bathymetry < -3.5 & !is.na(bathymetry))
grid38_yrs <- replicate_df(grid38, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists) 
grid38_yrs <- grid38_yrs %>%  
  dplyr::select(-c(DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
                   total_curr_spd, grainsize, sed_mob_ind,
                   month5, month6, month7, month8, month9, month10)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(X:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)

########### Make grids of subareas for detailed analysis #########
statgrids <- readRDS("GridPolys_DepthPruned_37Split.rds") %>% 
  filter(LFA == 34)
SMB <- st_union(st_buffer(filter(statgrids, GRID_NO %in% c("69", "81", "92")), 0))
OSMB <- st_union(st_buffer(filter(statgrids, GRID_NO %in% c("103", "114", "125", "126")), 0))
GLB <- st_union(st_buffer(filter(statgrids, GRID_NO %in% c("127", "140", "141", "156", "157")), 0))
JR <- st_union(st_buffer(filter(statgrids, GRID_NO %in% c("138", "139", "154", "155")), 0))
OUT <- st_union(st_buffer(filter(statgrids, GRID_NO %nin% c("69", "81", "92", "103", "114", "125", "126",
                                                            "127", "140", "141", "156", "157", "138", "139", "154", "155")), 0))
areas <- vect(list(vect(SMB), vect(OSMB), vect(GLB), vect(JR), vect(OUT)))
areas_sf <- st_as_sf(areas) %>% 
  mutate(name = c("SMB", "OSMB", "GLB", "JR", "OUT")) %>% 
  st_transform(crs_UTM20N)

#Make prediction grids, removing depths under 3.5m
# Saint Mary's Bay
gridSMB <- make_prediction_grid_areas(covars_list,
                                      density_m = 2000, sweptarea = 1,
                                      area = "SMB")
gridSMB_sf <- terra::vect(gridSMB, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], gridSMB_sf)[,2]
gridSMB <- cbind(gridSMB, temp23)
gridSMB <- filter(gridSMB, bathymetry < -3.5 & !is.na(bathymetry))
gridSMB_yrs <- replicate_df(gridSMB, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
gridSMB_yrs <- gridSMB_yrs %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(xmetres:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# Outer Saint Mary's Bay
gridOSMB <- make_prediction_grid_areas(covars_list,
                                       density_m = 2000, sweptarea = 1,
                                       area = "OSMB")
gridOSMB_sf <- terra::vect(gridOSMB, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], gridOSMB_sf)[,2]
gridOSMB <- cbind(gridOSMB, temp23)
gridOSMB <- filter(gridOSMB, bathymetry < -3.5 & !is.na(bathymetry))
gridOSMB_yrs <- replicate_df(gridOSMB, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
gridOSMB_yrs <- gridOSMB_yrs %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(xmetres:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# Jacquard's Ridge
gridJR <- make_prediction_grid_areas(covars_list,
                                     density_m = 2000, sweptarea = 1,
                                     area = "JR")
gridJR_sf <- terra::vect(gridJR, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], gridJR_sf)[,2]
gridJR <- cbind(gridJR, temp23)
gridJR <- filter(gridJR, bathymetry < -3.5 & !is.na(bathymetry))
gridJR_yrs <- replicate_df(gridJR, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
gridJR_yrs <- gridJR_yrs %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(xmetres:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# Greater Lobster Bay
gridGLB <- make_prediction_grid_areas(covars_list,
                                      density_m = 2000, sweptarea = 1,
                                      area = "GLB")
gridGLB_sf <- terra::vect(gridGLB, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], gridGLB_sf)[,2]
gridGLB <- cbind(gridGLB, temp23)
gridGLB <- filter(gridGLB, bathymetry < -3.5 & !is.na(bathymetry))
gridGLB_yrs <- replicate_df(gridGLB, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
gridGLB_yrs <- gridGLB_yrs %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(xmetres:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
# LFA 34 outside of the named areas
gridOUT <- make_prediction_grid_areas(covars_list,
                                      density_m = 2000, sweptarea = 1,
                                      area = "OUT")
gridOUT_sf <- terra::vect(gridOUT, geom = c("xmetres", "ymetres")) 
temp23 <- terra::extract(avg_temp_list[[19]], gridOUT_sf)[,2]
gridOUT <- cbind(gridOUT, temp23)
gridOUT <- filter(gridOUT, bathymetry < -3.5 & !is.na(bathymetry))
gridOUT_yrs <- replicate_df(gridOUT, "year", unique(lob_mod_dat$year)) %>% 
  mutate(year_fac = as_factor(year)) %>% 
  add_month_temps(list_of_lists)
gridOUT_yrs <- gridOUT_yrs %>% 
  select(-c(temp23, DMV, curvature_mean, curvature_max, curvature_min, VRM, SAPA, RIE, sed_mob_freq, comb_shear_vel,
            total_curr_spd, grainsize, sed_mob_ind)) %>% 
  mutate(logdepth = log(-bathymetry), sqrtslope = sqrt(slope), logadjSD = log(adj_SD)) %>% 
  mutate(logdepth_std = (logdepth-mean(lob_mod_dat$logdepth))/sd(lob_mod_dat$logdepth),
         broad_BPI_std = (broad_BPI-mean(lob_mod_dat$broad_BPI))/sd(lob_mod_dat$broad_BPI),
         fine_BPI_std = (fine_BPI-mean(lob_mod_dat$fine_BPI))/sd(lob_mod_dat$fine_BPI),
         sqrtslope_std = (sqrtslope-mean(lob_mod_dat$sqrtslope))/sd(lob_mod_dat$sqrtslope),
         curv_prof_std = curvature_profile/sd(lob_mod_dat$curvature_profile),
         curv_plan_std = curvature_planform/sd(lob_mod_dat$curvature_planform),
         curv_twist_std = curvature_twisting/sd(lob_mod_dat$curvature_twisting),
         logadjSD_std = (logadjSD-mean(lob_mod_dat$logadjSD))/sd(lob_mod_dat$logadjSD),
         temp_std = (month_temps-mean(lob_mod_dat$month_temps))/sd(lob_mod_dat$month_temps)) %>% 
  select(xmetres:day_of_year, area_km2:year_fac, logdepth_std:temp_std, seabed_dist)
