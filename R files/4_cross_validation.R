########################
### Cross-validation ###
########################

# Make the mesh
mesh <- make_mesh(lob_mod_dat, xy_cols = c("X", "Y"), cutoff = 5,
                  max.edge = c(40, 80))
plot(mesh)
barrier_mesh <- add_barrier_mesh(mesh, st_as_sf(map_fundy), proj_scaling = 1000, plot = TRUE)

# Cross-validation with 5-fold CV
# Run cross-validation in parallel
# Start by examining models with many covariates, and different random fields
plan(multisession, workers = 5)
set.seed(11)
# No random fields
norf_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                        s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                        logadjSD_std+s(seabed_dist)+temp_std,
                      family = nbinom2(), 
                      data = lob_mod_dat, mesh = barrier_mesh,
                      offset = "log_area",
                      spatial = "off", spatiotemporal = "off", 
                      parallel = TRUE, k_folds = 5)
set.seed(11)
# Only spatial random field
spat_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                        s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                        logadjSD_std+s(seabed_dist)+temp_std,
                      family = nbinom2(), 
                      data = lob_mod_dat, mesh = barrier_mesh,
                      offset = "log_area",
                      spatial = "on", spatiotemporal = "off", 
                      parallel = TRUE, k_folds = 5)
set.seed(11)
# Spatial and IID spatiotemporal random fields
iid_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                       s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                       logadjSD_std+s(seabed_dist)+temp_std,
                     family = nbinom2(), time = "year",
                     data = lob_mod_dat, mesh = barrier_mesh,
                     offset = "log_area",
                     spatial = "on", spatiotemporal = "iid",
                     share_range = FALSE, parallel = TRUE, 
                     k_folds = 5)
set.seed(11)
# Spatial and AR(1) spatiotemporal random fields
ar1_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                       s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                       logadjSD_std+s(seabed_dist)+temp_std,
                     family = nbinom2(), time = "year",
                     data = lob_mod_dat, mesh = barrier_mesh,
                     offset = "log_area",
                     spatial = "on", spatiotemporal = "ar1",
                     share_range = TRUE, parallel = TRUE, 
                     k_folds = 5)
set.seed(11)
# Spatial RF and spatially varying coefficient
svc_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                       s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                       logadjSD_std+s(seabed_dist)+temp_std,
                     family = nbinom2(), time = "year",
                     spatial_varying = ~survey,
                     data = lob_mod_dat, mesh = barrier_mesh,
                     offset = "log_area",
                     spatial = "on", spatiotemporal = "off",
                     share_range = TRUE, parallel = TRUE, 
                     k_folds = 5)
set.seed(11)
# Spatial RF, IID spatiotemporal RF, and spatially varying coefficient
svciid_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                          s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                          logadjSD_std+s(seabed_dist)+temp_std,
                        family = nbinom2(), time = "year",
                        spatial_varying = ~survey,
                        data = lob_mod_dat, mesh = barrier_mesh,
                        offset = "log_area",
                        spatial = "on", spatiotemporal = "iid",
                        share_range = TRUE, parallel = TRUE, 
                        k_folds = 5)
set.seed(11)
# Spatial RF, AR(1) spatiotemporal RF, and spatially varying coefficient
svcar1_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                          s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                          logadjSD_std+s(seabed_dist)+temp_std,
                        family = nbinom2(), time = "year",
                        spatial_varying = ~survey,
                        data = lob_mod_dat, mesh = barrier_mesh,
                        offset = "log_area",
                        spatial = "on", spatiotemporal = "ar1",
                        share_range = TRUE, parallel = TRUE, 
                        k_folds = 5)
plan(sequential)
# Check sum of log-likelihoods
norf_kcv$sum_loglik # -16222.73
spat_kcv$sum_loglik # -14515.7
iid_kcv$sum_loglik # -14366.85
ar1_kcv$sum_loglik # -14247.25
svc_kcv$sum_loglik # -14472.95
svciid_kcv$sum_loglik # -14318.27
svcar1_kcv$sum_loglik # -14215.4

# Now we choose the covariates by leaving out one at a time
plan(multisession, workers = 5)
set.seed(11)
# no temperature
notemp_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                          s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                          logadjSD_std+s(seabed_dist),
                        family = nbinom2(), time = "year",
                        spatial_varying = ~survey,
                        data = lob_mod_dat, mesh = barrier_mesh,
                        offset = "log_area",
                        spatial = "on", spatiotemporal = "ar1",
                        share_range = TRUE, parallel = TRUE, 
                        k_folds = 5)
set.seed(11)
# no curvature
nocurv_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                          s(broad_BPI_std)+sqrtslope_std+
                          logadjSD_std+s(seabed_dist)+temp_std,
                        family = nbinom2(), time = "year",
                        spatial_varying = ~survey,
                        data = lob_mod_dat, mesh = barrier_mesh,
                        offset = "log_area",
                        spatial = "on", spatiotemporal = "ar1",
                        share_range = TRUE, parallel = TRUE, 
                        k_folds = 5)
set.seed(11)
# no slope
noslope_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                           s(broad_BPI_std)+curv_prof_std+
                           logadjSD_std+s(seabed_dist)+temp_std,
                         family = nbinom2(), time = "year",
                         spatial_varying = ~survey,
                         data = lob_mod_dat, mesh = barrier_mesh,
                         offset = "log_area",
                         spatial = "on", spatiotemporal = "ar1",
                         share_range = TRUE, parallel = TRUE, 
                         k_folds = 5)
set.seed(11)
# no adjusted SD
nologadjsd_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                              s(broad_BPI_std)+sqrtslope_std+curv_prof_std
                            +s(seabed_dist)+temp_std,
                            family = nbinom2(), time = "year",
                            spatial_varying = ~survey,
                            data = lob_mod_dat, mesh = barrier_mesh,
                            offset = "log_area",
                            spatial = "on", spatiotemporal = "ar1",
                            share_range = TRUE, parallel = TRUE, 
                            k_folds = 5)
set.seed(11)
# No seabed disturbance
noseabeddist_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                                logadjSD_std+temp_std,
                              family = nbinom2(), time = "year",
                              spatial_varying = ~survey,
                              data = lob_mod_dat, mesh = barrier_mesh,
                              offset = "log_area",
                              spatial = "on", spatiotemporal = "ar1",
                              share_range = TRUE, parallel = TRUE, 
                              k_folds = 5)
set.seed(11)
# No broad BPI
nobpi_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)
                       +sqrtslope_std+curv_prof_std+
                         logadjSD_std+s(seabed_dist)+temp_std,
                       family = nbinom2(), time = "year",
                       spatial_varying = ~survey,
                       data = lob_mod_dat, mesh = barrier_mesh,
                       offset = "log_area",
                       spatial = "on", spatiotemporal = "ar1",
                       share_range = TRUE, parallel = TRUE, 
                       k_folds = 5)
set.seed(11)
# No moult variable
nomoult_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+s(logdepth_std)+
                           s(broad_BPI_std)+sqrtslope_std+curv_prof_std+
                           logadjSD_std+s(seabed_dist)+temp_std,
                         family = nbinom2(), time = "year",
                         spatial_varying = ~survey,
                         data = lob_mod_dat, mesh = barrier_mesh,
                         offset = "log_area",
                         spatial = "on", spatiotemporal = "ar1",
                         share_range = TRUE, parallel = TRUE, 
                         k_folds = 5)
plan(sequential)
# Check sum of log-likelihoods
notemp_kcv$sum_loglik # -14214.73
nocurv_kcv$sum_loglik # -14216.32
noslope_kcv$sum_loglik # -14215.72
nologadjsd_kcv$sum_loglik # -14215.15
noseabeddist_kcv$sum_loglik # -14214.14
nobpi_kcv$sum_loglik # -14227.12
nomoult_kcv$sum_loglik # -14217.16

# Now try different combinations of covariates
plan(multisession, workers = 5)
set.seed(11)
moultdepthbpislopecurv_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                          s(broad_BPI_std)+sqrtslope_std+curv_prof_std,
                                        family = nbinom2(), time = "year",
                                        spatial_varying = ~survey,
                                        data = lob_mod_dat, mesh = barrier_mesh,
                                        offset = "log_area",
                                        spatial = "on", spatiotemporal = "ar1",
                                        share_range = TRUE, parallel = TRUE, 
                                        k_folds = 5)
set.seed(11)
moultdepthbpislope_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                      s(broad_BPI_std)+sqrtslope_std,
                                    family = nbinom2(), time = "year",
                                    spatial_varying = ~survey,
                                    data = lob_mod_dat, mesh = barrier_mesh,
                                    offset = "log_area",
                                    spatial = "on", spatiotemporal = "ar1",
                                    share_range = TRUE, parallel = TRUE, 
                                    k_folds = 5)
set.seed(11)
moultdepthbpicurv_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                     s(broad_BPI_std)+curv_prof_std,
                                   family = nbinom2(), time = "year",
                                   spatial_varying = ~survey,
                                   data = lob_mod_dat, mesh = barrier_mesh,
                                   offset = "log_area",
                                   spatial = "on", spatiotemporal = "ar1",
                                   share_range = TRUE, parallel = TRUE, 
                                   k_folds = 5)
set.seed(11)
moultdepthslopecurv_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                       sqrtslope_std+curv_prof_std,
                                     family = nbinom2(), time = "year",
                                     spatial_varying = ~survey,
                                     data = lob_mod_dat, mesh = barrier_mesh,
                                     offset = "log_area",
                                     spatial = "on", spatiotemporal = "ar1",
                                     share_range = TRUE, parallel = TRUE, 
                                     k_folds = 5)
set.seed(11)
depthbpislopecurv_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+s(logdepth_std)+
                                     s(broad_BPI_std)+sqrtslope_std+curv_prof_std,
                                   family = nbinom2(), time = "year",
                                   spatial_varying = ~survey,
                                   data = lob_mod_dat, mesh = barrier_mesh,
                                   offset = "log_area",
                                   spatial = "on", spatiotemporal = "ar1",
                                   share_range = TRUE, parallel = TRUE, 
                                   k_folds = 5)
set.seed(11)
moultdepthbpi_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std)+
                                 s(broad_BPI_std),
                               family = nbinom2(), time = "year",
                               spatial_varying = ~survey,
                               data = lob_mod_dat, mesh = barrier_mesh,
                               offset = "log_area",
                               spatial = "on", spatiotemporal = "ar1",
                               share_range = TRUE, parallel = TRUE, 
                               k_folds = 5)
set.seed(11)
moultdepth_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+moult+s(logdepth_std),
                            family = nbinom2(), time = "year",
                            spatial_varying = ~survey,
                            data = lob_mod_dat, mesh = barrier_mesh,
                            offset = "log_area",
                            spatial = "on", spatiotemporal = "ar1",
                            share_range = TRUE, parallel = TRUE, 
                            k_folds = 5)
set.seed(11)
depthbpi_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+s(logdepth_std)+
                            s(broad_BPI_std),
                          family = nbinom2(), time = "year",
                          spatial_varying = ~survey,
                          data = lob_mod_dat, mesh = barrier_mesh,
                          offset = "log_area",
                          spatial = "on", spatiotemporal = "ar1",
                          share_range = TRUE, parallel = TRUE, 
                          k_folds = 5)
set.seed(11)
depth_kcv <- sdmTMB_cv(count ~ 0+year_fac+survey+s(logdepth_std),
                       family = nbinom2(), time = "year",
                       spatial_varying = ~survey,
                       data = lob_mod_dat, mesh = barrier_mesh,
                       offset = "log_area",
                       spatial = "on", spatiotemporal = "ar1",
                       share_range = TRUE, parallel = TRUE, 
                       k_folds = 5)
plan(sequential)
# Check sum of log-likelihoods
moultdepthbpislopecurv_kcv$sum_loglik # -14213.01
moultdepthbpislope_kcv$sum_loglik # -14214.16
moultdepthbpicurv_kcv$sum_loglik # -14212.97 # best
moultdepthslopecurv_kcv$sum_loglik # -14226.45
depthbpislopecurv_kcv$sum_loglik # -14220.25
moultdepthbpi_kcv$sum_loglik # -14213.99 # simplest of good models
moultdepth_kcv$sum_loglik # -14226.19
depthbpi_kcv$sum_loglik # -14220.83
depth_kcv$sum_loglik # -14233.58

# Keep survey, moult, logdepth, and broad BPI
