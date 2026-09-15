###################
### Simulations ###
###################

# Note: The simulation study code should ideally be run
# using a high-performance computer with multiple cores.

# Make prediction grid
sim_grid_yrs <- select(grid_yrs, X, Y,
                       year, year_fac, survey, log_area, day_of_year,
                       logdepth_std, broad_BPI_std, moult,
                       seabed_dist, temp_std)
# Make data frames with covariates at real survey locations
all_template <- select(lob_mod_dat, X, Y,
                       year, year_fac, survey, log_area, day_of_year,
                       logdepth_std, broad_BPI_std, moult,
                       seabed_dist, temp_std)
scallop_template <- filter(all_template, survey == "scallop")
lobster_template <- filter(all_template, survey == "lobster")

# Simulate negative binomial distributed data sets 100 times
nsims <- 100
registerDoFuture()
registerDoRNG(10)
future::plan(multisession, workers = 10)
sim_nb_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(stmod, all_template, type = "mle-mvn")

scallop_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(stmod, scallop_template, type = "mle-mvn")

lobster_list <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(stmod, lobster_template, type = "mle-mvn")

scallop_list_single <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(scallopmod, scallop_template, type = "mle-mvn")

lobster_list_single <- foreach(i = 1:nsims, .options.future = list(seed = TRUE)) %dopar%
  get_nb_sims(lobstermod, lobster_template, type = "mle-mvn")

future::plan(sequential)

# Indices need scaling, and scallop index needs catchability adjustment
scale = st_index$est[1]
fix = stmod$parlist$b_j[20]

scallopscale = scallop_index$est[1]
lobsterscale = lobster_index$est[1]
# For each simulated data set, fit different models and make indices
future::plan(multisession, workers = 10)
# Indices based on a negative binomial model
set.seed(8)
sim_nb_index_list <- future_lapply(sim_nb_list, get_pars_and_index, 
                                   mesh = barrier_mesh, grid = sim_grid_yrs, 
                                   area = 4, model = "combined",
                                   scale = scale, bc = TRUE,
                                   future.seed = TRUE)
# Scallop survey locations only
set.seed(8)
sim_scallop_index_list <- future_lapply(scallop_list, get_pars_and_index, 
                                        mesh = scallop_barrier_mesh, 
                                        grid = sim_grid_yrs[sim_grid_yrs$year != 2020, ], 
                                        area = 4, model = "scallop",
                                        scale = scale, fix = fix, bc = TRUE,
                                        future.seed = TRUE)
# Lobster survey locations only
set.seed(8)
sim_lobster_index_list <- future_lapply(lobster_list, get_pars_and_index, 
                                        mesh = lobster_barrier_mesh, 
                                        grid = sim_grid_yrs, 
                                        area = 4, model = "lobster",
                                        scale = scale, bc = TRUE,
                                        future.seed = TRUE)
future::plan(sequential)

# Plot the simulated indices against the true index for each scenario
p_int <- plot_sim_indices(st_index_scaled$est, sim_nb_index_list)+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p_lob <- plot_sim_indices(st_index_scaled$est, sim_lobster_index_list)+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous(limits = c(0, 4.5))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p_sca <- plot_sim_indices(c(st_index_scaled$est[1:15],NA,st_index_scaled$est[17:19]),
                          sim_scallop_index_list)+
  scale_x_continuous(breaks = seq(2005, 2023, by = 5))+
  scale_y_continuous(limits = c(0, 4.5))

cowplot::plot_grid(p_int, p_lob, p_sca, labels = c('(a)', '(b)', '(c)'), label_size = 12, ncol = 1)

# Try comparing against the scallop and lobster "true" indices
plot_sim_indices(c(scallop_index_scaled_single$est[1:15],NA,scallop_index_scaled_single$est[17:19]),
                 sim_scallop_index_list_single)
plot_sim_indices(lobster_index_scaled_single$est, sim_lobster_index_list_single)

# Get relative errors for each scenario
nb_rel_errors <- get_rel_errors(st_index_scaled$est, 
                                sim_nb_index_list, n = nsims)
scallop_rel_errors <- get_rel_errors(st_index_scaled$est, 
                                     sim_scallop_index_list, n=nsims)
lobster_rel_errors <- get_rel_errors(st_index_scaled$est, 
                                     sim_lobster_index_list, n = nsims)

# Calculate Mean Relative Errors for each year
nb_MREs <- round(rowMeans(nb_rel_errors), 3)
scallop_MREs <- round(rowMeans(scallop_rel_errors), 3)
lobster_MREs <- round(rowMeans(lobster_rel_errors), 3)
MREs <- data.frame(year = 2005:2023, nb = nb_MREs, 
                   scallop = scallop_MREs, lobster = lobster_MREs)

# Calculate Root Mean Square Errors for each year
nb_RMSEs <- get_RMSEs(st_index_scaled$est, 
                      sim_nb_index_list, n = nsims)
scallop_RMSEs <- get_RMSEs(st_index_scaled$est, 
                           sim_scallop_index_list, n = nsims)
lobster_RMSEs <- get_RMSEs(st_index_scaled$est, 
                           sim_lobster_index_list, n = nsims)
RMSEs <- data.frame(year = 2005:2023, nb = nb_RMSEs, 
                    scallop = scallop_RMSEs, lobster = lobster_RMSEs)

# Calculate 95% CI coverage rates for each year
nb_CRs <- get_CovRates(st_index_scaled$est, 
                       sim_nb_index_list, n = nsims)
scallop_CRs <- get_CovRates(st_index_scaled$est,
                            sim_scallop_index_list, n = nsims)
lobster_CRs <- get_CovRates(st_index_scaled$est, 
                            sim_lobster_index_list, n = nsims)
CRs <- data.frame(year = 2005:2023, nb = nb_CRs,
                  scallop = scallop_CRs, lobster = lobster_CRs)

# Calculate 95% CI average widths for each year
nb_widths <- get_AveWidths(sim_nb_index_list, n = nsims)
scallop_widths <- get_AveWidths(sim_scallop_index_list, n = nsims)
lobster_widths <- get_AveWidths(sim_lobster_index_list, n = nsims)
AveWidths <- data.frame(year = 2005:2023, nb = nb_widths, 
                        scallop = scallop_widths, lobster = lobster_widths)

# Make perfomance metric tables
kable(MREs, format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(RMSEs, format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(CRs, format = "latex", booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
kable(AveWidths, format = "latex", digits = 3, booktabs= TRUE,
      col.names = c("Year", "Both surveys",
                    "Scallop survey", "Lobster survey"))
# Plot performance metrics
# Begin by pivoting the tables of metrics to long format
MREs_long <- pivot_longer(MREs, cols = nb:lobster, 
                          names_to = "model", values_to = "MRE")
RMSEs_long <- pivot_longer(RMSEs, cols = nb:lobster, 
                           names_to = "model", values_to = "RMSE")
CRs_long <- pivot_longer(CRs, cols = nb:lobster, 
                         names_to = "model", values_to = "CR")
AWs_long <- pivot_longer(AveWidths, cols = nb:lobster, 
                         names_to = "model", values_to = "AW")
# Join them together
metrics <- MREs_long %>% full_join(RMSEs_long) %>% 
  full_join(CRs_long) %>% full_join(AWs_long) %>% 
  pivot_longer(cols = MRE:AW, names_to = "metric", values_to = "value") %>% 
  mutate(model = factor(model, levels = c("nb", "lobster", "scallop")), metric = as_factor(metric))
# Make labels
labels <- c(
  MRE = "Mean relative error",
  RMSE = "Root mean square error",
  CR = "Coverage rate",
  AW = "Average width"
)
# Plot for simulation study
ggplot(filter(metrics, model %in% c("nb", "scallop", "lobster")), aes(year, value, colour = model))+
  geom_line()+
  facet_wrap(facets = vars(metric), scales = "free_y",
             labeller = labeller(metric = labels))+
  scale_colour_discrete(name = "Model", labels = c("Integrated", "Lobster", "Scallop"))+
  labs(x = "Year", y = "Metric value")+
  theme_light()
