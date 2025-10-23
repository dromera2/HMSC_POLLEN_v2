library(targets)
library(tarchetypes)

source("packages.R")
source("functions.R")
source("arguments.R")

#As administrator, if ploblems, use "" (NOT RUN)
# Mklink /D dsclim X:/Data/dsclim_v1/netcdf
# Mklink /D soil X:/Data/SoilGrids

values <- expand_grid(i_target = i)

tvalues <- data.frame(start = starts, end = ends)

fvalues <- expand.grid(rcp = rcps, gcm = gcms, stringsAsFactors = FALSE)

tmapped <- tar_map(values = tvalues,
                   unlist = FALSE,
                   tar_target(map_env_data, prepare.env.data(start, end))
)

fmapped <- tar_map(values = fvalues,
                   unlist = FALSE,
                   tar_target(map_fut_data, prepare.future.data(41, 150, rcp, gcm))
)

mapped <- tar_map(values = values,
                  unlist = FALSE,
                  tar_target(input, prepare.data("over", gbif_data, complection_rast,
                                                 complection_data, 0.75, climate_data, soil_data,
                                                 TRUE, 500, i_target, myspecies)),
                  tar_target(slow_model, slow.calibrate.hmsc(input, thin, samples, transient,
                                                             nchains, myspecies, var_names,
                                                             coord_names, formula, nParallel)),
                  tar_target(single_model, single.hmsc.calibration(myspecies, input, thin, samples, transient, nchains, var_names,
                                                                   coord_names, formula, nParallel)),
                  tar_target(diagnostics, plot.diagnostics(slow_model, "over", i_target)),
                  tar_target(variance, computeVariancePartitioning(slow_model)),
                  tar_target(grouped_variance, computeVariancePartitioning(slow_model,
                                                                           group = group,
                                                                           groupnames = groupnames)),
                  tar_target(betas, get.betas(slow_model, myspecies, "over", i_target, var_names)),
                  tar_target(omegas, get.omegas(slow_model, myspecies)),
                  tar_target(omegas_df, omegas_to_df(omegas, myspecies, "over", i_target)),
                  tar_target(spatial_pred, spatial.prediction(climate_data, soil_data, myspecies, slow_model,
                                                              i_target, "over")),
                  tar_target(single_spatial_pred, single.spatial.prediction(climate_data, soil_data, myspecies, single_model,
                                                              i_target, "over")),
                  tar_target(slow_present_exp_eval, tuning.metrics(input, input, slow_model, myspecies, Yc=NULL, i_target, var_names)),
                  tar_target(single_slow_present_exp_eval, single.tuning.metrics(myspecies, single_model, input, i_target, var_names)),
                  tar_target(slow_present_pred_eval, tuning.predictive(input, slow_model, myspecies, i_target)),
                  tar_target(single_slow_present_pred_eval,  single.tuning.predictive(myspecies, single_model, input, i_target)),
                  #######
                  tar_target(timeseries_pred, timeseries.prediction(timeseries_data, var_names, coord_names, myspecies, slow_model, i_target, "over")),
                  tar_target(spatial_past, past.projection(env_data, soil_data, myspecies, slow_model, CM = TRUE, 160, i_target)),
                  tar_target(spatial_past_LGM, past.projection(env_data, soil_data, myspecies, slow_model, CM = TRUE, 1, i_target)),
                  tar_target(single_timeseries_pred, single.timeseries.prediction(timeseries_data, var_names, coord_names, myspecies,
                                                                                  single_model, i_target, "over")),
                  tar_target(single_past, single.past.projection(myspecies, env_data, soil_data, single_model, CM = TRUE, period = 160, i = i_target)),
                  tar_target(single_past_LGM, single.past.projection(myspecies, env_data, soil_data, single_model, CM = TRUE, period = 1, i = i_target)),
                  tar_target(past_pred_22_1ka, past.projection.loop(env_data, soil_data, myspecies, slow_model, seq(1, 220, by = 10), i_target)),
                  tar_target(single_past_pred_22_1ka, single.past.projection.loop(env_data, soil_data, myspecies, single_model, seq(1, 220, by = 10), i_target))
)

pvalues <- expand.grid(j = i, rcp = rcps, gcm = gcms, stringsAsFactors = FALSE) %>% 
  mutate(data = rlang::syms(paste("map_fut_data", rcp, str_replace_all(gcm, "-", "."), sep = "_"))) %>%
  mutate(model = rlang::syms(paste("slow_model", j, sep = "_")))

pmapped <- tar_map(values = pvalues,
                   unlist = FALSE,
                   names = c(j, rcp, gcm),
                   tar_target(fut_pred,
                              fut.projection(data,
                                             soil_data, 
                                             myspecies, 
                                             model,
                                             81:110, 
                                             j, 
                                             rcp, 
                                             gcm))
)


single_pvalues <- expand.grid(j = i, rcp = rcps, gcm = gcms, stringsAsFactors = FALSE) %>% 
  mutate(data = rlang::syms(paste("map_fut_data", rcp, str_replace_all(gcm, "-", "."), sep = "_"))) %>%
  mutate(model = rlang::syms(paste("single_model", j, sep = "_")))

single_pmapped <- tar_map(values = single_pvalues,
                   unlist = FALSE,
                   names = c(j, rcp, gcm),
                   tar_target(single_fut_pred,
                              single.fut.projection(data,
                                             soil_data, 
                                             myspecies, 
                                             model,
                                             81:110, 
                                             j, 
                                             rcp, 
                                             gcm))
)



list(
  tar_target(indf_data, bdvis::format_bdvis(gbif_data, config = conf)),
  tar_target(complection_data, completeness(indf_data, mask)),
  tar_target(complection_rast, rasterize_completeness(complection_data, "c", mask), format="qs"),
  tar_target(sprichness_rast, rasterize_completeness(complection_data, "Sobs", mask), format="qs"),
  tar_target(climate_data, get.climate.data(bioclimate, envirem, var_list), format="qs"),
  tar_target(soil_data, get.soil.data(pH, sand, mask), format="qs"),
  tar_target(i_target, i),
  tar_target(site_data, prepare.site.data(neotoma, climate_data, soil_data)),
  tmapped,
  tar_combine(
    env_data, tmapped[["map_env_data"]], command = list(!!!.x)
  ),
  tar_target(timeseries_data, prepare.timeseries.data(neotoma, env_data, soil_data)),
  mapped,
  tar_combine(
    all_grouped_var, mapped[["grouped_variance"]], command = list(!!!.x)
  ),
  tar_combine(
    all_var, mapped[["variance"]], command = list(!!!.x)
  ),
  tar_target(grouped_var_plot, plot.variance(all_grouped_var, i_target, 1:3, "grouped")),
  tar_target(var_plot, plot.variance(all_var, i_target, 1:8, "all")),
  tar_combine(
    omegas_all, mapped[["omegas_df"]]
  ),
  tar_target(omegas_plot, plot.omegas(omegas_all, myspecies)),
  tar_combine(
    all_pred, mapped[["spatial_pred"]]
  ),
  tar_target(spatial_map, prediction.map(all_pred, myspecies, "over", "spatial_pred")),
  tar_combine(
    single_all_pred, mapped[["single_spatial_pred"]]
  ),
  tar_target(single_spatial_map, prediction.map(single_all_pred, myspecies, "over", "single_spatial_pred")),
  tar_combine(all_slow_present_exp_eval, mapped[["slow_present_exp_eval"]]),
  tar_combine(all_single_slow_present_exp_eval, mapped[["single_slow_present_exp_eval"]]),
  tar_combine(all_slow_present_pred_eval, mapped[["slow_present_pred_eval"]]),
  tar_combine(all_single_slow_present_pred_eval, mapped[["single_slow_present_pred_eval"]]),
  tar_target(exp_present, present.eval(all_slow_present_exp_eval, all_single_slow_present_exp_eval, "present_exp")),
  tar_target(pred_present, present.eval(all_slow_present_pred_eval, all_single_slow_present_pred_eval, "present_pred")),
  tar_target(combined_eval, combine.eval(exp_present, pred_present)),
  fmapped,
  pmapped,
  single_pmapped,
  tar_combine(
    all_timeseries_pred, mapped[["timeseries_pred"]]),
  tar_combine(
    all_past_pred, mapped[["spatial_past"]]),
   tar_target(CLM_eval, neotoma.validation(neotoma_2_eval, mask, all_timeseries_pred)),
  tar_combine(
    all_single_timeseries_pred, mapped[["single_timeseries_pred"]]),
  tar_target(SDM_eval, neotoma.validation(neotoma_2_eval, mask, all_single_timeseries_pred)),
  tar_target(past_evaluation, eval.past(CLM_eval, SDM_eval)),
  tar_combine(
    all_single_past_pred, mapped[["single_past"]]),
  tar_combine(
    all_past_pred_LGM, mapped[["spatial_past_LGM"]]),
  tar_combine(
    all_single_past_pred_LGM, mapped[["single_past_LGM"]]),
  tar_target(th, get.threshold(all_pred, single_all_pred, gbif_data, myspecies, climate_data)), 
  tar_target(projection_map, past.map(all_past_pred, all_single_past_pred, gbif_data, myspecies, climate_data, "past_proj_MH", th)),
  tar_target(projection_map_LGM, past.map(all_past_pred_LGM, all_single_past_pred_LGM, gbif_data, myspecies, climate_data, "past_proj_LGM", th)),
  tar_target(present_map, present.map(all_pred, single_all_pred, gbif_data, myspecies, climate_data, "present_proj", th)),
  tar_target(refugia, refugia.map(present_map, projection_map, projection_map_LGM, "refugia")),
  ####another approach using prevalence as threshold
  # tar_target(prev_projection_map, past.map(all_past_pred, all_single_past_pred, gbif_data, myspecies, climate_data, "prev_past_proj_MH", prev_th)),
  # tar_target(prev_projection_map_LGM, past.map(all_past_pred_LGM, all_single_past_pred_LGM, gbif_data, myspecies, climate_data, "prev_past_proj_LGM", prev_th)),
  # tar_target(prev_present_map, present.map(all_pred, single_all_pred, gbif_data, myspecies, climate_data, "prev_present_proj", prev_th)),
  # tar_target(prev_refugia, refugia.map(prev_present_map, prev_projection_map, prev_projection_map_LGM, "prev_refugia")),
  ######future predictions
  tar_combine(all_fut_pred, pmapped["fut_pred"]),
  tar_combine(all_single_fut_pred, single_pmapped["single_fut_pred"]),
  tar_target(fut_map, future.map(all_fut_pred, "fut_proj")),
  tar_target(single_fut_map, future.map(all_single_fut_pred, "single_fut_proj")),
  tar_target(fut_subset, subset.future(fut_map, "subset_fut_v2")),
  tar_target(single_fut_subset, subset.future(single_fut_map, "single_subset_fut_v2")),
  #########continuous past projections
  tar_combine(all_past_pred_22_1ka, mapped[["past_pred_22_1ka"]], command = list(!!!.x)),
  tar_target(past_proj_22_1ka, do.call(bind_rows, all_past_pred_22_1ka)),
  tar_combine(all_single_past_pred_22_1ka, mapped[["single_past_pred_22_1ka"]], command = list(!!!.x)),
  tar_target(single_past_proj_22_1ka, do.call(bind_rows, all_single_past_pred_22_1ka)),
  tar_target(past_slope, map.slope(past_proj_22_1ka, CM = TRUE, th, myspecies, "masked_past_slope")),
  tar_target(single_past_slope, map.slope(single_past_proj_22_1ka, CM = FALSE, th, myspecies, "masked_single_past_slope")),
  ########composite plot
  tar_target(arrange_figure, plot.arrange(neotoma_2_eval, mask, all_timeseries_pred, complection_rast, indf_data)),
  tar_target(diff_figure, plot.diff(fut_subset, single_fut_subset, th, myspecies)),
  tar_target(combined_slope, combine.slope.map(past_slope, single_past_slope)),
  tar_target(sup_diff_figure, plot.diff.sup(fut_map, single_fut_map, th, myspecies))
  )
