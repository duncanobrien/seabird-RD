##########################################################################################
#### Preamble ####
##########################################################################################

set.seed(863)

library(readxl)
library(dplyr)
library(tidyr)
library(janitor)
library(poptrend)

lookup_vec <- c(
  year = "krykkje",
  year = "ar",
  year = "a_r",
  est_count = "antall",
  est_count = "antall_kyst",
  colony = "keysite",
  colony = "superloc_keysite",
  species = "artsnavn_no_art"
)

model_formula <-  est_count ~ trend(year, type = "index", tempRE = FALSE) + plot_id

#model_formula <-  est_count ~ trend(year, type = "index", tempRE = FALSE) + s(plot_id, bs = "re")

##########################################################################################
#### Load abundance data to be interpolated ####
##########################################################################################

abundances <- read.csv("Data/Abundances/SEAPOP_Abundances_forMultiStress_noDupl.csv") |>
  janitor::clean_names() |>
  dplyr::rename(dplyr::any_of(lookup_vec)) |>
  dplyr::mutate(
    species = dplyr::case_when(
      species == "Lunde" ~ "AtPu",
      species == "Krykkje" ~ "BlKi",
      TRUE ~ NA
    )
  ) |>
  subset((interaction(species, colony) %in% c("AtPu.Hjelmsøya","AtPu.Hornøya","BlKi.Hornøya"))) |>
  subset(!(interaction(species, colony) %in% c("AtPu.Hjelmsøya") & year <2015)) |> #drop time period of Hornoya puffins prior to standardised sampling
  dplyr::mutate(
    plot_id = interaction(
      loknr_lok,
      feltnr_lok,
      drop = TRUE,
      sep = "_"
    )
  )  |>
  dplyr::filter(!is.na(est_count)) |>
  dplyr::mutate(pop_id = paste0(interaction(species,colony))) 

##########################################################################################
#### Fit interpolation models ####
##########################################################################################

abundances_interp <- lapply(unique(abundances$pop_id),function(i){
  
  mod_data <- subset(abundances, pop_id == i) |>
    dplyr::mutate(plot_id = droplevels(plot_id))
  
  est_index <- ptrend(
    formula = model_formula,
    data = mod_data,
    nBoot = 1000,
    gamModel = TRUE
  )
  
  time_factor <- est_index$timeVarFac
  year_levels <- levels(est_index$gam$model[[time_factor]])
  plot_levels <- levels(mod_data$plot_id)
  
  prediction_grid <- expand_grid(
    year = as.numeric(year_levels),
    plot_id = plot_levels
  ) |>
    dplyr::mutate(plot_id = factor(plot_id, levels = plot_levels))
  
  prediction_grid[[time_factor]] <- factor(
    prediction_grid$year,
    levels = year_levels
  )
  
  prediction_grid$predicted_count <- as.numeric(
    predict(
      est_index$gam,
      newdata = prediction_grid,
      type = "response"
    )
  )
  
  #### Retain observations and impute missing counts
  imputed_counts <- prediction_grid |>
    dplyr::left_join(
      mod_data |> dplyr::select(year, plot_id, observed_count = est_count),
      by = c("year", "plot_id")
    ) |>
    dplyr::mutate(
      was_imputed = is.na(observed_count),
      imputed_count = dplyr::coalesce(observed_count, predicted_count),
      .before = 1
    ) |>
    dplyr::arrange(year, plot_id)
  
  #### Sum across plots
  annual_totals <- imputed_counts |>
    dplyr::reframe(observed_total = sum(observed_count, na.rm = TRUE),
                    total_imputed = sum(imputed_count),
                    total_fitted = sum(predicted_count),
                    n_plots = dplyr::n(),
                    n_plots_observed = sum(!was_imputed),
                    n_plots_imputed = sum(was_imputed),
                    .by = c(year)) |>
    dplyr::mutate(pop_id = i,.before = 1) |>
    tidyr::separate_wider_delim(pop_id,delim = ".",names = c("species","colony"))

  return(annual_totals)
}) |>
  do.call("rbind",args = _) |>
  #dplyr::mutate(total_imputed = round(total_imputed, digits = 0)) |> #convert to integer
  dplyr::select(species,colony,year,total_imputed)

##########################################################################################
#### Fit interpolation models (Rost kittiwakes) ####
##########################################################################################

## Define plots to use
vedoy_plots <- c("k_i", "k_ii",  "k_iv", "k_v", "k_xii") 
harbour_plots <- c("k_viii", "k_viii_2", "skrent_a", "skrent_b","k_ix", "k_x", "greger", "tyvsoy")

abundance_rost_blki <- read_excel("Data/Abundances/Krykkje overvåkingsdata Røst til Kate.xlsx") |>
  janitor::row_to_names(row_number = 2) |>
  janitor::clean_names() |>
  dplyr::slice(-1) |>
  dplyr::rename(year = na) |>
  dplyr::select(year, dplyr::all_of(c(vedoy_plots, harbour_plots))) |>
  dplyr::mutate(across(dplyr::everything(), ~ as.numeric(.x))) |>
  dplyr::filter(!is.na(year))|>
  subset(year >= 1988)

rost_long <- abundance_rost_blki |>
  tidyr::pivot_longer(
    cols = -year,
    names_to = "plot_id",
    values_to = "est_count"
  ) |>
  dplyr::mutate(
    sub_locality = if_else(
      plot_id %in% vedoy_plots,
      "vedoy",
      "harbour"
    )
  ) |>
  dplyr::filter(!is.na(est_count))

rost_interp <- lapply(unique(rost_long$sub_locality),function(i){
  
  mod_data <- subset(rost_long, sub_locality == i) |>
    dplyr::mutate(plot_id = droplevels(factor(plot_id)))
  
  est_index <- ptrend(
    formula = model_formula,
    data = mod_data |>
      {
        function(x) {
          if(i == "vedoy"){
            x |>
              subset(year <2020) #subset vedoy prior to 2020 as extinct from 2020 onwards
          }else{
            x
          }
        }}(),
    nBoot = 1000,
    gamModel = TRUE
  )
  
  time_factor <- est_index$timeVarFac
  year_levels <- levels(est_index$gam$model[[time_factor]])
  plot_levels <- levels(mod_data$plot_id)
  
  prediction_grid <- expand_grid(
    year = as.numeric(year_levels),
    plot_id = plot_levels
  ) |>
    dplyr::mutate(plot_id = factor(plot_id, levels = plot_levels))
  
  prediction_grid[[time_factor]] <- factor(
    prediction_grid$year,
    levels = year_levels
  )
  
  prediction_grid$predicted_count <- as.numeric(
    predict(
      est_index$gam,
      newdata = prediction_grid,
      type = "response"
    )
  )
  
  #### Retain observations and impute missing counts
  imputed_counts <- prediction_grid |>
    dplyr::left_join(
      mod_data |> dplyr::select(year, plot_id, observed_count = est_count),
      by = c("year", "plot_id")
    ) |>
    dplyr::mutate(
      was_imputed = is.na(observed_count),
      imputed_count = dplyr::coalesce(observed_count, predicted_count),
      .before = 1
    ) |>
    dplyr::arrange(year, plot_id)
  
  #### Sum across plots
  annual_totals <- imputed_counts |>
    dplyr::reframe( observed_total = sum(observed_count, na.rm = TRUE),
                    total_imputed = sum(imputed_count),
                    total_fitted = sum(predicted_count),
                    n_plots = dplyr::n(),
                    n_plots_observed = sum(!was_imputed),
                    n_plots_imputed = sum(was_imputed),
                    .by = c(year)) |>
    dplyr::mutate(sub_locality = i)
  
  return(annual_totals)
}) |>
  do.call("rbind",args = _) |>
  # Retain entirely unobserved years as NA, rather than zero.
  tidyr::complete(
    year = sort(unique(rost_long$year)),
    sub_locality = c("vedoy", "harbour")
  ) |>
  dplyr::mutate(dplyr::across(c(total_imputed,total_fitted), ~ dplyr::if_else(sub_locality == "vedoy" & is.na(.x),0,.x))) |>
  dplyr::select(year,sub_locality,total_imputed, n_plots_observed, n_plots_imputed) |>
  tidyr::pivot_wider(
    names_from = sub_locality,
    values_from = c(total_imputed, n_plots_observed, n_plots_imputed),
    names_glue = "{sub_locality}_{.value}"
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(n_vedoy = length(vedoy_plots),
                n_harbour = length(harbour_plots)) |>
  dplyr::mutate(
    total = sum(vedoy_total_imputed, harbour_total_imputed),
    n_total = sum(n_vedoy, n_harbour)
  ) |>
  dplyr::mutate(ref_sample = 17296 / 408) |> #true abundance divided by sampled abundance on Vedoy
  dplyr::mutate(
    vedoy_weight = (vedoy_total_imputed * (
      ref_sample * n_vedoy / ((ref_sample * n_vedoy) + n_harbour)
    )),
    harbour_weight = (harbour_total_imputed * (
      n_harbour / ((ref_sample * n_vedoy) + n_harbour)
    )),
    total_imputed = vedoy_weight + harbour_weight,
    nsites = sum(n_vedoy, n_harbour)
  ) |>
  dplyr::mutate(vedoy_contr = vedoy_weight / total_imputed,
                harbour_contr = harbour_weight / total_imputed) |>
  dplyr::ungroup() |>
  #dplyr::mutate(total_imputed = round(total_imputed, digits = 0)) |> #convert weighted average to integer
  dplyr::mutate(colony = "Røst", species = "BlKi") |>
  dplyr::select(species,colony,year,total_imputed)

##########################################################################################
#### Merge interpolations and save out ####
##########################################################################################

new_counts <- dplyr::bind_rows(abundances_interp,rost_interp)

write.csv(
  new_counts,
  "data/Abundances/newpops_poptrend_annual_imputed_totals.csv",
  row.names = FALSE
)
