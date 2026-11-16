##########################################################################################
#### Preamble ####
##########################################################################################

library(dplyr)
library(brms)
library(cmdstanr)
library(ggplot2)
library(tidybayes)
library(ggdist)
library(patchwork)

set.seed(333)

source("R/H_clean_data.R")
source("R/utilities/respDIV.R")
source("R/utilities/respDIVme.R")

model_files <- fs::dir_ls(path = "outputs/models", glob = "*.RDS") #list all data files

##########################################################################################
#### Prepare data ####
##########################################################################################

RD_data_behav <- dplyr::left_join(energetics_behav_prop,
                                  sst_daily |>
                                    subset(season == "breeding"), by = c("species","colony", "season","date","year"),
                                  relationship = "many-to-many") |>
  tidyr::drop_na(dplyr::contains("std")) |> #need to drop first to allow proper ranking
  dplyr::mutate(bird1 = as.numeric(gsub("\\D", "", individ_id)))   |>
  dplyr::mutate(bird1 = paste(species,colony, individ_id,sep = "_")) |>
  dplyr::mutate(julian = as.numeric(format(date, "%j"))) |>
  dplyr::arrange(julian) |>
  dplyr::mutate(day1 = base::match(julian, sort(unique(julian)))) |> #rank julian days
  dplyr::mutate(year1 = base::match(year, sort(unique(year)))) |> #rank years
  dplyr::arrange(bird1,year, date) |>
  dplyr::mutate(pop_id = paste(species,colony,sep = "_")) |>
  dplyr::rename(sst = mean_sst) #must rename to not interfere with marginaleffects internal naming conventions

RD_data_energy <- energetics |>
  subset(DEEkJ_sd >0) |>
  dplyr::mutate(bird1 = as.numeric(gsub("\\D", "", individ_id)))   |>
  dplyr::mutate(bird1 = paste(species,colony, individ_id,sep = "_")) |>
  dplyr::mutate(julian = as.numeric(format(date, "%j"))) |>
  dplyr::arrange(julian) |>
  dplyr::mutate(day1 = base::match(julian, sort(unique(julian)))) |> #rank julian days
  dplyr::mutate(year1 = base::match(year, sort(unique(year)))) |> #rank years
  dplyr::arrange(bird1,year1, date) |>
  dplyr::mutate(pop_id = paste(species,colony,sep = "_")) |>
  dplyr::rename(sst = mean_sst) #must rename to not interfere with marginaleffects internal naming conventions

RD_data_abun <- dplyr::left_join(abundances |>
                                   dplyr::select(species,colony,year,tot_count),
                                 sst_yr |>
                                   subset(season == "breeding"), by = c("species","colony","year")) |>
  dplyr::rename(sst = mean_sst) |>
  #dplyr::mutate(tot_count = log10(tot_count)) |>
  #dplyr::mutate(std_count = tot_count/nsites) |>
  dplyr::left_join(RD_data_behav |>
                     dplyr::select(species, colony, year, year1) |>
                     dplyr::distinct(),by = c("species","colony", "year"))

##########################################################################################
#### Prepare behaviourial data ####
##########################################################################################

RD_behav_atpu <- RD_data_behav |>
  subset(species == "AtPu") |>
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, tForage_std) |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

RD_behav_blki <- RD_data_behav |>
  subset(species == "BlKi") |>
  subset(!colony %in% c("Isfjorden", "Hjelmsoya")) |>
  subset(!bird1 %in% c("BlKi_Anda_NOS_6198759","BlKi_Anda_NOS_6198952", "BlKi_Hornoya_NOS_6222268", "BlKi_Sklinna_NOS_6211517")) |> #drop birds that have unexplained periods of zero foraging across all conditions
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, tForage_std) |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

mod_atpu_behav <- lapply(grep("atpu_behav",model_files, value = TRUE),readRDS) |>
  lapply(X = _,`[[`,"model")
mod_blki_behav <- lapply(grep("blki_behav",model_files, value = TRUE),readRDS) |>
  lapply(X = _,`[[`,"model")

behav_stability_blki <- lapply(seq_along(mod_blki_behav),function(x){

  abun_data <- RD_data_abun |>
    subset(colony == unique(RD_behav_blki[[x]]$colony) & species == unique(RD_behav_blki[[x]]$species)) |>
    dplyr::rename(sst_raw = sst) |>
    dplyr::mutate(sst = (sst_raw - unique(RD_behav_blki[[x]]$mean_sst))/unique(RD_behav_blki[[x]]$sd_sst)) |>
    #na.omit() |>
    dplyr::arrange(year)
  
  newdata_yr <- expand.grid("sst" = abun_data$sst,
                            "bird1" = as.character(unique(mod_blki_behav[[x]]$data[["bird1"]]))
                            )
  
  rd_blki_behav_yr <- respDiv(mod_blki_behav[[x]], newdata_yr, gradient = "sst", unit = "bird1", resp = "tForage_std", ndraws = 1000, type = "link", raw = TRUE)
  
  detrend_abun <- resid(lm(abun_data$tot_count ~ abun_data$year))
  
  list("summary" = data.frame(
    "dev" = EWSmethods::uniJI(abun_data[, c("year", "tot_count")], winsize = 100, scale = TRUE)$smap_J,
    "dev_detrend" = EWSmethods::uniJI(cbind(abun_data$year,detrend_abun), winsize = 100, scale = TRUE)$smap_J,
    "ts" = mean(abun_data$tot_count) / sd(abun_data$tot_count),
    "ts_detrend" = mean(abun_data$tot_count) / sd(detrend_abun),
    "median_dissim" = median(rd_blki_behav_yr$dissim$.dissimilarity),
    "median_diverg" = median(rd_blki_behav_yr$diverg$.divergence),
    "median_imbalance" = median(rd_blki_behav_yr$imbalance$.imbalance),
    "max_dissim" = max(rd_blki_behav_yr$dissim$.dissimilarity),
    "max_diverg" = max(rd_blki_behav_yr$diverg$.divergence),
    "max_imbalance" = max(rd_blki_behav_yr$imbalance$.imbalance),
    "sd_dissim" = sd(rd_blki_behav_yr$dissim$.dissimilarity),
    "sd_diverg" = sd(rd_blki_behav_yr$diverg$.divergence),
    "sd_imbalance" = sd(rd_blki_behav_yr$imbalance$.imbalance)
  ) |> 
    dplyr::mutate(species = unique(RD_behav_blki[[x]]$species),
                  colony = unique(RD_behav_blki[[x]]$colony)),
  "raw" = rd_blki_behav_yr |>
    lapply(function(.x){
      .x |>
        dplyr::mutate(species = unique(RD_behav_blki[[x]]$species),
                      colony = unique(RD_behav_blki[[x]]$colony))
    }))
}) 

behav_stability_atpu <- lapply(seq_along(mod_atpu_behav),function(x){
  
  abun_data <- RD_data_abun |>
    subset(colony == unique(RD_behav_atpu[[x]]$colony) & species == unique(RD_behav_atpu[[x]]$species)) |>
    dplyr::rename(sst_raw = sst) |>
    dplyr::mutate(sst = (sst_raw - unique(RD_behav_atpu[[x]]$mean_sst))/unique(RD_behav_atpu[[x]]$sd_sst))  |>
    #na.omit() |>
    dplyr::arrange(year)
  
  newdata_yr <- expand.grid("sst" = abun_data$sst,
                            "bird1" = as.character(unique(mod_atpu_behav[[x]]$data[["bird1"]]))
                            )
  
  rd_atpu_behav_yr <- respDiv(mod_atpu_behav[[x]], newdata_yr, gradient = "sst", unit = "bird1", resp = "tForage_std", ndraws = 1000, type = "link", raw = TRUE)
  
  detrend_abun <- resid(lm(abun_data$tot_count ~ abun_data$year))
  
  list("summary" = data.frame(
    "dev" = EWSmethods::uniJI(abun_data[, c("year", "tot_count")], winsize = 100, scale = TRUE)$smap_J,
    "dev_detrend" = EWSmethods::uniJI(cbind(abun_data$year,detrend_abun), winsize = 100, scale = TRUE)$smap_J,
    "ts" = mean(abun_data$tot_count) / sd(abun_data$tot_count),
    "ts_detrend" = mean(abun_data$tot_count) / sd(detrend_abun),
    "median_dissim" = median(rd_atpu_behav_yr$dissim$.dissimilarity),
    "median_diverg" = median(rd_atpu_behav_yr$diverg$.divergence),
    "median_imbalance" = median(rd_atpu_behav_yr$imbalance$.imbalance),
    "max_dissim" = max(rd_atpu_behav_yr$dissim$.dissimilarity),
    "max_diverg" = max(rd_atpu_behav_yr$diverg$.divergence),
    "max_imbalance" = max(rd_atpu_behav_yr$imbalance$.imbalance),
    "sd_dissim" = sd(rd_atpu_behav_yr$dissim$.dissimilarity),
    "sd_diverg" = sd(rd_atpu_behav_yr$diverg$.divergence),
    "sd_imbalance" = sd(rd_atpu_behav_yr$imbalance$.imbalance)
  ) |> 
    dplyr::mutate(species = unique(RD_behav_atpu[[x]]$species),
                  colony = unique(RD_behav_atpu[[x]]$colony)),
  "raw" = rd_atpu_behav_yr |>
    lapply(function(.x){
      .x |>
        dplyr::mutate(species = unique(RD_behav_atpu[[x]]$species),
                      colony = unique(RD_behav_atpu[[x]]$colony))
    }))
}) 

behav_stability_abun_atpu <- behav_stability_atpu |>
  lapply(X = , `[[`, "summary") |>
  do.call("rbind", args = _)
behav_stability_abun_blki <- behav_stability_blki |>
  lapply(X = , `[[`, "summary") |>
  do.call("rbind", args = _)

##########################################################################################
#### Prepare DEE data ####
##########################################################################################

RD_dee_atpu <-  RD_data_energy |>
  subset(species == "AtPu") |>
  subset(DEEkJ < 1500) |>
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, DEEkJ, DEEkJ_sd) |>
  dplyr::filter(bird1 != "AtPu_Hjelmsoya_NOS_5175702") |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

RD_dee_blki <- RD_data_energy |>
  subset(species == "BlKi") |>
  subset(!colony %in% c("Isfjorden", "Hjelmsoya")) |>
  subset(!bird1 %in% c("BlKi_Anda_NOS_6198759","BlKi_Anda_NOS_6198952", "BlKi_Hornoya_NOS_6222268", "BlKi_Sklinna_NOS_6211517")) |> #drop birds that have unexplained periods of zero foraging across all conditions
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, DEEkJ, DEEkJ_sd) |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

mod_atpu_dee <- lapply(grep("atpu_energy",model_files, value = TRUE),readRDS) |>
  lapply(X = _,`[[`,"model")
mod_blki_dee <- lapply(grep("blki_energy",model_files, value = TRUE),readRDS) |>
  lapply(X = _,`[[`,"model")

dee_stability_blki <- lapply(seq_along(mod_blki_dee),function(x){
  
  abun_data <- RD_data_abun |>
    subset(colony == unique(RD_dee_blki[[x]]$colony) & species == unique(RD_dee_blki[[x]]$species)) |>
    dplyr::rename(sst_raw = sst) |>
    dplyr::mutate(sst = (sst_raw - unique(RD_dee_blki[[x]]$mean_sst))/unique(RD_dee_blki[[x]]$sd_sst)) |>
    #na.omit() |>
    dplyr::arrange(year)
  
  newdata_yr <- expand.grid("sst" = abun_data$sst,
                            "bird1" = as.character(unique(mod_blki_dee[[x]]$data[["bird1"]]))
                            ) |>
    dplyr::mutate(DEEkJ_sd = 200)
  
  rd_blki_dee_yr <- respDivme(mod_blki_dee[[x]], newdata_yr, gradient = "sst", unit = "bird1", resp = "DEEkJ", ndraws = 1000, type = "link", raw = TRUE)
  
  detrend_abun <- resid(lm(abun_data$tot_count ~ abun_data$year))
  
  list("summary" = data.frame(
    "dev" = EWSmethods::uniJI(abun_data[, c("year", "tot_count")], winsize = 100, scale = FALSE)$smap_J,
    "dev_detrend" = EWSmethods::uniJI(cbind(abun_data$year,detrend_abun), winsize = 100, scale = FALSE)$smap_J,
    "ts" = mean(abun_data$tot_count) / sd(abun_data$tot_count),
    "ts_detrend" = mean(abun_data$tot_count) / sd(detrend_abun),
    "median_dissim" = median(rd_blki_dee_yr$dissim$.dissimilarity),
    "median_diverg" = median(rd_blki_dee_yr$diverg$.divergence),
    "median_imbalance" = median(rd_blki_dee_yr$imbalance$.imbalance),
    "max_dissim" = max(rd_blki_dee_yr$dissim$.dissimilarity),
    "max_diverg" = max(rd_blki_dee_yr$diverg$.divergence),
    "max_imbalance" = max(rd_blki_dee_yr$imbalance$.imbalance),
    "sd_dissim" = sd(rd_blki_dee_yr$dissim$.dissimilarity),
    "sd_diverg" = sd(rd_blki_dee_yr$diverg$.divergence),
    "sd_imbalance" = sd(rd_blki_dee_yr$imbalance$.imbalance)
  ) |> 
    dplyr::mutate(species = unique(RD_dee_blki[[x]]$species),
                  colony = unique(RD_dee_blki[[x]]$colony)),
  "raw" = rd_blki_dee_yr |>
    lapply(function(.x){
      .x |>
        dplyr::mutate(species = unique(RD_dee_blki[[x]]$species),
                      colony = unique(RD_dee_blki[[x]]$colony))
    }))
}) 

dee_stability_atpu <- lapply(seq_along(mod_atpu_dee),function(x){
  
  abun_data <- RD_data_abun |>
    subset(colony == unique(RD_dee_atpu[[x]]$colony) & species == unique(RD_dee_atpu[[x]]$species)) |>
    dplyr::rename(sst_raw = sst) |>
    dplyr::mutate(sst = (sst_raw - unique(RD_dee_atpu[[x]]$mean_sst))/unique(RD_dee_atpu[[x]]$sd_sst))  |>
    #na.omit() |>
    dplyr::arrange(year)
  
  newdata_yr <- expand.grid("sst" = abun_data$sst,
                            "bird1" = as.character(unique(mod_atpu_dee[[x]]$data[["bird1"]]))
                            ) |>
    dplyr::mutate(DEEkJ_sd = 150)
  
  rd_atpu_dee_yr <- respDivme(mod_atpu_dee[[x]], newdata_yr, gradient = "sst", unit = "bird1", resp = "DEEkJ", ndraws = 1000, type = "link", raw = TRUE)
  
  detrend_abun <- resid(lm(abun_data$tot_count ~ abun_data$year))
  list("summary" = data.frame(
    "dev" = EWSmethods::uniJI(abun_data[, c("year", "tot_count")], winsize = 100, scale = FALSE)$smap_J,
    "dev_detrend" = EWSmethods::uniJI(cbind(abun_data$year,detrend_abun), winsize = 100, scale = FALSE)$smap_J,
    "ts" = mean(abun_data$tot_count) / sd(abun_data$tot_count),
    "ts_detrend" = mean(abun_data$tot_count) / sd(detrend_abun),
    "median_dissim" = median(rd_atpu_dee_yr$dissim$.dissimilarity),
    "median_diverg" = median(rd_atpu_dee_yr$diverg$.divergence),
    "median_imbalance" = median(rd_atpu_dee_yr$imbalance$.imbalance),
    "max_dissim" = max(rd_atpu_dee_yr$dissim$.dissimilarity),
    "max_diverg" = max(rd_atpu_dee_yr$diverg$.divergence),
    "max_imbalance" = max(rd_atpu_dee_yr$imbalance$.imbalance),
    "sd_dissim" = sd(rd_atpu_dee_yr$dissim$.dissimilarity),
    "sd_diverg" = sd(rd_atpu_dee_yr$diverg$.divergence),
    "sd_imbalance" = sd(rd_atpu_dee_yr$imbalance$.imbalance)
  ) |> 
    dplyr::mutate(species = unique(RD_dee_atpu[[x]]$species),
                  colony = unique(RD_dee_atpu[[x]]$colony)),
  "raw" = rd_atpu_dee_yr |>
    lapply(function(.x){
      .x |>
        dplyr::mutate(species = unique(RD_dee_atpu[[x]]$species),
                      colony = unique(RD_dee_atpu[[x]]$colony))
    }))
}) 

dee_stability_abun_atpu <- dee_stability_atpu |>
  lapply(X = , `[[`, "summary") |>
  do.call("rbind", args = _)
dee_stability_abun_blki <- dee_stability_blki |>
  lapply(X = , `[[`, "summary") |>
  do.call("rbind", args = _)

##########################################################################################
#### Behaviour imbalance ####
##########################################################################################

behav_imbalance_atpu <- lapply(behav_stability_atpu,`[[`,"raw")  |>
  lapply(`[[`,"imbalance") |>
  do.call("rbind",args = _) |>
  dplyr::reframe(.imbalance = median(.imbalance),.by = c(species,colony,drawid))
behav_imbalance_blki <- lapply(behav_stability_blki,`[[`,"raw") |>
  lapply(`[[`,"imbalance") |>
  do.call("rbind",args = _) |>
  dplyr::reframe(.imbalance = median(.imbalance),.by = c(species,colony,drawid))

behav_imbalance_a <- lapply(1:1000, function(i) {
  mod_data <- subset(behav_imbalance_atpu, drawid == i) |>
    dplyr::left_join(behav_stability_abun_atpu, by = c("species","colony"))
  
  m_dev <- brms::brm(brms::bf(dev ~ .imbalance, family = "gaussian"),
                     data = mod_data,
                     prior = c(prior(normal(0, 1), class = "b"),
                               prior(normal(0, 1), class = "Intercept")
                               ,prior(exponential(10), class = "sigma")
                     ),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99),
                     backend = "cmdstanr",
                     silent = 2)
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .imbalance, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  m_ts <- brms::brm(brms::bf(ts ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  prior_sens <-  rbind(priorsense::powerscale_sensitivity(m_dev) |>
                         dplyr::mutate(response = "dev"),
                       priorsense::powerscale_sensitivity(m_dev_detrend) |>
                         dplyr::mutate(response = "dev_detrend"),
                       priorsense::powerscale_sensitivity(m_ts) |>
                         dplyr::mutate(response = "ts"),
                       priorsense::powerscale_sensitivity(m_ts_detrend) |>
                         dplyr::mutate(response = "ts_detrend")) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame("rho_dev" = cor(mod_data$.imbalance,mod_data$dev, method = "pearson"),
                            "rho_dev_detrend" = cor(mod_data$.imbalance,mod_data$dev_detrend, method = "pearson"),
                            "rho_ts" = cor(mod_data$.imbalance,mod_data$ts, method = "pearson"),
                            "rho_ts_detrend" = cor(mod_data$.imbalance,mod_data$ts_detrend, method = "pearson")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  out<- rbind(
    marginaleffects::avg_slopes(m_dev,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list("draws" = out,
              "rho" = correlation,
              "prior_sens" = prior_sens)
  )
})

behav_imbalance_b <- lapply(1:1000, function(i) {
  mod_data <- subset(behav_imbalance_blki, drawid == i) |>
    dplyr::left_join(behav_stability_abun_blki, by = c("species","colony"))
  
  m_dev <- brms::brm(brms::bf(dev ~ .imbalance, family = "gaussian"),
                     data = mod_data,
                     prior = c(prior(normal(0, 1), class = "b"),
                               prior(normal(0, 1), class = "Intercept")
                               ,prior(exponential(10), class = "sigma")
                     ),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99),
                     backend = "cmdstanr",
                     silent = 2)
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .imbalance, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  m_ts <- brms::brm(brms::bf(ts ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  prior_sens <-  rbind(priorsense::powerscale_sensitivity(m_dev) |>
                         dplyr::mutate(response = "dev"),
                       priorsense::powerscale_sensitivity(m_dev_detrend) |>
                         dplyr::mutate(response = "dev_detrend"),
                       priorsense::powerscale_sensitivity(m_ts) |>
                         dplyr::mutate(response = "ts"),
                       priorsense::powerscale_sensitivity(m_ts_detrend) |>
                         dplyr::mutate(response = "ts_detrend")) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame("rho_dev" = cor(mod_data$.imbalance,mod_data$dev, method = "pearson"),
                             "rho_dev_detrend" = cor(mod_data$.imbalance,mod_data$dev_detrend, method = "pearson"),
                            "rho_ts" = cor(mod_data$.imbalance,mod_data$ts, method = "pearson"),
                            "rho_ts_detrend" = cor(mod_data$.imbalance,mod_data$ts_detrend, method = "pearson")
                            ) |>
    dplyr::mutate(mc_sim = i)
  
  out<- rbind(
    marginaleffects::avg_slopes(m_dev,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list("draws" = out,
              "rho" = correlation,
              "prior_sens" = prior_sens)
  )
})

saveRDS(behav_imbalance_a,"outputs/stability/behav_imbalance_atpu.RDS")
saveRDS(behav_imbalance_b,"outputs/stability/behav_imbalance_blki.RDS")

##########################################################################################
#### Behaviour dissimilarity ####
##########################################################################################

behav_dissim_atpu <- lapply(behav_stability_atpu,`[[`,"raw")  |>
  lapply(`[[`,"dissim") |>
  do.call("rbind",args = _) |>
  dplyr::reframe(.dissimilarity = median(.dissimilarity),.by = c(species,colony,drawid))
behav_dissim_blki <- lapply(behav_stability_blki,`[[`,"raw")  |>
  lapply(`[[`,"dissim") |>
  do.call("rbind",args = _) |>
  dplyr::reframe(.dissimilarity = median(.dissimilarity),.by = c(species,colony,drawid))

behav_dissim_a <- lapply(1:1000, function(i) {
  mod_data <- subset(behav_dissim_atpu, drawid == i) |>
    dplyr::left_join(behav_stability_abun_atpu, by = c("species","colony"))
  
  m_dev <- brms::brm(brms::bf(dev ~ .dissimilarity, family = "gaussian"),
                     data = mod_data,
                     prior = c(prior(normal(0, 1), class = "b"),
                               prior(normal(0, 1), class = "Intercept")
                               ,prior(exponential(10), class = "sigma")
                     ),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99),
                     backend = "cmdstanr",
                     silent = 2)
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .dissimilarity, family = "gaussian"),
                     data = mod_data,
                     prior = c(prior(normal(0, 1), class = "b"),
                               prior(normal(0, 1), class = "Intercept")
                               ,prior(exponential(10), class = "sigma")
                     ),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99),
                     backend = "cmdstanr",
                     silent = 2)
  
  m_ts <- brms::brm(brms::bf(ts ~ .dissimilarity, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .dissimilarity, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  prior_sens <-  rbind(priorsense::powerscale_sensitivity(m_dev) |>
                         dplyr::mutate(response = "dev"),
                       priorsense::powerscale_sensitivity(m_dev_detrend) |>
                         dplyr::mutate(response = "dev_detrend"),
                       priorsense::powerscale_sensitivity(m_ts) |>
                         dplyr::mutate(response = "ts"),
                       priorsense::powerscale_sensitivity(m_ts_detrend) |>
                         dplyr::mutate(response = "ts_detrend")) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame("rho_dev" = cor(mod_data$.dissimilarity,mod_data$dev, method = "pearson"),
                             "rho_dev_detrend" = cor(mod_data$.dissimilarity,mod_data$dev_detrend, method = "pearson"),
                            "rho_ts" = cor(mod_data$.dissimilarity,mod_data$ts, method = "pearson"),
                            "rho_ts_detrend" = cor(mod_data$.dissimilarity,mod_data$ts_detrend, method = "pearson")) |>
    dplyr::mutate(mc_sim = i)
  
  out<- rbind(
    marginaleffects::avg_slopes(m_dev,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list("draws" = out,
              "rho" = correlation,
              "prior_sens" = prior_sens)
  )
})

behav_dissim_b <- lapply(1:1000, function(i) {
  mod_data <- subset(behav_dissim_blki, drawid == i)  |>
    dplyr::left_join(behav_stability_abun_blki, by = c("species","colony"))
  
  m_dev <- brms::brm(brms::bf(dev ~ .dissimilarity, family = "gaussian"),
                     data = mod_data,
                     prior = c(prior(normal(0, 1), class = "b"),
                               prior(normal(0, 1), class = "Intercept")
                               ,prior(exponential(10), class = "sigma")
                     ),
                     iter = 2000, warmup = 1000, chains = 4, cores = 4,
                     control = list(adapt_delta = 0.99),
                     backend = "cmdstanr",
                     silent = 2)
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .dissimilarity, family = "gaussian"),
                         data = mod_data,
                         prior = c(prior(normal(0, 1), class = "b"),
                                   prior(normal(0, 1), class = "Intercept"),
                                   prior(exponential(10), class = "sigma")),
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         control = list(adapt_delta = 0.99),
                         backend = "cmdstanr",
                         silent = 2)
  
  m_ts <- brms::brm(brms::bf(ts ~ .dissimilarity, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .dissimilarity, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  prior_sens <-  rbind(priorsense::powerscale_sensitivity(m_dev) |>
                         dplyr::mutate(response = "dev"),
                       priorsense::powerscale_sensitivity(m_dev_detrend) |>
                         dplyr::mutate(response = "dev_detrend"),
                       priorsense::powerscale_sensitivity(m_ts) |>
                         dplyr::mutate(response = "ts"),
                       priorsense::powerscale_sensitivity(m_ts_detrend) |>
                         dplyr::mutate(response = "ts_detrend")) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame("rho_dev" = cor(mod_data$.dissimilarity,mod_data$dev, method = "pearson"),
                             "rho_dev_detrend" = cor(mod_data$.dissimilarity,mod_data$dev_detrend, method = "pearson"),
                            "rho_ts" = cor(mod_data$.dissimilarity,mod_data$ts, method = "pearson"),
                            "rho_ts_detrend" = cor(mod_data$.dissimilarity,mod_data$ts_detrend, method = "pearson")) |>
    dplyr::mutate(mc_sim = i)
  
  out<- rbind(
    marginaleffects::avg_slopes(m_dev,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list("draws" = out,
              "rho" = correlation,
              "prior_sens" = prior_sens)
  )
})

saveRDS(behav_dissim_a,"outputs/stability/behav_dissim_atpu.RDS")
saveRDS(behav_dissim_b,"outputs/stability/behav_dissim_blki.RDS")

##########################################################################################
#### DEE imbalance ####
##########################################################################################

dee_imbalance_atpu <- lapply(dee_stability_atpu, `[[`, "raw")  |>
  lapply(`[[`, "imbalance") |>
  do.call("rbind", args = _) |>
  dplyr::reframe(.imbalance = median(.imbalance),
                 .by = c(species, colony, drawid))
dee_imbalance_blki <- lapply(dee_stability_blki, `[[`, "raw") |>
  lapply(`[[`, "imbalance") |>
  do.call("rbind", args = _) |>
  dplyr::reframe(.imbalance = median(.imbalance),
                 .by = c(species, colony, drawid))

dee_imbalance_a <- lapply(1:1000, function(i) {
  mod_data <- subset(dee_imbalance_atpu, drawid == i) |>
    dplyr::left_join(dee_stability_abun_atpu, by = c("species", "colony"))
  
  m_dev <- brms::brm(
    brms::bf(dev ~ .imbalance, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept")
      ,
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )

  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .imbalance, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)

  m_ts <- brms::brm(
    brms::bf(ts ~ .imbalance, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  prior_sens <-  rbind(
    priorsense::powerscale_sensitivity(m_dev) |>
      dplyr::mutate(response = "dev"),
    priorsense::powerscale_sensitivity(m_dev_detrend) |>
      dplyr::mutate(response = "dev_detrend"),
    priorsense::powerscale_sensitivity(m_ts) |>
      dplyr::mutate(response = "ts"),
    priorsense::powerscale_sensitivity(m_ts_detrend) |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame(
    "rho_dev" = cor(mod_data$.imbalance, mod_data$dev, method = "pearson"),
     "rho_dev_detrend" = cor(mod_data$.imbalance,mod_data$dev_detrend, method = "pearson"),
    "rho_ts" = cor(mod_data$.imbalance, mod_data$ts, method = "pearson"),
    "rho_ts_detrend" = cor(mod_data$.imbalance, mod_data$ts_detrend, method = "pearson")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  out <- rbind(
    marginaleffects::avg_slopes(m_dev, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list(
    "draws" = out,
    "rho" = correlation,
    "prior_sens" = prior_sens
  ))
})

dee_imbalance_b <- lapply(1:1000, function(i) {
  mod_data <- subset(dee_imbalance_blki, drawid == i)  |>
    dplyr::left_join(dee_stability_abun_blki, by = c("species", "colony"))
  
  m_dev <- brms::brm(
    brms::bf(dev ~ .imbalance, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept")
      ,
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .imbalance, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
 
  m_ts <- brms::brm(
    brms::bf(ts ~ .imbalance, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .imbalance, family = "gaussian"),
                    data = mod_data,
                    prior = c(prior(normal(0, 1), class = "b"),
                              prior(normal(0, 1), class = "Intercept"),
                              prior(exponential(10), class = "sigma")),
                    iter = 2000, warmup = 1000, chains = 4, cores = 4,
                    control = list(adapt_delta = 0.99),
                    backend = "cmdstanr",
                    silent = 2)
  
  prior_sens <-  rbind(
    priorsense::powerscale_sensitivity(m_dev) |>
      dplyr::mutate(response = "dev"),
    priorsense::powerscale_sensitivity(m_dev_detrend) |>
      dplyr::mutate(response = "dev_log"),
    priorsense::powerscale_sensitivity(m_ts) |>
      dplyr::mutate(response = "ts"),
    priorsense::powerscale_sensitivity(m_ts_detrend) |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame(
    "rho_dev" = cor(mod_data$.imbalance, mod_data$dev, method = "pearson"),
     "rho_dev_detrend" = cor(mod_data$.imbalance,mod_data$dev_detrend, method = "pearson"),
    "rho_ts" = cor(mod_data$.imbalance, mod_data$ts, method = "pearson"),
    "rho_ts_detrend" = cor(mod_data$.imbalance, mod_data$ts_detrend, method = "pearson")
    
  ) |>
    dplyr::mutate(mc_sim = i)
  
  out <- rbind(
    marginaleffects::avg_slopes(m_dev, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".imbalance") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend, variables = ".imbalance") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list(
    "draws" = out,
    "rho" = correlation,
    "prior_sens" = prior_sens
  ))
})

saveRDS(dee_imbalance_a,"outputs/stability/dee_imbalance_atpu.RDS")
saveRDS(dee_imbalance_b,"outputs/stability/dee_imbalance_blki.RDS")

##########################################################################################
#### DEE dissimilarity ####
##########################################################################################

dee_dissim_atpu <- lapply(dee_stability_atpu, `[[`, "raw")  |>
  lapply(`[[`, "dissim") |>
  do.call("rbind", args = _) |>
  dplyr::reframe(
    .dissimilarity = median(.dissimilarity),
    .by = c(species, colony, drawid)
  )
dee_dissim_blki <- lapply(dee_stability_blki, `[[`, "raw")  |>
  lapply(`[[`, "dissim") |>
  do.call("rbind", args = _) |>
  dplyr::reframe(
    .dissimilarity = median(.dissimilarity),
    .by = c(species, colony, drawid)
  )

dee_dissim_a <- lapply(1:1000, function(i) {
  mod_data <- subset(dee_dissim_atpu, drawid == i) |>
    dplyr::left_join(dee_stability_abun_atpu, by = c("species", "colony"))
  
  m_dev <- brms::brm(
    brms::bf(dev ~ .dissimilarity, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept")
      ,
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .dissimilarity, family = "gaussian"),
                         data = mod_data,
                         prior = c(prior(normal(0, 1), class = "b"),
                                   prior(normal(0, 1), class = "Intercept"),
                                   prior(exponential(10), class = "sigma")),
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         control = list(adapt_delta = 0.99),
                         backend = "cmdstanr",
                         silent = 2)
 
  m_ts <- brms::brm(
    brms::bf(ts ~ .dissimilarity, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .dissimilarity, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  prior_sens <-  rbind(
    priorsense::powerscale_sensitivity(m_dev) |>
      dplyr::mutate(response = "dev"),
    priorsense::powerscale_sensitivity(m_dev_detrend) |>
      dplyr::mutate(response = "dev_detrend"),
    priorsense::powerscale_sensitivity(m_ts) |>
      dplyr::mutate(response = "ts"),
    priorsense::powerscale_sensitivity(m_ts_detrend) |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame(
    "rho_dev" = cor(mod_data$.dissimilarity, mod_data$dev, method = "pearson"),
     "rho_dev_detrend" = cor(mod_data$.dissimilarity,mod_data$dev_detrend, method = "pearson"),
    "rho_ts" = cor(mod_data$.dissimilarity, mod_data$ts, method = "pearson"),
    "rho_ts_detrend" = cor(mod_data$.dissimilarity, mod_data$ts_detrend, method = "pearson")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  out <- rbind(
    marginaleffects::avg_slopes(m_dev, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list(
    "draws" = out,
    "rho" = correlation,
    "prior_sens" = prior_sens
  ))
})

dee_dissim_b <- lapply(1:1000, function(i) {
  mod_data <- subset(dee_dissim_blki, drawid == i)  |>
    dplyr::left_join(dee_stability_abun_blki, by = c("species", "colony"))
  
  m_dev <- brms::brm(
    brms::bf(dev ~ .dissimilarity, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept")
      ,
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_dev_detrend <- brms::brm(brms::bf(dev_detrend ~ .dissimilarity, family = "gaussian"),
                         data = mod_data,
                         prior = c(prior(normal(0, 1), class = "b"),
                                   prior(normal(0, 1), class = "Intercept"),
                                   prior(exponential(10), class = "sigma")),
                         iter = 2000, warmup = 1000, chains = 4, cores = 4,
                         control = list(adapt_delta = 0.99),
                         backend = "cmdstanr",
                         silent = 2)

  m_ts <- brms::brm(
    brms::bf(ts ~ .dissimilarity, family = "gaussian"),
    data = mod_data,
    prior = c(
      prior(normal(0, 1), class = "b"),
      prior(normal(0, 1), class = "Intercept"),
      prior(exponential(10), class = "sigma")
    ),
    iter = 2000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.99),
    backend = "cmdstanr",
    silent = 2
  )
  
  m_ts_detrend <- brms::brm(brms::bf(ts_detrend ~ .dissimilarity, family = "gaussian"),
                        data = mod_data,
                        prior = c(prior(normal(0, 1), class = "b"),
                                  prior(normal(0, 1), class = "Intercept"),
                                  prior(exponential(10), class = "sigma")),
                        iter = 2000, warmup = 1000, chains = 4, cores = 4,
                        control = list(adapt_delta = 0.99),
                        backend = "cmdstanr",
                        silent = 2)
  
  prior_sens <-  rbind(
    priorsense::powerscale_sensitivity(m_dev) |>
      dplyr::mutate(response = "dev"),
    priorsense::powerscale_sensitivity(m_dev_detrend) |>
      dplyr::mutate(response = "dev_detrend"),
    priorsense::powerscale_sensitivity(m_ts) |>
      dplyr::mutate(response = "ts"),
    priorsense::powerscale_sensitivity(m_ts_detrend) |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  correlation <- data.frame(
    "rho_dev" = cor(mod_data$.dissimilarity, mod_data$dev, method = "pearson"),
     "rho_dev_detrend" = cor(mod_data$.dissimilarity,mod_data$dev_detrend, method = "pearson"),
    "rho_ts" = cor(mod_data$.dissimilarity, mod_data$ts, method = "pearson"),
    "rho_ts_detrend" = cor(mod_data$.dissimilarity, mod_data$ts_detrend, method = "pearson")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  out <- rbind(
    marginaleffects::avg_slopes(m_dev, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "dev"),
    marginaleffects::avg_slopes(m_dev_detrend,variables = ".dissimilarity") |>
      marginaleffects::get_draws()|>
      dplyr::mutate(response = "dev_detrend"),
    marginaleffects::avg_slopes(m_ts, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts"),
    marginaleffects::avg_slopes(m_ts_detrend, variables = ".dissimilarity") |>
      marginaleffects::get_draws() |>
      dplyr::mutate(response = "ts_detrend")
  ) |>
    dplyr::mutate(mc_sim = i)
  
  return(list(
    "draws" = out,
    "rho" = correlation,
    "prior_sens" = prior_sens
  ))
})

saveRDS(dee_dissim_a,"outputs/stability/dee_dissim_atpu.RDS")
saveRDS(dee_dissim_b,"outputs/stability/dee_dissim_blki.RDS")

