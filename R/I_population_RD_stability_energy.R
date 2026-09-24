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

source("R/G_clean_data.R")
source("R/utilities/respDIVme.R")

##########################################################################################
#### Prepare data ####
##########################################################################################

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

bform_energy <- brms::bf(DEEkJ  | mi(DEEkJ_sd) ~ s(sst, bs="tp", k = 10) + t2(sst, bird1, bs= c("tp","re")) + (1|day1) + (1|year1))

##########################################################################################
#### Puffin models ####
##########################################################################################

RD_atpu <- RD_data_energy |>
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

priors_energy_atpu <- c(
  prior(normal(0, 3), class = "b"),
  prior(normal(log(800), 2), class ="Intercept"),
  prior(exponential(10), class = "shape"),
  prior(exponential(1), class = "sds")
  ,prior(exponential(1), class = "sd")
)

lapply(seq_along(RD_atpu)[4:5],function(x){
  
  inits  <- list(Yl = RD_atpu[[x]] |>
                   dplyr::pull(DEEkJ))
  
  m1 <- brms::brm(bform_energy,
                  family = Gamma(link = "log"),
                  data = RD_atpu[[x]] |>
                    dplyr::mutate(sst = (sst - unique(mean_sst))/unique(sd_sst)) |>
                    dplyr::mutate(year1 = factor(year1)),
                  chains = 4,
                  cores = 4,
                  prior = priors_energy_atpu,
                  iter = 2000,
                  warmup = 1000,
                  control = list(adapt_delta=0.999, max_treedepth = 15),
                  save_pars = save_pars(latent = FALSE),
                  init = rep(list(inits),4),
                  backend = "cmdstanr")
  
  newdata <- expand.grid("sst" = seq(min(m1$data[["sst"]]), max(m1$data[["sst"]]), length.out = 100),
                         "bird1" = as.character(unique(m1$data[["bird1"]]))
                         #,"year1" = unique(m1$data[["year1"]])
  ) |>
    dplyr::mutate(DEEkJ_sd = 150)
  
  rd_atpu_energy <- respDivme(m1, newdata, gradient = "sst", unit = "bird1", resp = "DEEkJ", error = "DEEkJ_sd", ndraws = 1000, type = "link", raw = TRUE)
  
  epreds_atpu_energy <- rd_atpu_energy[["epreds"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  derivs_atpu_energy <- rd_atpu_energy[["derivs"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  dissim_atpu_energy <- rd_atpu_energy[["dissim"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  diverg_atpu_energy <- rd_atpu_energy[["diverg"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  imbalance_atpu_energy <- rd_atpu_energy[["imbalance"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  
  sst_info <- data.frame("mean_sst" = unique(RD_atpu[[x]]$mean_sst),
                         "sd_sst" = unique(RD_atpu[[x]]$sd_sst)) |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))

  saveRDS(list("model" = m1,
               "sst_info" = sst_info,
               "epreds" = epreds_atpu_energy,
               "derivs" = derivs_atpu_energy,
               "dissim" = dissim_atpu_energy,
               "diverg" = diverg_atpu_energy,
               "imbalance" = imbalance_atpu_energy),paste0("outputs/models/atpu_energy_",x,".RDS"))
})

##########################################################################################
#### Kittiwake models ####
##########################################################################################

RD_blki <- RD_data_energy |>
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

priors_energy_blki <- c(
  prior(normal(0, 3), class = "b"),
  prior(normal(log(700), 2), class = "Intercept"),
  prior(exponential(100), class = "shape"),
  prior(exponential(1), class = "sds")
  ,prior(exponential(1), class = "sd")
)

lapply(seq_along(RD_blki)[1],function(x){
  
  inits  <- list(Yl = RD_blki[[x]] |>
                   dplyr::pull(DEEkJ))
  
  m1 <- brms::brm(bform_energy,
                  family = Gamma(link = "log"),
                  data = RD_blki[[x]] |>
                    dplyr::mutate(sst = (sst - unique(mean_sst))/unique(sd_sst)) |>
                    dplyr::mutate(year1 = factor(year1)),
                  chains = 4,
                  cores = 4,
                  prior = priors_energy_blki,
                  iter = 2000,
                  warmup = 1000,
                  control = list(adapt_delta=0.999, max_treedepth = 17),
                  #sample_prior = "only",
                  init = rep(list(inits),4),
                  save_pars = save_pars(latent = FALSE),
                  backend = "cmdstanr")
  
  newdata <- expand.grid("sst" = seq(min(m1$data[["sst"]]), max(m1$data[["sst"]]), length.out = 100),
                         "bird1" = as.character(unique(m1$data[["bird1"]]))
                         #,"year1" = unique(m1$data[["year1"]])
  ) |>
    dplyr::mutate(DEEkJ_sd = 200)
  
  rd_blki_energy <- respDivme(m1, newdata, gradient = "sst", unit = "bird1", resp = "DEEkJ", ndraws = 1000, type = "link", raw = TRUE)
  
  epreds_blki_energy <- rd_blki_energy[["epreds"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  derivs_blki_energy <- rd_blki_energy[["derivs"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  dissim_blki_energy <- rd_blki_energy[["dissim"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  diverg_blki_energy <- rd_blki_energy[["diverg"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  imbalance_blki_energy <- rd_blki_energy[["imbalance"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  
  sst_info <- data.frame("mean_sst" = unique(RD_blki[[x]]$mean_sst),
                         "sd_sst" = unique(RD_blki[[x]]$sd_sst)) |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  
  saveRDS(list("model" = m1,
               "sst_info" = sst_info,
               "epreds" = epreds_blki_energy,
               "derivs" = derivs_blki_energy,
               "dissim" = dissim_blki_energy,
               "diverg" = diverg_blki_energy,
               "imbalance" = imbalance_blki_energy),paste0("outputs/models/blki_energy_",x,".RDS"))
})
