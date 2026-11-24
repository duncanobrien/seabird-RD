##########################################################################################
#### Preamble ####
##########################################################################################

library(dplyr)
library(brms)
library(ordbetareg)
library(cmdstanr)
library(ggplot2)
library(tidybayes)
library(ggdist)
library(patchwork)

set.seed(333)

source("R/H_clean_data.R")
source("R/utilities/respDIV.R")

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

priors_behav <- c(
  prior(exponential(1), class = "sds")
  ,prior(exponential(1), class = "sd")
)

bform_behav <- bf(tForage_std ~ s(sst, bs="tp", k = 5) + t2(sst, bird1, bs= c("tp","re"), k = 5)  + (1 | day1) + (1|year1))

##########################################################################################
#### Puffin models ####
##########################################################################################

RD_atpu <- RD_data_behav |>
  subset(species == "AtPu") |>
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, tForage_std) |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

lapply(seq_along(RD_atpu),function(x){
  
  m1 <- ordbetareg::ordbetareg(bform_behav,
                               data = RD_atpu[[x]] |>
                                 dplyr::mutate(sst = (sst - unique(mean_sst))/unique(sd_sst)) |>
                                 dplyr::mutate(year1 = factor(year1)),
                               chains = 4,
                               cores = 4,
                               extra_prior = priors_behav,
                               iter = 2000,
                               warmup = 1000,
                               control = list(adapt_delta=0.99, max_treedepth = 12),
                               backend = "cmdstanr")
  
  newdata <- expand.grid("sst" = seq(min(m1$data[["sst"]]), max(m1$data[["sst"]]), length.out = 100),
                         "bird1" = as.character(unique(m1$data[["bird1"]]))
                         #,"year1" = unique(m1$data[["year1"]])
                         )

  rd_atpu_behav <- respDiv(m1, newdata, gradient = "sst", unit = "bird1", resp = "tForage_std", ndraws = 1000, type = "link", raw = TRUE)
  
  epreds_atpu_behav <- rd_atpu_behav[["epreds"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  derivs_atpu_behav <- rd_atpu_behav[["derivs"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  dissim_atpu_behav <- rd_atpu_behav[["dissim"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  diverg_atpu_behav <- rd_atpu_behav[["diverg"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  imbalance_atpu_behav <- rd_atpu_behav[["imbalance"]] |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  
  sst_info <- data.frame("mean_sst" = unique(RD_atpu[[x]]$mean_sst),
                         "sd_sst" = unique(RD_atpu[[x]]$sd_sst)) |>
    dplyr::mutate(species = unique(RD_atpu[[x]]$species),
                  colony = unique(RD_atpu[[x]]$colony))
  
  saveRDS(list("model" = m1,
               "sst_info" = sst_info,
               "epreds" = epreds_atpu_behav,
               "derivs" = derivs_atpu_behav,
               "dissim" = dissim_atpu_behav,
               "diverg" = diverg_atpu_behav,
               "imbalance" = imbalance_atpu_behav),paste0("outputs/models/atpu_behav_",x,".RDS"))
})

##########################################################################################
#### Kittiwake models ####
##########################################################################################

RD_blki <- RD_data_behav |>
  subset(species == "BlKi") |>
  subset(!colony %in% c("Isfjorden", "Hjelmsoya")) |>
  subset(!bird1 %in% c("BlKi_Anda_BlKi_62","BlKi_Anda_BlKi_370", "BlKi_Hornoya_BlKi_90", "BlKi_Sklinna_BlKi_556")) |> #drop birds that have unexplained periods of zero foraging across all conditions
  dplyr::mutate(mean_sst = mean(sst),
                sd_sst = sd(sst),.by = colony) |>
  dplyr::filter(dplyr::row_number()==1,.by = c(date, bird1)) |>
  dplyr::filter(length(unique(year)) >= 3,.by = bird1) |>
  dplyr::filter(dplyr::n_distinct(bird1)>1, .by = colony) |>
  dplyr::select(species, colony, year, date, sst, mean_sst, sd_sst, bird1, year1, day1, tForage_std) |>
  dplyr::group_by(species,colony) |>
  dplyr::group_split() 

lapply(seq_along(RD_blki),function(x){
  m1 <- ordbetareg::ordbetareg(bform_behav,
                               data = RD_blki[[x]] |>
                                 dplyr::mutate(sst = (sst - unique(mean_sst))/unique(sd_sst)) |>
                                 dplyr::mutate(year1 = factor(year1)),
                               chains = 4,
                               cores = 4,
                               extra_prior = priors_behav,
                               iter = 2000,
                               warmup = 1000,
                               control = list(adapt_delta=0.99, max_treedepth = 15),
                               backend = "cmdstanr")
  
  newdata <- expand.grid("sst" = seq(min(m1$data[["sst"]]), max(m1$data[["sst"]]), length.out = 100),
                         "bird1" = as.character(unique(m1$data[["bird1"]]))
                         #,"year1" = unique(m1$data[["year1"]])
  )
  
  rd_blki_behav <- respDiv(m1, newdata, gradient = "sst", unit = "bird1", resp = "tForage_std", ndraws = 1000, type = "link", raw = TRUE)
  
  epreds_blki_behav <- rd_blki_behav[["epreds"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  derivs_blki_behav <- rd_blki_behav[["derivs"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  dissim_blki_behav <- rd_blki_behav[["dissim"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  diverg_blki_behav <- rd_blki_behav[["diverg"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  imbalance_blki_behav <- rd_blki_behav[["imbalance"]] |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  
  sst_info <- data.frame("mean_sst" = unique(RD_blki[[x]]$mean_sst),
                         "sd_sst" = unique(RD_blki[[x]]$sd_sst)) |>
    dplyr::mutate(species = unique(RD_blki[[x]]$species),
                  colony = unique(RD_blki[[x]]$colony))
  
  saveRDS(list("model" = m1,
               "sst_info" = sst_info,
               "epreds" = epreds_blki_behav,
               "derivs" = derivs_blki_behav,
               "dissim" = dissim_blki_behav,
               "diverg" = diverg_blki_behav,
               "imbalance" = imbalance_blki_behav),paste0("outputs/models/blki_behav_",x,".RDS"))
  
})
