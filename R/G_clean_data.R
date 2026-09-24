##########################################################################################
#### Preamble ####
##########################################################################################

library(dplyr)
library(tidyr)
library(fs)
library(tibble)
library(janitor)
library(lubridate)
library(readxl)

##########################################################################################
#### Abundance ####
##########################################################################################

lookup_vec <- c(
  year = "krykkje",
  year = "ar",
  year = "a_r",
  est_count = "antall",
  est_count = "antall_kyst",
  site_number = "feltnr_lok",
  site_number = "dellokalitet_nr",
  colony = "keysite",
  colony = "superloc_keysite",
  species = "artsnavn_no_art"
)

abundance_isf <- readxl::read_excel("data/Abundances/Isfjorden_BlKi_count.xlsx") |>
  janitor::row_to_names(row_number = 1) |> #remove excess row
  janitor::clean_names() |>
  dplyr::rename(dplyr::any_of(lookup_vec)) |> #rename columns using the `lookup_vec`
  dplyr::mutate(year = as.numeric(year)) |>
  dplyr::mutate(colony = "Isfjorden", species =  "BlKi") |> #using file name, label colony and species
  dplyr::mutate(est_count = as.numeric(est_count), #convert from character to numeric
                site_number = 1) |> #add missing number of subsite - assumed singular subsite
  dplyr::reframe(
    tot_count = est_count,
    nsites = length(1:dplyr::n()),
    .by = c(species, colony, year)
  )

abundance_rost <- readxl::read_excel("data/Abundances/Krykkje overvåkingsdata Røst til Kate.xlsx") |>
  janitor::row_to_names(row_number = 2) |> #set colnames as plot names
  janitor::clean_names() |>
  dplyr::slice(-1) |> #drop irrelevant row
  dplyr::rename(year = na) |>
  dplyr::select(year, k_i, k_ii, k_iii, k_iv, k_v, k_xii, k_xiii, k_viii, k_viii_2, skrent_a, skrent_b, k_ix, k_x, greger, tyvsoy) |> #only select estimates at certain plots
  dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(.x)))|>
  dplyr::rowwise() |>
  dplyr::mutate(vedoy = sum(k_i, k_ii, k_iii, k_iv, k_v, k_xii, na.rm = TRUE),
                n_vedoy = length(!is.na(c(k_i, k_ii, k_iii, k_iv, k_v, k_xii))),
                harbour = sum(k_viii, k_viii_2, skrent_a, skrent_b, k_ix, k_x, greger, tyvsoy, na.rm = TRUE),
                n_harbour = length(!is.na(c(k_viii, k_viii_2, skrent_a, skrent_b, k_ix, k_x, greger, tyvsoy)))) |>
  dplyr::mutate(total = sum(vedoy,harbour),
                n_total = sum(n_vedoy,n_harbour)) |>
  dplyr::mutate(ref_sample = 17296/408) |> #true abundance divided by sampled abundance on Vedoy
  dplyr::mutate(vedoy_weight = (vedoy*(ref_sample*n_vedoy/(ref_sample*n_vedoy+n_harbour))),
                harbour_weight = (harbour*(n_harbour/(ref_sample*n_vedoy+n_harbour))),
                tot_count = vedoy_weight + harbour_weight,
                nsites = sum(n_vedoy,n_harbour)) |>
  dplyr::mutate(vedoy_contr = vedoy_weight/tot_count,
                harbour_contr = harbour_weight/tot_count) |>
  dplyr::ungroup() |>
  dplyr::mutate(tot_count = round(tot_count, digits = 0)) |> #convert weighted average to integer
  dplyr::mutate(colony = "Røst", species = "BlKi") |>
  subset(year >= 1988) #subset following Tycho's advice

abundances <- read.csv("data/Abundances/SEAPOP_Abundances_forMultiStress_noDupl.csv") |>
  janitor::clean_names() |>
  dplyr::rename(dplyr::any_of(lookup_vec)) |> #rename columns using the `lookup_vec`
  dplyr::mutate(
    species = dplyr::case_when(
      species == "Lunde" ~ "AtPu",
      species == "Krykkje" ~ "BlKi",
      TRUE ~ NA
    )
  ) |>
  dplyr::reframe(
    tot_count = sum(est_count, na.rm = TRUE),
    nsites = length(1:dplyr::n()),
    .by = c(species, colony, year)
  ) |> #(?) sum counts across subsites (estimate spatial area too)
  subset(interaction(species, colony) != "BlKi.Røst") |> #drop Rost kittiwakes as we're using a different dataset
  dplyr::bind_rows(abundance_isf, abundance_rost |>
                     dplyr::select(species,colony,year,tot_count,nsites)) |> #add colonies with separate data files
  dplyr::mutate(
    colony = dplyr::case_when(
      colony == "Hjelmsøya" ~ "Hjelmsoya",
      colony == "Hornøya" ~ "Hornoya",
      colony == "Røst" ~ "Rost",
      colony == "Runde and Ålesund" ~ "Runde",
      colony == "Bjørnøya" ~ "Bjornoya",
      TRUE ~ colony
    )
  )


##########################################################################################
#### Behaviour/Energetics ####
##########################################################################################

energetics_files <- fs::dir_ls(path = "Data/Energetics/", regexp = "OutputDat_Energetics.*(AtPu|BlKi)\\.RDS")

energetics <- do.call("rbind", lapply(
  energetics_files,
  FUN = function(file) {
    readRDS(file) |>
      dplyr::mutate(year = as.numeric(format(date, format = "%Y"))) |>
      dplyr::mutate(species = strsplit(file, "_")[[1]][3]) |>
      dplyr::mutate(species = gsub("\\.RDS", "", species))
  }))
  
#energetics_behav_files <- fs::dir_ls(path = "Data/Behaviour/", regexp = "OutputDat_Behav_Energy.*(AtPu|BlKi).*(light)\\.RDS")
energetics_behav_files <- fs::dir_ls(path = "data/Behaviour/", regexp = "OutputDat_Behav_Energy.*(AtPu|BlKi)\\.RDS")

energetics_behav_prop <- do.call("rbind", lapply(
  energetics_behav_files,
  FUN = function(file) {
    readRDS(file) |>
      dplyr::mutate(year = as.numeric(format(date, format = "%Y"))) |>
      dplyr::mutate(species = strsplit(file, "_")[[1]][4]) |>
      dplyr::mutate(species = gsub("\\.RDS", "", species))
  })) |>
  dplyr::mutate(colony = dplyr::case_when(
    colony %in% "Hjelmsøya" ~ "Hjelmsoya",
    colony %in% "Hornøya" ~ "Hornoya",
    colony %in% "Røst" ~ "Rost",
    colony %in% "Runde and Ålesund" ~ "Runde",
    colony %in% "Bjørnøya" ~ "Bjornoya",
    TRUE ~ colony)) |>
  dplyr::mutate(tForage_std = tForage/tDaylight) |>
  dplyr::filter(tForage_std <=1.0)
  
##########################################################################################
#### Environmentals ####
##########################################################################################

sst_daily <- readRDS("data/SST/OutputDat_mean_SST_daily.RDS") |>
  dplyr::filter(kernel == season) |>
  dplyr::select(-season) |>
  dplyr::rename(season = kernel)

sst_yr <- sst_daily |>
  dplyr::reframe(mean_sst = mean(mean_sst),.by = c(species,colony,season,year))
