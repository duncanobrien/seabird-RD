##########################################################################################
#### Preamble ####
##########################################################################################

library(data.table)
library(magrittr)
data.table::setDTthreads(threads = 1)

sst_daily <- readRDS("data/SST/OutputDat_mean_SST_daily.RDS") %>%
  as.data.table() %>%
  .[season == kernel,] %>%
  .[, season := NULL] %>%
  data.table::setnames("kernel", "season")

airtemp_daily <- readRDS("data/TAS/OutputDat_mean_airtemp_daily.RDS") %>%
  as.data.table() %>%
  subset(season_date == "breeding")

c_parameters <-  list("AtPu" = seq(1.5, 3, 0.1), "BlKi" = seq(1.5, 3, 0.1)) #coefficient for adjusting leg tucking
c1_parameters <- list("AtPu" = seq(106, 176, 1), "BlKi" = seq(2, 5.7, 0.1)) #cost of flight
c2_parameters <- list("AtPu" = 0, "BlKi" = seq(0.3, 3.8, 0.1)) #cost of foraging
c3_parameters <- list("AtPu" = seq(3.1, 15.3, 0.1), "BlKi" = seq(0.05, 1.1, 0.1)) #cost of being on land
c4_parameters <- list("AtPu" = seq(3.1, 15.3, 0.1), "BlKi" = seq(0.1, 2.8, 0.1)) # cost of resting on the water

##########################################################################################
#### Estimate puffin energy budgets ####
##########################################################################################

AtPu_energetics_env <- readRDS("data/Behaviour/OutputDat_Behav_Energy_AtPu.RDS") %>%
  .[, species := "AtPu"] %>%
  .[, colony := data.table::fcase(
    colony %in% "Hjelmsøya", 
    "Hjelmsoya",
    colony %in% "Hornøya",
    "Hornoya",
    colony %in% "Røst",
    "Rost",
    colony %in% "Runde and Ålesund",
    "Runde",
    colony %in% "Bjørnøya",
    "Bjornoya",
    default = colony
  )] %>%
  merge(subset(sst_daily, season == "breeding",
               select = c(date, mean_sst, species, colony, season)),
        all.x = TRUE,
        by = c("species", "colony", "date", "season")) %>%
  merge(
    subset(
      airtemp_daily,
      select = c(date, airtemp, species, colony, season_date),
      all.x = TRUE
    ) %>%
      data.table::setnames("season_date", "season"),
    by = c("species", "colony", "date", "season")
  ) %>%
  .[, `:=`(
    LCT_water = 14.18, #lower critical temperature in water
    LCT_air = 5.72, #lower critical temperature in air
    TC_water = 2.75, #thermal conductivity of water
    TC_air = 0.0282, #thermal conductivity of air
    convf = 1 / 803^0.689,
    #803 = costdivider
    weightG = 395,
    Beta_active = 118,  #intercept of active metabolic rate at 0°C
    Beta_rest = 72.2  #intercept of resting metabolic rate at 0°C
  )]

AtPu_iter <- lapply(1:1000,function(i){
  
  out <- copy(AtPu_energetics_env) %>%
    .[, `:=`(
      c = sample(c_parameters[["AtPu"]], 1),
      c1 = sample(c1_parameters[["AtPu"]], 1),
      c2 = sample(c2_parameters[["AtPu"]], 1),
      c3 = sample(c3_parameters[["AtPu"]], 1),
      c4 = sample(c4_parameters[["AtPu"]], 1)
    ), by = c("species", "colony", "individ_id")] %>%
    .[, `:=`(
      flightCoef_kj = c1 * 3.6,
      restCoef2_kj = c4 * 3.6,
      landCoef_kj = c3 * 3.6
    )] %>%
    .[, `:=`(
      flightConstant = (flightCoef_kj * convf) * weightG^0.689,
      activeConstant = (Beta_active * convf) * weightG^0.689,
      restConstant1 = (Beta_rest * convf) * weightG^0.689,
      restConstant2 = (restCoef2_kj * convf) * weightG^0.689,
      TC_water_Constant = (TC_water * convf) * weightG^0.689,
      TC_air_Constant = (TC_air * 20.1 / 1000) * 819.3 / (819.3^0.689) * weightG^0.689,
      landConstant = (landCoef_kj * convf) * weightG^0.689
    )] %>%
    .[, `:=`(
      restConstant_adjust = restConstant1 - LCT_water * TC_water_Constant,
      # So that rest is equal to land when sst > LCT
      activeConstant_adjust = activeConstant - LCT_water * TC_water_Constant,
      beta_land = landConstant + LCT_air * TC_air_Constant
    )] %>%
    .[, `:=`(
      DEEkJ_active = data.table::fifelse(
        mean_sst <= LCT_water,
        (activeConstant - mean_sst * TC_water_Constant) * tActive,
        activeConstant_adjust * tActive
      ),
      DEEkJ_rest = data.table::fifelse(
        mean_sst <= LCT_water,
        (restConstant1 - TC_water_Constant * mean_sst) * tRestWater,
        restConstant_adjust * tRestWater
      ),
      DEEkJ_flight = flightConstant * tFlight,
      DEEkJ_forage = 0,
      DEEkJ_restland = data.table::fifelse(
        airtemp <= LCT_air,
        (beta_land - TC_air_Constant * airtemp) * tLand,
        landConstant * tLand
      )    
    )] %>%
    .[,  DEEkJ := DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland] %>%
    subset(select = c(species, colony, individ_id, date, mean_sst, airtemp, DEEkJ))
  
  return(out)
  
}) %>%
  data.table::rbindlist() %>%
  .[,.(
    DEEkJ = mean(DEEkJ),
    DEEkJ_sd = sd(DEEkJ)
  ), by = c("species", "colony", "individ_id", "date", "mean_sst", "airtemp")]

saveRDS(AtPu_iter, "data/Energetics/OutputDat_Energetics_AtPu.RDS")

##########################################################################################
#### Estimate kittiwake energy budgets ####
##########################################################################################

BlKi_energetics_env <- readRDS("data/Behaviour/OutputDat_Behav_Energy_BlKi.RDS") %>%
  .[, species := "BlKi"] %>%
  .[, colony := data.table::fcase(
    colony %in% "Hjelmsøya", 
    "Hjelmsoya",
    colony %in% "Hornøya",
    "Hornoya",
    colony %in% "Røst",
    "Rost",
    colony %in% "Runde and Ålesund",
    "Runde",
    colony %in% "Bjørnøya",
    "Bjornoya",
    default = colony
  )] %>%
  merge(subset(sst_daily, season == "breeding",
               select = c(date, mean_sst, species, colony, season)),
        by = c("species", "colony", "date", "season")) %>%
  merge(
    subset(
      airtemp_daily,
      select = c(date, airtemp, species, colony, season_date)
    ) %>%
      data.table::setnames("season_date", "season"),
    by = c("species", "colony", "date", "season")
  )

BlKi_iter <- lapply(1:1000,function(i){
  c1 <- sample(c1_parameters[["BlKi"]], 1)
  c2 <- mean(sample(c2_parameters[["BlKi"]], 1), sample(seq(2, 5.7, 0.1), 1))
  c3 <- sample(c3_parameters[["BlKi"]], 1)
  c4 <- sample(c4_parameters[["BlKi"]], 1)
  
  out <- copy(BlKi_energetics_env) %>%
    .[, `:=`(
      LCT_water = 12.5, #lower critical temperature in water
      LCT_air = 4.5, #lower critical temperature in air
      TC_water = 0.0997, #thermal conductivity of water
      TC_air = 0.0466, #thermal conductivity of air
      c1 = c1,
      c2 = c2,
      c3 = c3,
      c4 = c4,
      restCoef = c4 / 24,
      forageCoef = c2 / 24,
      landCoef = c3 / 24,
      flightCoef = c1 / 24,
      weightG = 392,
      cf = 20.1,
      Beta_active = 0, #intercept of active metabolic rate at 0°C
      Beta_rest = 1.87 #intercept of resting metabolic rate at 0°C
    )] %>%
    .[, `:=`(
      flightConstantx = ((flightCoef * 450) / 450^0.717) * weightG^0.717,
      restConstant2x = ((restCoef * 450) / 450^0.717) * weightG^0.717,
      forageConstantx = ((forageCoef * 450) / 450^0.717) * weightG^0.717,
      landConstantx = ((landCoef * 450) / 450^0.717) * weightG^0.717,
      betax = (((Beta_rest * cf / 1000) * 365) / 365^0.717) * weightG^0.717,
      TCx_water = (((TC_water * cf / 1000) * 365) / 365^0.717) * weightG^0.717,
      TCx_air = (((TC_air * cf / 1000) * 365) / 365^0.717) * weightG^0.717
    )] %>%
    .[, `:=`(
      restConstant1x = (LCT_water * TCx_water + restConstant2x),
      beta_land = landConstantx + LCT_air * TCx_air
    )] %>%
    .[, `:=`(
      DEEkJ_active = 0,
      DEEkJ_rest = data.table::fifelse(
        mean_sst <= LCT_water,
        (restConstant1x - TCx_water * mean_sst) * tRestWater,
        restConstant2x * tRestWater
      ),
      DEEkJ_flight = flightConstantx * tFlight,
      DEEkJ_forage = forageConstantx * tForage,
      DEEkJ_restland = data.table::fifelse(
        airtemp <= LCT_air,
        (beta_land - airtemp * TCx_air) * tLand,
        landConstantx * tLand
      )
    )] %>%
    .[,  DEEkJ := DEEkJ_active + DEEkJ_rest + DEEkJ_flight + DEEkJ_restland] %>%
    subset(select = c(species, colony, individ_id, date, mean_sst, airtemp, DEEkJ))
  
  
  return(out)
  
}) %>%
  data.table::rbindlist() %>%
  .[,.(
    DEEkJ = mean(DEEkJ),
    DEEkJ_sd = sd(DEEkJ)
  ), by = c("species", "colony", "individ_id", "date", "mean_sst", "airtemp")]

saveRDS(BlKi_iter,"data/Energetics/OutputDat_Energetics_BlKi.RDS")
