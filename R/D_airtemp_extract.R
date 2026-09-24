##########################################################################################
#### Preamble ####
##########################################################################################
#https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=overview
#https://doi.org/10.24381/cds.4991cf48
library(terra)
library(foreach)
library(dplyr)
library(lubridate)

arrival_dates <- readRDS("data/Behaviour/BreedingSeasonStartEndDates.rds") |>
  subset(select = c(Species, Colony, ArrMean, DepMean)) |>
  setNames(c(
    "species",
    "colony",
    "arrival_julian_mean",
    "departure_julian_mean"
  )) |>
  dplyr::mutate(species = dplyr::if_else(species %in% "Fratercula_arctica", "AtPu", "BlKi")) |> #rename species and colony for consistency
  dplyr::mutate(
    colony = dplyr::case_when(
      colony == "Hjelmsøya" ~ "Hjelmsoya",
      colony == "Hornøya" ~ "Hornoya",
      colony == "Røst" ~ "Rost",
      colony == "Runde and Alesund" ~ "Runde",
      colony == "Bjørnøya" ~ "Bjornoya",
      TRUE ~ colony
    )
  )  |> #clean colony names to be consistent across dataframes
  dplyr::mutate(
    arrival_date_mean = as.Date(arrival_julian_mean, origin = "2024-12-31"),
    #assume mean day is for a non-leapyear (2025 isn't a leap year) and that Jan 1st is day 1 (2024-12-31 is day 0)
    arrival_month_mean = lubridate::month(arrival_date_mean),
    arrival_day_mean = lubridate::day(arrival_date_mean),
    departure_date_mean = as.Date(departure_julian_mean, origin = "2024-12-31"),
    departure_month_mean = lubridate::month(departure_date_mean),
    departure_day_mean = lubridate::day(departure_date_mean)
  )

colony_locs <- rbind(
  readRDS("data/Behaviour/puffin_position_data.RDS") |>
    dplyr::select(species, colony, col_lon, col_lat),
  readRDS("data/Behaviour/kittiwake_position_data.RDS") |>
    dplyr::select(species, colony, col_lon, col_lat)
) |> #load bird location data which also contains colony coordinates
  dplyr::distinct() |> #drop duplicated rows
  dplyr::mutate(species = dplyr::if_else(species == "Atlantic puffin", "AtPu", "BlKi")) |>
  dplyr::mutate(
    colony = dplyr::case_when(
      colony == "Hjelmsøya" ~ "Hjelmsoya",
      colony == "Hornøya" ~ "Hornoya",
      colony == "Røst" ~ "Rost",
      colony == "Runde and Ålesund" ~ "Runde",
      colony == "Bjørnøya" ~ "Bjornoya",
      TRUE ~ colony
    )
  )  |> #clean colony names to be consistent across dataframes
  dplyr::rename(lon = col_lon, lat = col_lat) |>
  dplyr::mutate(ID = 1:dplyr::n()) |>
  dplyr::left_join(arrival_dates, by = c("species", "colony"))

##########################################################################################
#### Load surface temperature (Kelvin @2m; var = 2m_temperature) and crop ####
##########################################################################################

airtemp_files <- fs::dir_ls(
  path ="Data/TAS",
  glob = "*.nc",
  recurse = TRUE
) #list all .nc data files in the rotated folder

cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

airtemp <- foreach::foreach(
  i = seq_along(airtemp_files),
  .combine = "rbind",
  .packages = c("terra", "dplyr", "tidyr"),
  .export = c("colony_locs")
) %dopar% {

  rast <- terra::rast(airtemp_files[i]) #load rotated raster
  terra::values(rast) <- terra::values(rast) - 273.15 #convert from Kelvin to Celsius

  colony_vect <- terra::vect(colony_locs, geom = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84")

  prj <- terra::crs(rast) #extract raster coordinate system
  pol <- terra::project(colony_vect, prj)

  airtemp_daily <- terra::extract(x = rast, y = pol) |>
    tidyr::pivot_longer(-ID, names_to = NULL, values_to = "airtemp")

  daily_airtemp <- data.frame(
    "ID" = airtemp_daily$ID,
    "date" = terra::time(rast),
    "airtemp" = airtemp_daily$airtemp
  )

  rownames(daily_airtemp) <- NULL

  return(
    daily_airtemp |>
      dplyr::mutate(date = as.Date(date)) |> #convert character date to Date class
      dplyr::arrange(date) |>
      dplyr::left_join(colony_locs, by = "ID")
  )

}

parallel::stopCluster(cl)

##########################################################################################
#### Summarise surface temperature to year ####
##########################################################################################

airtemp_daily <- airtemp |>
  dplyr::mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  ) |>
  dplyr::mutate(season_date = ifelse( #assign breeding season based upon days and months
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  ),
  .by = c("species","colony"))

airtemp_yr <- airtemp_daily |>
  dplyr::mutate(year = lubridate::year(date),
                 month = lubridate::month(date),
                 day = lubridate::day(date)) |>
  dplyr::mutate(season = dplyr::if_else(
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  )) |>
  subset(season == "breeding") |>
  dplyr::reframe(
    mean_airtemp = mean(airtemp, na.rm = TRUE),
    .by = c(species, colony, year, season)
  ) |>
  dplyr::select(species, colony, season, year, mean_airtemp)

saveRDS(airtemp_daily,"data/TAS/OutputDat_mean_airtemp_daily.RDS")
write.csv(airtemp_yr,"data/TAS/OutputDat_mean_airtemp.csv", row.names = FALSE)
