##########################################################################################
#### Preamble ####
##########################################################################################

library(dplyr)
library(arrow)
library(ggplot2)
library(sp)
library(adehabitatHR)

my_wkt = 'PROJCS["Lambert_Azimuthal_Equal_Area",
 GEOGCS["GCS_WGS_1984",
  DATUM["D_WGS_1984",
   SPHEROID["WGS_1984",6378137.0,298.257223563]],
  PRIMEM["Greenwich",0.0],
  UNIT["Degree",0.0174532925199433]],
 PROJECTION["Lambert_Azimuthal_Equal_Area"],
 PARAMETER["False_Easting",0.0],
 PARAMETER["False_Northing",0.0],
 PARAMETER["Central_Meridian",22.0],
 PARAMETER["Latitude_Of_Origin",70.0],
 UNIT["Meter",1.0]]'

arrival_dates <- readRDS("data/Behaviour/BreedingSeasonStartEndDates.rds") |>
  subset(select = c(Species, Colony, ArrMean, DepMean)) |>
  setNames(c(
    "species",
    "colony",
    "arrival_julian_mean",
    "departure_julian_mean"
  )) |>
  dplyr::mutate(species = ifelse(species %in% "Fratercula_arctica", "AtPu", "BlKi")) |>
  dplyr::mutate(
    colony = dplyr::case_when(
      colony %in% "Hjelmsoya" ~ "Hjelmsøya",
      colony %in% "Hornoya" ~ "Hornøya",
      colony %in% "Rost" ~ "Røst",
      colony %in% "Runde and Alesund" ~ "Runde and Ålesund",
      colony %in% "Bjornoya" ~ "Bjørnøya",
      TRUE ~ colony
    )
  ) |> #need to standardise colony names with other files
  dplyr::mutate(arrival_date_mean = as.Date(arrival_julian_mean, origin = "2024-12-31")) |> #assume mean day is for a non-leapyear (2025 isn't a leap year) and that Jan 1st is day 1 (2024-12-31 is day 0)
  dplyr::mutate(arrival_date_mean = dplyr::if_else(
    colony == "Isfjorden",
    as.Date("2024-05-01"),
    arrival_date_mean
  )) |>
  dplyr::mutate(
    arrival_month_mean = lubridate::month(arrival_date_mean),
    arrival_day_mean = lubridate::day(arrival_date_mean)
  ) |>
  dplyr::mutate(departure_date_mean = as.Date(departure_julian_mean, origin = "2024-12-31")) |>
  dplyr::mutate(
    departure_month_mean = lubridate::month(departure_date_mean),
    departure_day_mean = lubridate::day(departure_date_mean)
  )

AtPu_locs <- arrow::read_parquet("data/Behaviour/Layton_MULTISTRESS_IRMA_positional_data_2026-05-22.gz.parquet") |>
  subset(species == "Atlantic puffin") |>
  subset(
    colony %in% c(
      "Anda",
      "Hjelmsøya",
      "Hornøya",
      "Røst",
      "Runde and Ålesund"
    )
  ) |>
  dplyr::mutate(date = as.Date(date_time),
                time = format(date_time, format = "%H:%M:%S")) |>
  dplyr::mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  ) |>
  dplyr::left_join(subset(
    arrival_dates,
    species %in% "AtPu",
    select = c(
      colony,
      arrival_day_mean,
      arrival_month_mean,
      departure_day_mean,
      departure_month_mean
    )
  ),
  by = c("colony")) |>
  dplyr::mutate(season = ifelse(
    #assign breeding season based upon days and months
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  )) |>
  dplyr::mutate(
    species = "AtPu",
    colony_season = paste(colony, season, sep = "_"),
    colony_season_year = paste(colony, season, year, sep = "_")
  ) |>
  subset(year < 2024)

BlKi_locs <- arrow::read_parquet("data/Behaviour/Layton_MULTISTRESS_IRMA_positional_data_2026-05-22.gz.parquet") |>
  subset(species == "Black-legged kittiwake") |>
  subset(
    colony %in% c(
      "Anda",
      "Bjørnøya",
      "Hjelmsøya",
      "Hornøya",
      "Isfjorden",
      "Røst",
      "Sklinna",
      "Runde and Ålesund"
    )
  ) |>
  dplyr::mutate(date = as.Date(date_time),
                time = format(date_time, format = "%H:%M:%S")) |>
  dplyr::mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  ) |>
  dplyr::left_join(subset(
    arrival_dates,
    species %in% "AtPu",
    select = c(
      colony,
      arrival_day_mean,
      arrival_month_mean,
      departure_day_mean,
      departure_month_mean
    )
  ),
  by = c("colony")) |>
  dplyr::mutate(season = ifelse(
    #assign breeding season based upon days and months
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  )) |>
  dplyr::mutate(
    species = "BlKi",
    colony_season = paste(colony, season, sep = "_"),
    colony_season_year = paste(colony, season, year, sep = "_")
  ) |>
  subset(year < 2024)

##########################################################################################
#### Estimate 50% distribution kernels ####
##########################################################################################

kernels <- do.call("rbind", lapply(list(AtPu_locs,BlKi_locs), function(x) {

  coords <- sp::SpatialPointsDataFrame(coords = x |>
                                         dplyr::select(lon, lat),
                                       data = x |>
                                         dplyr::select(colony_season),
                                       proj4string = sp::CRS('+proj=longlat +datum=WGS84 +no_defs +type=crs')) |> # define projection of original dataset
    sp::spTransform(coords, CRSobj = sp::CRS(my_wkt)) # reproject

  kerneldata <- adehabitatHR::kernelUD(coords, h = 'href', grid = 1000, extent = 1.0) #apply utilization distribution kernel estimator

  kd <- adehabitatHR::getverticeshr(kerneldata, 50) |> #50% quantiles following defintion of 'core area'
    ggplot2::fortify() |>
    dplyr::mutate(
      species = unique(x$species)) |>
    tidyr::separate_wider_delim(id, delim = "_", names = c("colony","season"))

  return(kd)

}))

##########################################################################################
#### Save out ####
##########################################################################################

write.csv(kernels,"data/Behaviour/OutputDat_Kernel_95.csv", row.names = FALSE)
