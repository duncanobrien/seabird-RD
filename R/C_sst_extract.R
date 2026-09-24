##########################################################################################
#### Preamble ####
##########################################################################################
#https://cds.climate.copernicus.eu/datasets/derived-era5-single-levels-daily-statistics?tab=overview
#https://doi.org/10.24381/cds.4991cf48
library(terra)
library(foreach)
library(dplyr)

foraging_polys <- read.csv("data/Behaviour/OutputDat_Kernel.csv") |>
  dplyr::group_by(species, colony, season) |>
  dplyr::group_split()

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

sea_pol <- terra::vect("/Users/ul20791/Desktop/Academia/RA/Repositories/multidim-wp1/Data/SST/GOaS_v1_20211214/goas_v01.shp") #load global ocean shapefiles so land can be masked out of foraging kernels
sea_pol <- terra::subset(
  sea_pol,
  sea_pol$name %in% c("Arctic Ocean", "North Atlantic Ocean")
)

##########################################################################################
#### Load SST (Kelvin; var = sea_surface_temperature) and crop ####
##########################################################################################

sst_files <- fs::dir_ls(
  path =  "data/SST",
  glob = "*.nc",
  recurse = TRUE
) #list all .nc data files in the rotated folder

loop_grid <- expand.grid("raster" = seq_along(sst_files),
                         "poly" = seq_along(foraging_polys)) #define combinations of rasters and range polygons

merge_polys <- TRUE #boolean to indicate whether separate polygons should be pooled or not

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

sst <- foreach::foreach(
  i = 1:NROW(loop_grid),
  .combine = "rbind",
  .packages = c("terra", "dplyr", "tidyr"),
  .export = c("merge_polys")
) %dopar% {
  rast <- terra::rast(sst_files[loop_grid[i, 1]]) #load rotated raster

  poly <- foraging_polys[[loop_grid[i, 2]]]

  p <- cbind(id = as.numeric(factor(poly$group)),
             x = poly$long,
             y = poly$lat)
  pp <- terra::vect(p, "polygons", crs = my_wkt) #create polygon from kernel coordinates and project into longlat coordinates
  
  prj <- terra::crs(rast) #extract raster coordinate system
  pol <- terra::project(pp, prj) #project polygon(s) to a shared crs
  
  spol <- terra::project(sea_pol, prj) #project polygon(s) to a shared crs
  
  masked_pol <- terra::intersect(spol, pol)

  if (isFALSE(merge_polys)) {
    terra::values(rast) <- terra::values(rast) - 273.15 #convert from Kelvin to Celsius
    sst_mean <- terra::extract(rast, pol, fun = mean) |>
      tidyr::pivot_longer(-ID, names_to = NULL, values_to = "mean_sst") #estimate mean/sd/cv for the polygon(s)
    sst_sd <- terra::extract(rast, pol, fun = sd) |>
      tidyr::pivot_longer(-ID, names_to = NULL, values_to = "sd_sst")
    sst_cv <- terra::extract(
      rast,
      pol,
      fun = function(x) {
        sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
      }
    ) |>
      tidyr::pivot_longer(-ID, names_to = NULL, values_to = "cv_sst")

    daily_sst <- data.frame(
      "polygon" = sst_mean$ID,
      "date" = terra::time(rast),
      "mean_sst" = sst_mean$mean_sst,
      "sd_sst" = sst_sd$sd_sst,
      "cv_sst" = sst_cv$cv_sst
    )
  } else if (isTRUE(merge_polys)) {
    rcrop <- terra::crop(rast, masked_pol, mask = FALSE) #mask the raster to the polygon(s)
    rmask <- terra::mask(rcrop, masked_pol, touches = FALSE)
    poly_ts <-  terra::values(rmask) - 273.15 #convert from Kelvin to Celsius
    
    mean_poly <- apply(poly_ts, 2, mean, na.rm = TRUE) #estimate mean/sd/cv across the merged polygon(s)
    sd_poly <-  apply(poly_ts, 2, sd, na.rm = TRUE)
    cv_poly <-  apply(poly_ts, 2, function(x) {
      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
    })

    daily_sst <- data.frame(
      "date" = terra::time(rmask),
      "mean_sst" = mean_poly,
      "sd_sst" = sd_poly,
      "cv_sst" = cv_poly
    )
  }

  rownames(daily_sst) <- NULL

  return(
    daily_sst |>
      dplyr::mutate(
        date = as.Date(date),
        species = unique(poly$species),
        colony = unique(poly$colony),
        kernel = unique(poly$season)
      ) |> #convert character date to Date class
      dplyr::arrange(date)
  )
} |>
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

parallel::stopCluster(cl)

##########################################################################################
#### Summarise SST to year ####
##########################################################################################

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
    #arrival_date_mean = dplyr::if_else(colony == "Isfjorden", as.Date("2024-05-01"),arrival_date_mean),
    arrival_month_mean = lubridate::month(arrival_date_mean),
    arrival_day_mean = lubridate::day(arrival_date_mean),
    departure_date_mean = as.Date(departure_julian_mean, origin = "2024-12-31"),
    departure_month_mean = lubridate::month(departure_date_mean),
    departure_day_mean = lubridate::day(departure_date_mean)
  )

sst_daily <- sst |>
  dplyr::mutate(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  ) |>
  dplyr::left_join(
    subset(
      arrival_dates,
      select = c(species, colony, arrival_day_mean, arrival_month_mean, departure_day_mean, departure_month_mean)
    ),
    by = c("species","colony")) |>
  dplyr::mutate(season = ifelse( #assign breeding season based upon days and months
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  ),
  .by = c("species","colony"))

sst_yr <- sst_daily |>
  dplyr::filter(kernel == season) |>
  dplyr::mutate(year_season = dplyr::if_else(((month > departure_month_mean) |
                                                (month == departure_month_mean & day > departure_day_mean)),
                                             year + 1,
                                             year),
                .by = c(species, colony, kernel)) |> #as the non-breeding season spans January, if post breeding season, need to associate it with the following year
  dplyr::reframe(
    mean_sst = mean(mean_sst, na.rm = TRUE),
    sd_sst = mean(sd_sst, na.rm = TRUE),
    cv_sst = mean(cv_sst, na.rm = TRUE),
    .by = c(species, colony, kernel, year_season)
  ) |>
  dplyr::rename(season = kernel, year = year_season)

saveRDS(sst_daily,"data/SST/OutputDat_mean_SST_daily.RDS")
write.csv(sst_yr,"data/SST/OutputDat_mean_SST.csv", row.names = FALSE)
