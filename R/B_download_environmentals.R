##########################################################################################
#### Preamble ####
##########################################################################################

library(ecmwfr)

loop_grid <- expand.grid("month" = 3:9,
                         "year" = 1970:2024)

##########################################################################################
#### Download SST ####
##########################################################################################

sst_requests <- lapply(c(1:NROW(loop_grid)),function(i){
  yr = loop_grid$year[i]
  month = ifelse(nchar(loop_grid$month[i]) >1,paste0(loop_grid$month[i]),paste0("0",loop_grid$month[i]))  
  list(
    dataset_short_name = "derived-era5-single-levels-daily-statistics",
    product_type = "reanalysis",
    variable = "sea_surface_temperature",
    year = yr,
    month = month,
    day = c("01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31") ,
    time_zone = "utc+00:00",
    frequency = "1_hourly",
    daily_statistic = "daily_mean",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(90, -180, 0, 180),
    target = paste0("ERA5_SST_",yr,"_",month,".nc")
  )
})

sst_files <- ecmwfr::wf_request_batch(sst_requests,
                                 workers = 5,
                                 path  = "data/SST",       # store data in current working directory
                                 retry = 60,
                                 time_out =  3600*3
)

##########################################################################################
#### Download TAS ####
##########################################################################################

tas_requests <- lapply(1:NROW(loop_grid),function(i){
  yr = loop_grid$year[i]
  month = ifelse(nchar(loop_grid$month[i]) >1,paste0(loop_grid$month[i]),paste0("0",loop_grid$month[i]))  
  list(
    dataset_short_name = "derived-era5-single-levels-daily-statistics",
    product_type = "reanalysis",
    variable = "2m_temperature",
    year = yr,
    month = month,
    day = c("01", "02", "03",
            "04", "05", "06",
            "07", "08", "09",
            "10", "11", "12",
            "13", "14", "15",
            "16", "17", "18",
            "19", "20", "21",
            "22", "23", "24",
            "25", "26", "27",
            "28", "29", "30",
            "31") ,
    time_zone = "utc+00:00",
    frequency = "1_hourly",
    daily_statistic = "daily_mean",
    data_format = "netcdf",
    download_format = "unarchived",
    area = c(90, -180, 0, 180),
    target = paste0("ERA5_TAS_",yr,"_",month,".nc")
  )
})

tas_files <- ecmwfr::wf_request_batch(tas_requests,
                                      workers = 5,
                                      path  = "data/TAS",       # store data in current working directory
                                      retry = 60,
                                      time_out =  3600*3
)

