##########################################################################################
#### Preamble ####
##########################################################################################

library(data.table)
data.table::setDTthreads(threads = 1)
library(magrittr)
library(lubridate)
library(suncalc)

arrival_dates <- readRDS("data/Behaviour/BreedingSeasonStartEndDates.rds") %>%
  data.table::as.data.table() %>%
  subset(select = c(Species, Colony, ArrMean, DepMean)) %>%
  data.table::setnames(., c("species", "colony", "arrival_julian_mean","departure_julian_mean")) %>%
  .[, species := data.table::fifelse(species %in% "Fratercula_arctica", "AtPu", "BlKi")] %>% #rename species and colony for consistency
  .[,colony :=   data.table::fcase(
    colony %in% "Hjelmsoya","Hjelmsøya",
    colony %in% "Hornoya","Hornøya",
    colony %in% "Rost","Røst",
    colony %in% "Runde and Alesund","Runde and Ålesund",
    colony %in% "Bjornoya","Bjørnøya",
    default = colony)
  ] %>% #need to standardise colony names with other files
  .[, arrival_date_mean := as.Date(arrival_julian_mean, origin = "2024-12-31")] %>% #assume mean day is for a non-leapyear (2025 isn't a leap year) and that Jan 1st is day 1 (2024-12-31 is day 0)
  .[, `:=`(
    arrival_month_mean = lubridate::month(arrival_date_mean),
    arrival_day_mean = lubridate::day(arrival_date_mean)
  )] %>%
  .[, departure_date_mean := as.Date(departure_julian_mean, origin = "2024-12-31")] %>% #assume mean day is for a non-leapyear (2025 isn't a leap year) and that Jan 1st is day 1 (2024-12-31 is day 0)
  .[, `:=`(
    departure_month_mean = lubridate::month(departure_date_mean),
    departure_day_mean = lubridate::day(departure_date_mean)
  )]

##########################################################################################
#### Process puffin data ####
##########################################################################################

AtPu_light <- readRDS("data/Behaviour/puffin_light_data.RDS") %>% #puffin data
  data.table::as.data.table()
AtPu_wet <- readRDS("data/Behaviour/puffin_immersion_data.RDS") %>% #puffin immersion data
  data.table::as.data.table() %>%
  merge(
    #4 individuals not present in AtPu_light
    subset(
      AtPu_light,
      select = c(
        colony,
        session_id,
        individ_id,
        date_time,
        year_tracked,
        std_light
      )
    ),
    by = c(
      "session_id",
      "date_time",
      "year_tracked",
      "individ_id",
      "colony"
    ),
    all.x = TRUE
  ) %>%
  merge(
    subset(
      arrival_dates,
      species %in% "AtPu",
      select = c(colony, arrival_day_mean, arrival_month_mean, departure_day_mean, departure_month_mean)
    ),
    by = c("colony"),
    all.x = TRUE
  ) %>%
  .[, `:=`(date = as.Date(date_time),
           time = format(date_time, format = "%H:%M:%S"))] %>%
  .[, `:=`(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  )] %>%
  .[, season := data.table::fifelse(
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  )] %>% #classify dates as breeding vs non-breeding and subset to breeding season
  subset(season == "breeding") 

AtPu_colonies <- readRDS("data/Behaviour/puffin_position_data.RDS") %>% #puffin colony coordinates
  subset(select = c(colony, col_lon, col_lat)) %>%
  unique()

AtPu_breeding_dates <- data.table::copy(AtPu_wet) %>%
  merge(
    AtPu_colonies,
    by = c("colony"),
    all.x = TRUE,
    sort = FALSE
  ) %>%
  subset(season %in% "breeding",
         select = c("colony", "date", "col_lat", "col_lon")) %>%
  unique() %>% #find unique combinations of colony, day and coordinates
  data.table::setnames(c("col_lat", "col_lon"), c("lat", "lon")) %>% #rename day to date to match `suncalc::getSunlightTimes` requirements
  .[, as.list(cbind(
    "colony" = colony,
    suncalc::getSunlightTimes(
      data = .,
      tz = "UTC",
      keep = c("sunrise", "sunset")
    )
  ))] #estimate sunrise and sunset times

AtPu_wet <- AtPu_wet %>%
  merge(
    subset(AtPu_breeding_dates, select = c(colony, date, sunrise, sunset)),
    by = c("colony", "date"),
    all.x = TRUE
  ) %>%
  .[, isdaylight :=
      data.table::fcase(
        #assign whether sampling time is during daylight hours
        (season %in% "breeding" &
           (
             is.na(sunrise) == TRUE |
               (date_time >= sunrise - 3600 &
                  date_time <= sunset + 3600)
           )),
        TRUE,
        (season %in% "nonbreeding" & std_light > 0),
        TRUE,
        default = FALSE
      )] %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, `:=`(
    isdry = std_conductivity <= 0.02,
    isdark = std_light <= 0
  )] %>% #assign whether dry and dark based upon our species specific threshold
  .[, behaviour :=
      data.table::fcase(
        ((!isdark | (is.na(isdark) & isdaylight)) & std_conductivity > 0.02 & std_conductivity < 0.98),
        "foraging",
        #"foraging" = daylight and immersion >2% dry and > 2% wet - works all year
        (std_conductivity >= 0.98),
        "restwater"
        #resting on sea surface = daylight and immersion >=98% wet - works all year
      )]

##########################################################################################
#### Puffin energetic classification ####
##########################################################################################

#estimate bouts of dry activity, intially assuming this is "flight"
AtPu_flightbouts <- data.table::copy(AtPu_wet) %>%
  subset(isdry) %>%
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  #find time difference between activity measurements
  .[, timediff := data.table::nafill(timediff, type = "const", fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>% #only keep the recordings at the start of a bout...
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

AtPu_flightbout_lengths <- data.table::copy(AtPu_wet) %>%
  merge(AtPu_flightbouts, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo := data.table::nafill(BoutNo, "locf"), by = c("colony", "individ_id", "session_id","isdry")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(isdry, length(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo")]
#calculate length of time (in minutes) spent in a dry bout

AtPu_darkness_boutNo <- data.table::copy(AtPu_flightbout_lengths) %>%
  subset(!isdaylight) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  #find time difference between activity measurements
  .[, timediff := data.table::nafill(timediff, type = "const", fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>% #only keep the recordings at the start of a bout...
  .[, BoutNo_dark := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo_dark)) 

AtPu_darkness_propDry <-  data.table::copy(AtPu_wet) %>%
  merge(AtPu_darkness_boutNo, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo_dark := data.table::nafill(BoutNo_dark, type = "locf"), by = c("colony", "individ_id", "session_id","isdaylight")] %>% #fill missing bouts with nearest BoutNo
  .[, `:=`(totDry = sum(1-std_conductivity)*10,
           darknessLength = data.table::uniqueN(date_time)*10), by =  c("colony", "individ_id", "session_id", "isdaylight", "BoutNo_dark")] %>%
  # .[, `:=`(totDry = sum(isdry)*10,
  #          darknessLength = data.table::uniqueN(date_time)*10), by =  c("colony", "individ_id", "session_id", "isdaylight", "BoutNo_dark")] %>%
  .[, propDarkDry := data.table::fifelse(!isdaylight, totDry/darknessLength, NA)] %>%
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo_dark, propDarkDry))

AtPu_activityAdjust1 <- data.table::copy(AtPu_flightbout_lengths) %>%
  merge(
    AtPu_darkness_propDry,
    by = c("colony", "individ_id", "session_id", "date_time"),
    all.x = TRUE
  ) %>%
  .[, behaviour := data.table::fifelse(
    #is.na(behaviour) &
    !isdaylight & propDarkDry == 1.0 & isdry,
    #isdark  & propDarkDry > 0.9 & isdry,
    "restland",
    behaviour
  ), by = c("colony", "individ_id", "session_id", "BoutNo", "BoutNo_dark")]

AtPu_FlightBouts2 <- data.table::copy(AtPu_activityAdjust1) %>%
  subset(isdry) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, type = "const", fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>% #only keep the recordings at the start of a bout...
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) 

AtPu_FlightBoutLengths2 <- data.table::copy(AtPu_activityAdjust1) %>%
  subset(select = -c(BoutNo)) %>%
  merge(AtPu_FlightBouts2, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo := data.table::nafill(BoutNo, type = "locf"), by = c("colony", "individ_id", "session_id", "isdry")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(isdry, data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo")]

AtPu_activityAdjust2 <- data.table::copy(AtPu_FlightBoutLengths2) %>%
  .[, behaviour := data.table::fcase(
    (!isdark | is.na(isdark) & isdaylight) & isdry & flightLengthMins >= 90,
    "restland",
    !isdark & isdry & flightLengthMins < 90,
    "flight",
    is.na(isdark) & isdaylight & isdry & flightLengthMins < 90,
    sample(c("flight", "restwater"), 1),
    default = behaviour
  ),
  by = c("colony", "individ_id", "session_id", "BoutNo")] %>%
  # .[, behaviour := data.table::fifelse(
  #   #is.na(behaviour) & 
  #     !isdark & isdry,
  #   data.table::fifelse(
  #       flightLengthMins >= 90, #...and is longer than max flight time....
  #       "restland", #...assign as restland...
  #       sample(c("flight", "restwater"), 1) #if shorter than max flight time, assign as flight or restwater
  #     ),
  #   behaviour
  # ), by = c("colony", "individ_id", "session_id", "BoutNo")] %>%
  .[, behaviour := data.table::fifelse(
    is.na(behaviour) & isdaylight & std_light < 1 & std_conductivity < 0.25,
    "in_burrow",
    behaviour)] %>% #in burrow following Annette
  .[,  behaviour := data.table::fifelse(behaviour %in% "in_burrow" &
                                          #if bird is currently in a burrow
                                          (
                                            (
                                              data.table::shift(session_id, n = 1, type = "lag") %in% session_id &
                                                data.table::shift(behaviour, n = 1, type = "lag") %in% "in_burrow" #%in% required over == due to presence of NAs
                                            ) |
                                              #and that bird was in a burrow at the previous sample time
                                              
                                              (
                                                data.table::shift(session_id, n = 1, type = "lead") %in% session_id &
                                                  data.table::shift(behaviour, n = 1, type = "lead") %in% "in_burrow"
                                              ) #or is in a burrow at the next sample time
                                          ), "restland", #assign to being at the 'colony'
                                        behaviour), by = c("individ_id")] %>%
  .[, behaviour := data.table::fifelse(
    is.na(behaviour),
    "restwater",
    behaviour
  )]

AtPu_dry_lengths <- data.table::copy(AtPu_activityAdjust2) %>%
  .[, island := (behaviour == "restland")] %>%
  subset(island) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, gap := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  #find time difference between activity measurements
  .[, gap := data.table::nafill(gap, type = "const", fill = 11)] %>%
  subset(gap > 10) %>%
  .[, LandBoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>%
  subset(select = c(colony, individ_id, session_id, date_time, LandBoutNo)) 

# Sensitivity analysis Scenario
AtPu_reAllocate <- data.table::copy(AtPu_activityAdjust2) %>%
  merge(AtPu_dry_lengths, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, LandBoutNo := data.table::nafill(LandBoutNo, "locf"), by = c("colony", "individ_id", "session_id", "behaviour")] %>% #fill missing bouts with nearest BoutNo
  .[, LandBoutNo := data.table::fifelse(behaviour == "restland", LandBoutNo, NA)] %>%
  .[, DurationLandMins := data.table::fifelse(behaviour == "restland", data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "LandBoutNo")] %>%
  .[, `:=`(
    FirstLand =  data.table::fifelse(behaviour %in% "restland" & !data.table::shift(behaviour %in% "restland", n = 1, type = "lag"), TRUE, FALSE),
    PrevAct =  data.table::fifelse(behaviour %in% "restland" & !data.table::shift(behaviour %in% "restland", n = 1, type = "lag"), data.table::shift(behaviour, n = 1, type = "lag"), NA)
  ), by = c("colony", "individ_id", "session_id")] %>%
  .[, LagMins := data.table::fifelse(behaviour %in% "restland" & !PrevAct %in% "flight" & FirstLand, sample(seq(0,90,10), size = 1, replace = TRUE),0, na = 0), by = c("colony", "individ_id", "session_id", "LandBoutNo")] %>%
  #randomly generate a time period of flight that a bird could hypothetically fly prior to landing
  .[, LagMinsAdj := data.table::fifelse(LagMins >= DurationLandMins & FirstLand, DurationLandMins-10, LagMins), by = c("colony", "individ_id", "session_id", "LandBoutNo")] %>%
  #adjust `LagMins` so that it doesn't exceed the time spent on land
  .[, fill_flight := any(LagMinsAdj >0) & rowid(behaviour) %in% 1:(LagMinsAdj[1]/10), by = c("colony", "individ_id", "session_id", "LandBoutNo")] %>%
  #replace LagMinsAdj future `restland`s with `flight` 
  .[, behaviour := data.table::fifelse(fill_flight, "flight", behaviour)]


AtPu_FlightBouts_lastcheck <- data.table::copy(AtPu_reAllocate) %>%
  subset(behaviour == "flight") %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, type = "const", fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>%
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

AtPu_FlightBoutLengths_lastcheck <- data.table::copy(AtPu_reAllocate) %>%
  subset(select = -c(BoutNo)) %>%
  merge(AtPu_FlightBouts_lastcheck, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo := data.table::nafill(BoutNo, type = "locf"), by = c("colony", "individ_id", "session_id", "behaviour")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(all(behaviour == "flight"), data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "behaviour")]

AtPu_FlightBoutLengths_lastcheck_reclassify <- data.table::copy(AtPu_FlightBoutLengths_lastcheck) %>%
  .[, behaviour := data.table::fifelse(behaviour == "flight" & flightLengthMins > 90 & (!isdark | is.na(isdark) & isdaylight), "restland", behaviour)]

AtPu_FlightBouts_lastcheck2 <- data.table::copy(AtPu_FlightBoutLengths_lastcheck_reclassify) %>%
  subset(behaviour == "flight") %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, type = "const", fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>%
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

AtPu_FlightBoutLengths_lastcheck2 <- data.table::copy(AtPu_FlightBoutLengths_lastcheck_reclassify) %>%
  subset(select = -c(BoutNo)) %>%
  merge(AtPu_FlightBouts_lastcheck, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo := data.table::nafill(BoutNo, "locf"), by = c("colony", "individ_id", "session_id", "behaviour")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(all(behaviour == "flight"), data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo")]

AtPu_FlightBoutLengths_lastcheck_reclassify2 <- data.table::copy(AtPu_FlightBoutLengths_lastcheck2) %>%
  .[, index := data.table::rowid(BoutNo), by = c("colony", "individ_id", "session_id")] %>%
  .[, Activity_change := data.table::fifelse(behaviour == "flight" & flightLengthMins > 90, TRUE, FALSE)] %>%
  .[, Activity_change_duration := data.table::fifelse(Activity_change, (flightLengthMins - 90)/10, 0)] %>%
  .[, behaviour := data.table::fifelse(Activity_change & index <= ceiling(Activity_change_duration), "restwater", behaviour)] %>%
  subset(select = -c(Activity_change, Activity_change_duration))

AtPu_FlightBouts_lastcheck3 <- data.table::copy(AtPu_FlightBoutLengths_lastcheck_reclassify2) %>%
  subset(behaviour == "flight") %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := round(as.numeric(difftime(date_time, date_time_lag, unit=c("mins")))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>%
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

AtPu_FlightBoutLengths_lastcheck3 <- data.table::copy(AtPu_FlightBoutLengths_lastcheck_reclassify2) %>%
  subset(select = -c(BoutNo)) %>%
  merge(AtPu_FlightBouts_lastcheck3, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>% #arrange data by time so we can make lag/lead comparisons
  .[, BoutNo := data.table::nafill(BoutNo, "locf"), by = c("colony", "individ_id", "session_id", "behaviour")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(all(behaviour == "flight"), data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo")] %>%
  .[, Duration := data.table::uniqueN(date_time)*10, by = c("colony", "individ_id", "session_id", "date")] %>%
  subset(Duration >= 1430)

AtPu_daily <- data.table::copy(AtPu_FlightBoutLengths_lastcheck3) %>%
  .[, .(
    tDaylight =  (10 / 60) * sum(isdaylight),
    tForage = sum(behaviour %in% "foraging", na.rm = TRUE) * (10 /
                                                                60),
    tRestWater = sum(behaviour %in% "restwater", na.rm = TRUE) * (10 /
                                                                    60),
    tLand = sum(behaviour %in% "restland", na.rm = TRUE) * (10 /
                                                              60),
    tFlight = sum(behaviour %in% "flight", na.rm = TRUE) * (10 /
                                                              60),
    # tActive = sum(behaviour %in% c("foraging","flight"), na.rm = TRUE) * (10 /
    #                                                         60),
    tActive = sum(behaviour %in% c("foraging"), na.rm = TRUE) * (10 /
                                                                   60)
  ), by = c("colony",
            "individ_id",
            "latest_sex",
            "season",
            "date",
            "sunrise",
            "sunset")] %>%
  .[, tTotal := rowSums(.SD), .SDcols = c("tForage", "tRestWater", "tLand", "tFlight")] %>%
  .[, tDaylight := replace(
    tDaylight,
    (
      season == "breeding" &
        is.na(sunrise) == FALSE & is.na(sunset) == FALSE
    ),
    na.omit(as.numeric(difftime(
      sunset, sunrise, units = "hours"
    )))
  )] %>%
  .[order(date), .SD, by = c("colony", "individ_id")] %>%
  subset(select = c("colony","individ_id","season","date","tDaylight","tForage","tRestWater","tLand","tFlight","tActive","tTotal"))

saveRDS(AtPu_daily, "data/Behaviour/OutputDat_Behav_Energy_AtPu.RDS")

##########################################################################################
#### Process Kittiwake data ####
##########################################################################################

BlKi_light <- readRDS("data/Behaviour/kittiwake_light_data.RDS") %>% #kittiwake data
  data.table::as.data.table()
BlKi_wet <- readRDS("data/Behaviour/kittiwake_immersion_data.RDS") %>% #kittiwake immersion data
  data.table::as.data.table() %>%
  merge(
    #all individuals present in BlKi
    subset(
      BlKi_light,
      select = c(
        colony,
        session_id,
        individ_id,
        date_time,
        year_tracked,
        std_light
      )
    ),
    by = c(
      "session_id",
      "date_time",
      "year_tracked",
      "individ_id",
      "colony"
    ),
    all.x = TRUE
  ) %>%
  merge(
    subset(
      arrival_dates,
      species %in% "BlKi",
      select = c(colony, arrival_day_mean, arrival_month_mean, departure_day_mean, departure_month_mean)
    ),
    by = c("colony"),
    all.x = TRUE
  ) %>%
  .[, `:=`(date = as.Date(date_time),
           time = format(date_time, format = "%H:%M:%S"))] %>%
  .[, `:=`(
    year = lubridate::year(date),
    month = lubridate::month(date),
    day = lubridate::day(date)
  )] %>%
  .[, season := data.table::fifelse(
    (month > arrival_month_mean & month < departure_month_mean) |
      (month == arrival_month_mean & day >= arrival_day_mean) |
      (month == departure_month_mean & day <= departure_day_mean),
    "breeding",
    "nonbreeding"
  )] %>% #classify dates as breeding vs non-breeding and subset to breeding season
  subset(season == "breeding") 

BlKi_colonies <- readRDS("data/Behaviour/kittiwake_position_data.RDS") %>% #puffin colony coordinates
  subset(select = c(colony, col_lon, col_lat)) %>%
  unique()

BlKi_breeding_dates <- data.table::copy(BlKi_wet) %>%
  merge(
    BlKi_colonies,
    by = c("colony"),
    all.x = TRUE,
    sort = FALSE
  ) %>%
  subset(season %in% "breeding",
         select = c("colony", "date", "col_lat", "col_lon")) %>%
  unique() %>% #find unique combinations of colony, day and coordinates
  data.table::setnames(c("col_lat", "col_lon"), c("lat", "lon")) %>% #rename day to date to match `getSunlightTimes` requirements
  .[, as.list(cbind(
    "colony" = colony,
    suncalc::getSunlightTimes(
      data = .,
      tz = "UTC",
      keep = c("sunrise", "sunset")
    )
  ))] #estimate sunrise and sunset times

BlKi_wet <- BlKi_wet %>%
  merge(
    subset(BlKi_breeding_dates, select = c(colony, date, sunrise, sunset)),
    by = c("colony", "date"),
    all.x = TRUE
  ) %>%
  .[, isdaylight :=
      data.table::fcase(
        #assign whether sampling time is during daylight hours
        (season %in% "breeding" &
           (
             is.na(sunrise) == TRUE |
               (date_time >= sunrise - 3600 &
                  date_time <= sunset + 3600) #(?)what about circumstances when std_light > 0
           )),
        TRUE,
        (season %in% "nonbreeding" & std_light > 0),
        TRUE,
        default = FALSE
      )] %>%
  data.table::setorder("colony", "individ_id", "session_id", "date_time") %>%
  .[, `:=`(
    isdry = std_conductivity == 0,
    isdark = std_light <= 0
  )] %>%
  .[, isdark := data.table::fcase(is.na(isdark) & !isdaylight,
                                  TRUE, 
                                  is.na(isdark) & isdaylight,
                                  FALSE,
                                  default = isdark)] %>%
  .[, behaviour :=
      data.table::fcase(
        (std_light > 0 &
           std_conductivity > 0 & std_conductivity < 1),
        "foraging",
        # foraging = daylight and immersion >0% dry and > 0% wet - works all year
        (std_conductivity == 1),
        "restwater"
        #resting on sea surface = immersion 100% wet - works all year
      )]

##########################################################################################
#### Kittiwake energetic classification ####
##########################################################################################

BlKi_flightbouts <- data.table::copy(BlKi_wet) %>%
  subset(isdry) %>%
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := as.numeric(difftime(date_time, date_time_lag, unit=c("mins"))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, fill = 0)] %>%
  subset(timediff == 0 | timediff > 10) %>%
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>%
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

BlKi_flightbout_lengths <- data.table::copy(BlKi_wet) %>%
  merge(BlKi_flightbouts, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  .[, BoutNo := data.table::nafill(BoutNo, "locf"), by = "isdry"] %>%
  .[, flightLengthMins := data.table::fifelse(isdry,length(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo")]

BlKi_activityAdjust1 <- data.table::copy(BlKi_flightbout_lengths) %>%
  .[, behaviour := data.table::fcase(
    is.na(behaviour) & isdry & flightLengthMins <= 240,
    "flight",
    is.na(behaviour) & isdry & flightLengthMins > 240,
    "restland",
    default = behaviour
  )]

BlKi_activityAdjust2_reallocate <- data.table::copy(BlKi_activityAdjust1) %>%
  .[, `:=`(
    FirstLand =  data.table::fifelse(behaviour %in% "restland" & !data.table::shift(behaviour %in% "restland", n = 1, type = "lag"), TRUE, FALSE),
    PrevAct =  data.table::fifelse(behaviour %in% "restland" & !data.table::shift(behaviour %in% "restland", n = 1, type = "lag"), data.table::shift(behaviour, n = 1, type = "lag"), NA)
  ), by = c("colony", "individ_id", "session_id")] %>%
  .[, LandBoutNo := cumsum(FirstLand), by = c("colony", "individ_id", "session_id")] %>%
  .[, LandBoutNo := data.table::fifelse(behaviour == "restland", LandBoutNo, NA)] %>%
  .[, DurationLandMins :=  data.table::fifelse(all(behaviour == "restland"),  data.table::uniqueN(date_time)*10, NA), by = c("colony", "individ_id", "session_id", "LandBoutNo")] %>%
  .[, LagMins := data.table::fifelse(behaviour %in% "restland" & !PrevAct %in% "flight" & FirstLand, sample(seq(0,240,10), size = 1, replace = TRUE),0), by = c("colony", "individ_id", "session_id")] %>%
  .[, LagMinsAdj := data.table::fifelse(LagMins >= DurationLandMins & FirstLand, DurationLandMins-10, LagMins), by = c("colony", "individ_id", "session_id")] %>%
  .[, LandBoutRow := data.table::fifelse(behaviour %in% "restland", rowid(LandBoutNo)*10, 0), by = c("colony", "individ_id", "session_id")] %>%
  .[, behaviour := data.table::fifelse(
    behaviour == "restland" &
      LagMinsAdj[1] > 0 &
      LandBoutRow <= LagMinsAdj[1] & !is.na(LandBoutRow) & 
      PrevAct[1] %in% c(0),
    "flight",
    behaviour
  ), by =  c("colony", "individ_id", "session_id", "LandBoutNo")]

BlKi_FlightBouts_2 <- data.table::copy(BlKi_activityAdjust2_reallocate) %>%
  subset(behaviour == "flight") %>%
  .[, date_time_lag := data.table::shift(date_time, n = 1, type = "lag"),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := as.numeric(difftime(date_time, date_time_lag, unit=c("mins"))),  by = c("colony", "individ_id", "session_id")] %>%
  .[, timediff := data.table::nafill(timediff, fill = 0)] %>% #if no lagged date_time, set difference to 0
  subset(timediff == 0 | timediff > 10) %>%
  .[, BoutNo := seq_len(.N),  by = c("colony", "individ_id", "session_id")] %>% #...and assign all following dry activities to the same bout
  subset(select = c(colony, individ_id, session_id, date_time, BoutNo)) %>%
  unique()

BlKi_FlightBoutLengths_final <- data.table::copy(BlKi_activityAdjust2_reallocate) %>%
  subset(select = -c(BoutNo)) %>%
  merge(BlKi_FlightBouts_2, by = c("colony", "individ_id", "session_id", "date_time"),  all.x = TRUE) %>%
  .[, BoutNo := data.table::nafill(BoutNo, "locf"), by = c("colony", "individ_id", "session_id", "behaviour")] %>% #fill missing bouts with nearest BoutNo
  .[, flightLengthMins := data.table::fifelse(all(behaviour == "flight"), data.table::uniqueN(date_time)*10, 0), by = c("colony", "individ_id", "session_id", "BoutNo", "behaviour")] %>%
  .[, Duration := data.table::uniqueN(date_time)*10, by = c("colony", "individ_id", "session_id", "date")] %>%
  subset(Duration >= 1430)

BlKi_daily <- data.table::copy(BlKi_FlightBoutLengths_final) %>%
  .[, .(
    tDaylight =  (10 / 60) * sum(isdaylight),
    tForage = sum(behaviour %in% "foraging", na.rm = TRUE) * (10 /
                                                                60),
    tRestWater = sum(behaviour %in% "restwater", na.rm = TRUE) * (10 /
                                                                    60),
    tLand = sum(behaviour %in% "restland", na.rm = TRUE) * (10 /
                                                              60),
    tFlight = sum(behaviour %in% "flight", na.rm = TRUE) * (10 /
                                                              60),
    tActive = sum(behaviour %in% c("foraging"), na.rm = TRUE) * (10 /
                                                                   60)
  ), by = c("colony",
            "individ_id",
            "latest_sex",
            "season",
            "date",
            "sunrise",
            "sunset")] %>%
  .[, tTotal := rowSums(.SD), .SDcols = c("tForage", "tRestWater", "tLand", "tFlight")] %>%
  .[, tDaylight := replace(
    tDaylight,
    (
      season == "breeding" &
        is.na(sunrise) == FALSE & is.na(sunset) == FALSE
    ),
    na.omit(as.numeric(difftime(
      sunset, sunrise, units = "hours"
    )))
  )] %>%
  .[order(date), .SD, by = c("colony", "individ_id")]  %>%
  subset(select = c("colony","individ_id","season","date","tDaylight","tForage","tRestWater","tLand","tFlight","tActive","tTotal"))

saveRDS(BlKi_daily, "data/Behaviour/OutputDat_Behav_Energy_BlKi.RDS")
