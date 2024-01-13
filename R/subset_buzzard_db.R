## Derive datasets from buzzard_db, only rerun if new raw data to include

rm(list = ls())
library(tidyverse)
library(magrittr)

load("../../00-Raw/RData/buzzard_db.RData")

df <- filter(buzzard_db$annual_estimates, year > 2006)
range(df$total.attempts)
lm(total.attempts ~ year, df) %>% summary()

## load resigthings
resights <- 
  filter(buzzard_db$resights,
         Ring %in% buzzard_db$ring_db$Ring, # known Ring
         stringr::str_detect(ID, "wingtag"), # with wing-tag
         !is.na(DaysSinceHatching), # known age
         Date < as.Date("2021-04-30"),
         Year <= 2021) # restrict to complete years


## Collect all individuals ringed as juveniles before 2019
## Remove: Individuals died as chicks
## Remove: Individuals without metal ring
ringing_data <- 
  filter(buzzard_db$ring_db,
         Age == "juvenile", ## only locally born
         last == 1,
         !is.na(Residual_weight),
         !is.na(Rank), # brood rank recorded
         #!is.na(Lprev), # endo parasites scored
         Year %in% 2007:2020 # excluding cohort 2021
  ) %>% 
  mutate(year = Year) %>% 
  mutate(Sex = case_when(
    !is.na(Sex) ~ Sex,
    is.na(Sex) ~ Sex_guess 
  )) %>%  
  ## annual effects
  left_join(., buzzard_db$annual_estimates, by = "year") %>% 
  mutate(BCI = Residual_weight %>% as.numeric()) %>% 
  # mutate(Hatch = scale(Hatch) %>% as.numeric()) %>% 
  mutate(Age_days = round(Age_days) %>% as.numeric()) %>% 
  mutate(Sex = ifelse(Sex == 1, "Male", "Female")) %>% 
  mutate(Sex = as.factor(Sex)) %>% 
  mutate(Morph = as.factor(Morph)) %>% 
  mutate(Tag = case_when(
    stringr::str_detect(ID, "wingtag") ~ 1,
    !stringr::str_detect(ID, "wingtag") ~ 0,
    TRUE ~ 0)) %>% 
  subset(., select = c(Ring, Year, Sex, Morph,
                       Rank, Brood_Size, Hatch, Age_days,
                       BCI, Lprev, Tag, Brood_ID,
                       Dead, DateDeath,
                       CauseDeath)) %>% 
  filter(., !CauseDeath %in% c("Mercy killing", "Fell from bag", "Jumper"))
#table(ringing_data$CauseDeath)

resights <- filter(resights, Ring %in% ringing_data$Ring) %>% 
  subset(., select = c(
    Ring, Date, Year, Month, Day, HatchDate, DaysSinceHatching, study_area, Dist2Nest))

readr::write_csv(ringing_data, file = "data/ringing_data.csv", append = F)
readr::write_csv(resights, file =  "data/resighting_data.csv", append = F)

## NAO
djfm <- readr::read_delim(
  file = "https://climatedataguide.ucar.edu/sites/default/files/nao_pc_djfm.txt",
  skip = 1, delim = " ", 
  col_names = c("year", "DJFM")) %>% 
  apply(., 2, trimws) %>% 
  as.data.frame()

djf <- readr::read_delim(
  file = "https://climatedataguide.ucar.edu/sites/default/files/nao_pc_djf.txt",
  skip = 1, delim = " ", 
  col_names = c("year", "DJF")) %>% 
  apply(., 2, trimws) %>% 
  as.data.frame()

mam <- readr::read_delim(
  file = "https://climatedataguide.ucar.edu/sites/default/files/nao_pc_mam.txt",
  skip = 1, delim = " ", 
  col_names = c("year", "MAM")) %>% 
  apply(., 2, trimws) %>% 
  as.data.frame()

jja <- readr::read_delim(
  file = "https://climatedataguide.ucar.edu/sites/default/files/nao_pc_jja.txt",
  skip = 1, delim = " ", 
  col_names = c("year", "JJA")) %>% 
  apply(., 2, trimws) %>% 
  as.data.frame() %>% 
  rbind(., data.frame(year = 2021, JJA = NA))

son <- readr::read_delim(
  file = "https://climatedataguide.ucar.edu/sites/default/files/nao_pc_son.txt",
  skip = 1, delim = " ", 
  col_names = c("year", "SON")) %>% 
  apply(., 2, trimws) %>% 
  as.data.frame() %>% 
  rbind(., data.frame(year = 2021, SON = NA))

nao <- cbind(djfm, DJF = djf$DJF) %>% 
  filter(., year %in% 2007:2021)
#head(nao)
write.table(nao, file = "data/nao.txt", append = F, row.names = F)

## Read climate data in a radius of 60 km around the core of the study area 
## -----------------------------------------------------------------------------
stations <- rdwd::nearbyStations(lat = 52.100, lon = 8.46, radius = 60, mindate = "2021-03-31", quiet = T) 
stations <- stations %>% filter(., res == "daily", var == "kl", von_datum <= as.Date("2007-12-01"), Stationsname != "Bueckeburg")
(stations <- unique.data.frame(stations[, c("Stations_id", "Stationsname", "Stationshoehe", "geoBreite", "geoLaenge")]))

links <- c(unlist(rdwd::selectDWD(name = stations$Stationsname, res = 'daily', var = 'kl', per = 'h')), 
           unlist(rdwd::selectDWD(name = stations$Stationsname, res = 'daily', var = 'kl', per = 'r')))

climate.data <-  lapply(links, function(x) rdwd::dataDWD(x, quiet = T, force = T, overwrite = T, sleep = .05)) %>% 
  do.call('rbind', .) %>% 
  left_join(., stations, by = c("STATIONS_ID" = "Stations_id"))

## format climate.data 
## -----------------------------------------------------------------------------
climate.data[['Month']] <- months(climate.data[['MESS_DATUM']])
climate.data[['Year']] <- substring(climate.data[['MESS_DATUM']],1,4)
climate.data[['Day_of_year']] <- strftime(climate.data[['MESS_DATUM']], format = '%j')
climate.data[['Day']] <- substring(climate.data[['MESS_DATUM']],9,10)

## prepare data frame and omit all NAs
## -----------------------------------------------------------------------------
climate.data <- data.frame(
  station = climate.data$Stationsname,
  lat = climate.data$geoBreite,
  lon = climate.data$geoLaenge,
  elevation = climate.data$Stationshoehe,
  date = as.Date(climate.data$MESS_DATUM),
  month = climate.data$Month,
  year = as.numeric(climate.data$Year),
  day.year = climate.data$Day_of_year,
  day = climate.data$Day,
  precip = climate.data$RSK,
  snow = ifelse(climate.data$SHK_TAG == 0, 0, 1),
  temp = climate.data$TMK
  # cloud = climate.data$NM,
  # humid = climate.data$UPM,
  # max.air = climate.data$TXK,
  # min.air = climate.data$TNK,
  # min.ground = climate.data$TGK,
  # sun = climate.data$SDK
) %>%
  unique.data.frame() %>% 
  filter(., year %in% 2007:2021, month %in% month.name[c(12, 1:3)]) %>% 
  filter(., date %in% as.Date("2007-12-01") : as.Date("2021-03-31"))

table(climate.data$station)

climate.data.raw <- climate.data
climate.data <- na.omit(climate.data)

cat("Missing dat (%) = ", 100 - (nrow(climate.data)/nrow(climate.data.raw)*100))
climate.data <- climate.data.raw

## order month 
## -----------------------------------------------------------------------------
climate.data$month <- factor(climate.data$month, month.name[c(12, 1:3)])
climate.data$time <- ifelse(climate.data$month %in% month.name[1:3], climate.data$year - 1, climate.data$year)

## any missing days?
## -----------------------------------------------------------------------------
expected <- sapply(2007:2020, function(x) {
  as.Date(paste0(x,"-12-01")):as.Date(paste0(x+1,"-03-31"))
}) %>% unlist()
realised <- unique(climate.data$date)
all(expected %in% realised)

## Compute daily means 
## -----------------------------------------------------------------------------
mean_values <- climate.data %>% 
  group_by(date) %>% 
  summarise(
    time = unique(time),
    snow = mean(snow, na.rm = T),
    precip = mean(precip, na.rm = T),
    temp = mean(temp, na.rm = T))

## compute seasonal averages
## -----------------------------------------------------------------------------
seasons <- mean_values %>% 
  group_by(time) %>% 
  summarise(snow = sum(snow),
            precip = sum(precip),
            temp = mean(temp))

climate.data <- list(raw = climate.data, season = seasons)
saveRDS(climate.data, "data/climate.data.RDS", compress = "bzip2")
#save(climate.data, file = "data/climate.data.RData")

