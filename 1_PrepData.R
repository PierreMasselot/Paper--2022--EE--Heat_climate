#################################################################
#
#           AUDACE - Summer mortality prediction
#               Part 1: data preparation
#
#################################################################

#---------------------------------------
#  Daily data
#---------------------------------------

dlist <- lapply(region, function(reg){
  # Load data
  daily_dat <- read.table(sprintf("Data/%sdata.csv", reg), 
    header = T, sep = ",")
  
  # Create time variables
  daily_dat$Date <- as.Date(daily_dat$Date, 
    tryFormats = c("%d/%m/%Y", "%Y-%m-%d"))
  daily_dat$mon <- as.numeric(format(daily_dat$Date, "%m"))
  daily_dat$year <- as.numeric(format(daily_dat$Date, "%Y"))
  daily_dat$doy <- as.numeric(format(daily_dat$Date, "%j"))
  daily_dat$dow <- as.factor(weekdays(daily_dat$Date))
  
  # Select hot months
  daily_dat <- daily_dat[daily_dat$mon %in% hot_months,]
  
  # Number of years
  years <- unique(daily_dat$year)
  ny <- length(years)

  # Daily mean temperature
  daily_dat$Tmean <- with(daily_dat, (Tmin + Tmax) / 2)
  
  # Output
  daily_dat
})
names(dlist) <- region

# Obtain unique years
years <- sort(unique(unlist(lapply(dlist, "[[", "year"))))
ny <- length(years)
  
#---------------------------------------
#  Climate data
#---------------------------------------

# Read climate data
climate_dat <- read.table("Data/climat.csv", header = T, sep = ",")

# Lag climate
lagged_climate <- lapply(climate_dat[,climate_labs], Lag, clag)

# Select years
whichlines <- (climate_dat$Y %in% years) & (climate_dat$M == min(hot_months))
fclimate_dat <- lapply(lagged_climate, "[", whichlines,)
for (j in seq_along(fclimate_dat)){
  colnames(fclimate_dat[[j]]) <- paste(mon_vec, year_vec, sep = "-")
}

#----- Compute annual summer temperature

# Compute average by year
aggtemp <- lapply(dlist, function(d) 
  aggregate(Tmean ~ year, data = d, FUN = mean))

# Merge them
summertemp <- merge(aggtemp[[1]], aggtemp[[2]], by = "year", all = T,
  suffixes = sprintf("_%s", names(aggtemp)))

# Sort
summertemp <- summertemp[order(summertemp$year),]

#---------------------------------------
#  Final data.frame
#---------------------------------------

datatab <- c(Year = list(years), fclimate_dat,
  s = list(clag - 0.5), summertemp[,-1])

#---------------------------------------
#  Geographical data
#---------------------------------------

#----- Metropolitan areas boundaries

# Locate and unzip shapefiles
source <- paste0("C:/Users/PierreMasselot/Documents/Recherche", 
  "/2012-2020 - INRS/2019-2020 - AUDACE/Données/geo/SHP.zip")
unzip(zipfile = source, exdir = getwd())

# Load city boundaries
citylims <- st_read(paste0("SHP/munic_s.shp"))

# load metropolitan areas boudaries
cmlims <- st_read(paste0("SHP/comet_s.shp"))

# Remove temporary files
unlink("SHP", recursive = T)

# Intersect them
geolims <- list()
for (i in 1:nrow(cmlims)) geolims[[i]] <- st_intersection(citylims, cmlims[i,])

# Centroids
geocenter <- st_centroid(cmlims)
