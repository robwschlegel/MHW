#############################################################################
## This script does:
# 1. loads annual MHW/ MCS files;
# 2. calculates event count, length and mean intensity for each coast for both datasets;
# 3. Load all events and extract top three MHWs/ MCSs;
# 4. Extracts top three MHW/ MCS events for each coastal section and type of data;
# 5. Extracts top 1 MHW/ MCS events for each coastal section and type of data;
# 6. Calculate beginning and end dates for largest events;
# 7. Calculate co-occurrence between coastal sections

#############################################################################
## DEPENDS ON:
require(zoo); require(plyr); require(stringr); library(lubridate)
source("setupParams/theme.R")
# "graph/eventsPlots2.R" # This script calculates the co-occurrence rates for sites
# "data/metaData2.csv"
# All of the files from Eric
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "data/allAnnualResults.csv"
# "data/allAnnualSSTResults.csv"
# "data/mhw3.csv"
# "data/mcs3.csv"
# "data/mhwSST3.csv"
# "data/mcsSST3.csv"
# "data/mhw1dates.csv"
# "data/mcs1dates.csv"
#############################################################################

#############################################################################
## 1. loads annual and event MHW/ MCS files and metadata2

# First specify coastal sections
wc <- c("Hout Bay", "Kommetjie", "Port Nolloth", "Sea Point")
sc <- c("Fish Hoek", "Gordons Bay", "Hamburg", "Hermanus", "Humewood", "Knysna", 
        "Mossel Bay", "Muizenberg", "Pollock Beach", "Tsitsikamma West", 
        "Storms River Mouth", "Tsitsikamma East", "Ystervarkpunt")
ec <- c("Eastern Beach", "Nahoon Beach", "Orient Beach", "Sodwana")

# Then specify directories for loading
# Annual values
dir1 <- paste(getwd(), "/data/MHW/annual/", sep = "")
dir2 <- paste(getwd(), "/data/MCS/annual/", sep = "")
dir3 <- paste(getwd(), "/data/MHW/SST annual/", sep = "")
dir4 <- paste(getwd(), "/data/MCS/SST annual/", sep = "")
#  Individual events
dir5 <- paste(getwd(), "/data/MHW/events/", sep = "")
dir6 <- paste(getwd(), "/data/MCS/events/", sep = "")
dir7 <- paste(getwd(), "/data/MHW/SST events/", sep = "")
dir8 <- paste(getwd(), "/data/MCS/SST events/", sep = "")

# Load metadata
load("data/metaData2.Rdata")

# Load annual event stats
annualLoad <- function(dir) {
  fname1 <-  dir(dir, full.names = TRUE)
  fname2 <-  dir(dir, full.names = FALSE)
  pf <- str_sub(fname2, -20, -1)
  siteNames1 <-  unlist(strsplit(dir(dir, full.names = FALSE), pf[1]))
  siteNames1 <-  str_replace_all(siteNames1, "_", " ") # parse names for site column
  dat <- data.frame()
  for(i in 1:length(fname1)) { # A shameful for loop... in order to label sites correctly
    x <- read.csv(fname1[i], header = TRUE, skip = 2)
    x$site <-  siteNames1[i]
    if(x$site[1] %in% wc) {
      x$coast <- "wc"
    } else if(x$site[1] %in% sc) {
      x$coast <- "sc"
    } else if(x$site[1] %in% ec) {
      x$coast <- "ec"
    }
    x <- x[,c(25:26,1:24)]
    dat <- rbind(dat, x)
  }
  dat$coast <- factor(dat$coast, levels = c("wc", "sc", "ec"))
  return(dat)
}

mhwAnnual <- annualLoad(dir1)
mcsAnnual <- annualLoad(dir2)
mhwAnnualSST <- annualLoad(dir3)
mcsAnnualSST <- annualLoad(dir4)

#############################################################################
## 2. calculates event count, length and mean intensity for each coast for both datasets

# Function used for calculations
resultsAnnualCoastal <- function(mhw1, mcs1){ # To be used with "annual" data frames only
  results <- data.frame(coast = as.factor("All"), 
                        mhw_count = round(mean(mhw1[, 4], na.rm = T), 1), 
                        mhw_length = round(mean(mhw1[, 9], na.rm = T), 1),
                        mhw_intensity = round(mean(mhw1[, 7], na.rm = T), 2),
                        mcs_count = round(mean(mcs1[, 4], na.rm = T), 1), 
                        mcs_length = round(mean(mcs1[, 9], na.rm = T), 1),
                        mcs_intensity = round(mean(mcs1[, 7], na.rm = T), 2))
  for(i in 1:length(levels(mhw1$coast))){
    mhw2 <- droplevels(subset(mhw1, coast == levels(mhw1$coast)[i]))
    mcs2 <- droplevels(subset(mcs1, coast == levels(mcs1$coast)[i]))
    z <- data.frame(coast = levels(mhw1$coast)[i], 
                    mhw_count = round(mean(mhw2[, 4], na.rm = T), 1), 
                    mhw_length = round(mean(mhw2[, 9], na.rm = T), 1),
                    mhw_intensity = round(mean(mhw2[, 7], na.rm = T), 2),
                    mcs_count = round(mean(mcs2[, 4], na.rm = T), 1), 
                    mcs_length = round(mean(mcs2[, 9], na.rm = T), 1),
                    mcs_intensity = round(mean(mcs2[, 7], na.rm = T), 2))
    results <- rbind(results, z)
  }
  return(results)
}

allAnnualResults <- resultsAnnualCoastal(mhwAnnual, mcsAnnual)
write.csv(allAnnualResults, "data/allAnnualResults.csv")
allAnnualSSTResults <- resultsAnnualCoastal(mhwAnnualSST, mcsAnnualSST)
write.csv(allAnnualSSTResults, "data/allAnnualSSTResults.csv")

#############################################################################
## 3. Load all events and extract top three MHWs/ MCSs

##### 
#These hashes are here so that this chunk can be collapsed efficiently
colNames <- c(
  "eventNo", # Event number,
  "yearStrt", # Start year,
  "monthStrt", # Start month,
  "dayStrt", # Start day,
  "yearPk", # Peak year,
  "monthPk", # Peak month,
  "dayPk", # Peak day,
  "yearEnd", # End year,
  "monthEnd", # End month,
  "dayEnd", # End day,
  "duration", # Duration [days],
  "intMax", # Maximum intensity [deg C],
  "intMean", # Mean intensity [deg C],
  "intCum", # Cumulative intensity [deg C x days],
  "intVar", # Intensity variability [deg C],
  "onsetRt", # Rate of onset [deg C / days],
  "declRt", # Rate of decline [deg C / days],
  "relThreshMax", # Maximum intensity (rel. thresh.) [deg C],
  "relThreshMean", # Mean intensity (rel. thresh.) [deg C],
  "relThreshCum", # Cumulative intensity (rel. thresh.) [deg C x days],
  "relThreshVar", # Intensity variability (rel. thresh.) [deg C],
  "absIntMax", # Maximum intensity (absolute) [deg C],
  "absIntMean", # Mean intensity (absolute) [deg C],
  "absIntCum", # Cumulative intensity (absolute) [deg C x days],
  "absIntVar", # Intensity variability (absolute) [deg C],
  "normIntMax", # Maximum intensity (normalized) [unitless],
  "normIntMean" # Mean intensity (normalized) [unitless]
)
#####

# Load events
eventLoad <- function(dir) {
  fname1 <-  dir(dir, full.names = TRUE)
  fname2 <-  dir(dir, full.names = FALSE)
  pf <- str_sub(fname2, -20, -1)
  siteNames1 <-  unlist(strsplit(dir(dir, pattern = "data.events", full.names = FALSE), pf[1]))
  siteNames1 <-  str_replace_all(siteNames1, "_", " ") # parse names for site column
  dat <- data.frame()
  for(i in 1:length(fname1)) { # A shameful for loop... in order to label sites correctly
    x <- read.csv(fname1[i], header = FALSE, skip = 2, sep = ",",
                  col.names = colNames)
    x$site <-  siteNames1[i]
    x$coast <- metaData2$coast[metaData2$site == x$site[1]]
    x$date <-  as.Date(paste(x$yearStrt, x$monthStrt, x$dayStrt, sep = "-"))
    x$month <- floor_date(as.Date(paste(x$yearStrt, x$monthStrt, x$dayStrt, sep = "-")), "month")
    x$lon <- metaData2$lon[metaData2$site == x$site[1]]
    x$lat <- metaData2$lat[metaData2$site == x$site[1]]
    dat <- rbind(dat, x)
  }
  return(dat)
}

mhwEvent <- eventLoad(dir5)
mcsEvent <- eventLoad(dir6)
mhwEventSST <- eventLoad(dir7)
mcsEventSST <- eventLoad(dir8)

# The top three events per site
eventLoadn <- function(dir, nCum = 5) {
  fname1 = dir(dir, pattern = "data.events", full.names = TRUE)
  fname2 = dir(dir, pattern = "data.events", full.names = FALSE)
  pf <- str_sub(fname2, -20, -1)
  l1 = llply(fname1, read.csv, skip = 2, col.names = colNames) # loops avoided!
  l1 = llply(l1, arrange, -abs(intCum)) # sort by intCum in descending order
  top_n = function(x) head(x, nCum)
  df1 = ldply(l1, top_n) # pick the n highest ones
  siteNames1 = unlist(strsplit(dir(dir, pattern = "data.events", full.names = FALSE), pf[1]))
  siteNames1 = str_replace_all(siteNames1, "_", " ") # parse names for site column
  df1$site = rep(siteNames1, each = nCum)
  df1$month = floor_date(as.Date(paste(df1$yearStrt, df1$monthStrt, df1$dayStrt, sep = "-")), "month")
  return(df1)
}

# in situ
mhwn <- eventLoadn(dir5, nCum = 3) # produce a data frame with mhw data...
mhwn$event <- rep(1:3, length(levels(as.factor(mhwn$site)))) # Make sure to correct the rep() to match nCum
mhwn$type <- "insitu"
mcsn <- eventLoadn(dir6, nCum = 3) # produce a data frame with mcs data...
mcsn$event <- rep(1:3, length(levels(as.factor(mcsn$site))))
mcsn$type <- "insitu"
# SST
mhwnSST <- eventLoadn(dir7, nCum = 3)
mhwnSST$event <- rep(1:3, length(levels(as.factor(mhwnSST$site))))
mhwnSST$type <- "sst"
mcsnSST <- eventLoadn(dir8, nCum = 3)
mcsnSST$event <- rep(1:3, length(levels(as.factor(mcsnSST$site))))
mcsnSST$type <- "sst"

#############################################################################
## 4. Extracts top three MHW/ MCS events for each coastal section and type of data

# The extracting function
topCoastn <- function(dat, nCum){
  wcoast <- arrange(droplevels(dat[dat$site %in% wc,]), -abs(intCum))
  wcoast$coast <- "wc"
  scoast <- arrange(droplevels(dat[dat$site %in% sc,]), -abs(intCum))
  scoast$coast <- "sc"
  ecoast <- arrange(droplevels(dat[dat$site %in% ec,]), -abs(intCum))
  ecoast$coast <- "ec"
  result <- rbind(head(wcoast, nCum), head(scoast, nCum), head(ecoast, nCum))
  return(result)
}

mhw3 <- topCoastn(mhwn, 3)
mhw3 <- mhw3[,c(28:29,4,30:32,11,13:14)]
write.csv(mhw3, "data/mhw3.csv")
mcs3 <- topCoastn(mcsn, 3)
mcs3 <- mcs3[,c(28:29,4,30:32,11,13:14)]
write.csv(mcs3, "data/mcs3.csv")
mhwSST3 <- topCoastn(mhwnSST, 3)
mhwSST3 <- mhwSST3[,c(28:29,4,30:32,11,13:14)]
write.csv(mhwSST3, "data/mhwSST3.csv")
mcsSST3 <- topCoastn(mcsnSST, 3)
mcsSST3 <- mcsSST3[,c(28:29,4,30:32,11,13:14)]
write.csv(mcsSST3, "data/mcsSST3.csv")

#############################################################################
## 5. Extracts top 1 MHW/ MCS events for each coastal section and type of data

mhw1 <- topCoastn(mhwn, 1)
mcs1 <- topCoastn(mcsn, 1)
mhwSST1 <- topCoastn(mhwnSST, 1)
mcsSST1 <- topCoastn(mcsnSST, 1)

#############################################################################
## 6. Calculate beginning and end dates for largest events

eventDates <- function(dat){
  dat$start.date <- as.Date(paste(dat$yearStrt, dat$monthStrt, dat$dayStrt, sep = "/"), "%Y/%m/%d")
  dat$end.date <- as.Date(paste(dat$yearEnd, dat$monthEnd, dat$dayEnd, sep = "/"), "%Y/%m/%d")
  dat <- dat[,c(28,31:34,11,13:14)]
  return(dat)
}

mhw1dates <- eventDates(mhw1)
write.csv(mhw1dates, "data/mhw1dates.csv")
mcs1dates <- eventDates(mcs1)
write.csv(mcs1dates, "data/mcs1dates.csv")
mhwSST1dates <- eventDates(mhwSST1)
mcsSST1dates <- eventDates(mcsSST1)

#############################################################################
# ### Test draft of what figure 4 will look like
# # load data
# load("data/SACTNdaily_v4.0.Rdata")
# load("data/OISSTdaily.Rdata")
# 
# # Subset data
# mhw1ts <- droplevels(subset(SACTNdaily_v4.0, site %in% mhw1dates$site))[,c(1,3:4)] # Remove columns to match OISST
# mcs1ts <- droplevels(subset(SACTNdaily_v4.0, site %in% mcs1dates$site))[,c(1,3:4)]
# mhwSST1ts <- droplevels(subset(OISSTdaily, site %in% mhw1dates$site)) # in situ sites used intentionally to select OISST
# mcsSST1ts <- droplevels(subset(OISSTdaily, site %in% mcs1dates$site))
# 
# # Subset dates
# subsetDates <- function(dat1, dat2) {
#   dat3 <- data.frame()
#   for(i in 1:length(levels(dat1$site))){
#     dat4 <- droplevels(subset(dat1, site == levels(dat1$site)[i]))
#     dat5 <- droplevels(subset(dat2, site == levels(dat1$site)[i]))
#     dates <- seq(dat5$start.date, dat5$end.date, 1)
#     dat6 <- droplevels(subset(dat4, as.Date(date) %in% dates))
#     dat6$coast <- dat5$coast
#     dat3 <- rbind(dat3, dat6) 
#   }
#   return(dat3)
# }
# 
# mhw1ts <- subsetDates(mhw1ts, mhw1dates)
# mhw1ts$type <- as.factor("in situ")
# mhw1ts$event <- as.factor("mhw")
# mcs1ts <- subsetDates(mcs1ts, mcs1dates)
# mcs1ts$type <- as.factor("in situ")
# mcs1ts$event <- as.factor("mcs")
# mhwSST1ts <- subsetDates(mhwSST1ts, mhw1dates) # The in situ dates are used here intentionally to extract OISST data
# mhwSST1ts$type <- as.factor("OISST")
# mhwSST1ts$event <- as.factor("mhw")
# mcsSST1ts <- subsetDates(mcsSST1ts, mcs1dates)
# mcsSST1ts$type <- as.factor("OISST")
# mcsSST1ts$event <- as.factor("mcs")
# 
# # Combine for plotting
# all1ts <- rbind(mhw1ts, mcs1ts, mhwSST1ts, mcsSST1ts)
# # Remove east coast sites
# all1ts <- droplevels(all1ts[all1ts$coast != "ec", ])
# 
# # Create figure
# test <- ggplot(data = all1ts, aes(x = date, y = temp)) + bw_update +
#   geom_line() +
#   facet_wrap(event + coast ~ type, ncol = 2, scale = "free_x")
# test
# ggsave("graph/figure4demo.pdf")

#############################################################################
## 7.Calculate co-occurrence between coastal sections

## NB: The way in which this is calculated for the coast needs to be considered carefully...
    ## The following function is too biased for the south coast.
    ## Rather using the indivudal co-occurrence rates and meaning them by coast.

# cooccurrenceCoast <- function(dat1, dat2, lag = seq(2,14,2)){
#   dat3 <- data.frame()
#   direction <- c("b","x","a")
#   coasts <- factor(matrix(c("west", "south", "east")), levels = c("west", "south", "east"))
#   for(i in 1:length(levels(coasts))) {
#     x1 <- droplevels(subset(dat1, coast == coasts[i]))
#     x2 <- droplevels(subset(dat2, coast == coasts[i]))
#     x1 <- x1[x1$yearStrt >= min(x2$yearStrt), ] # Subset x so that dates match up
#     x1 <- x1[x1$yearStrt <= max(x2$yearStrt), ]
#     x2 <- x2[x2$yearStrt >= min(x1$yearStrt), ]
#     x2 <- x2[x2$yearStrt <= max(x1$yearStrt), ]
#     for(j in 1:length(lag)){
#       for(k in 1:length(seq(0.0,1,0.1))){
#         for(l in 1:length(direction)){
#           x1.1 <- x1[x1$intCum >= quantile(x1$intCum, probs = seq(0.0,1,0.1)[k]),]
#           x2.1 <- x2[x2$intCum >= quantile(x2$intCum, probs = seq(0.0,1,0.1)[k]),]
#           y <- 0
#           #x3 <- data.frame() # For test purposes to see which events match up
#           for(m in 1:nrow(x1.1)) {
#             x1.2 <- x1.1$date[m]
#             if(direction[l] == "b"){
#               x1.3 <- seq((x1.2 - days(lag[j])), x1.2, 1)
#             } else if(direction[l] == "x"){
#               x1.3 <- seq((x1.2 - days(lag[j])), (x1.2 + days(lag[j])), 1)
#             } else if (direction[l] == "a") {
#               x1.3 <- seq(x1.2, (x1.2 + days(lag[j])), 1)
#             }
#             x2.2 <- droplevels(subset(x2.1, date %in% x1.3))
#             y <- y + nrow(x2.2)
#           }
#           z <- data.frame(coast = x1$coast[1], lon = x1$lon[1], lat = x1$lat[1],
#                           lag = lag[j], quantile = seq(0.0,1,0.1)[k], direction = direction[l],
#                           insitu = nrow(x1.1), OISST = nrow(x2.1),
#                           cooccurrence= y, proportion = y/nrow(x1.1))
#           dat3 <- rbind(dat3, z)
#         }
#       }
#     }
#   }
#   return(dat3)
# }

#mhwCO0 <- cooccurrence(mhw, mhwSST, lag = 0) # Test to see which happen on exact same day
#MHW
# mhwCoastCO <- cooccurrenceCoast(mhwEvent, mhwEventSST) # This takes several minutes to run... cursed for loops...
# write.csv(mhwCoastCO, "data/mhwCoastCO.csv", row.names = F)
# mhwCoastCO <- read.csv("data/mhwCoastCO.csv")
#MCS
# mcsCoastCO <- cooccurrenceCoast(mcsEvent, mcsEventSST)
# write.csv(mcsCoastCO, "data/mcsCoastCO.csv", row.names = F)
# mcsCoastCO <- read.csv("data/mcsCoastCO.csv")

mhwCoastCO <- read.csv("data/mhwCO.csv")
mcsCoastCO <- read.csv("data/mcsCO.csv")

# Add coast column
mhwCoastCO$coast <- "x"
mhwCoastCO$coast[mhwCoastCO$site %in% wc] <- "west"
mhwCoastCO$coast[mhwCoastCO$site %in% sc] <- "south"
mhwCoastCO$coast[mhwCoastCO$site %in% ec] <- "east"

mcsCoastCO$coast <- "x"
mcsCoastCO$coast[mcsCoastCO$site %in% wc] <- "west"
mcsCoastCO$coast[mcsCoastCO$site %in% sc] <- "south"
mcsCoastCO$coast[mcsCoastCO$site %in% ec] <- "east"

mhwCoastCO <- ddply(mhwCoastCO, .(coast, lag, quantile, direction), summarize, 
                    proportion = mean(proportion, na.rm = TRUE))
write.csv(mhwCoastCO, "data/mhwCoastCO.csv")

mcsCoastCO <- ddply(mcsCoastCO, .(coast, lag, quantile, direction), summarize, 
                    proportion = mean(proportion, na.rm = TRUE))
write.csv(mcsCoastCO, "data/mcsCoastCO.csv")

# Some exploratory states
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="b"])
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="x"])
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="a"])

mean(mcsCoastCO$proportion[mcsCoastCO$direction =="b"])
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="x"])
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="a"])

#############################################################################
## 8. Additional analyses/ computations

# Quantify the occurrence of the top three MHWs and MCSs tper coast
  # This can be used to infer climate change

