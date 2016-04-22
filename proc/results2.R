#############################################################################
## This script does:
# 1. load annual MHW/ MCS files and clean meta-data for publication;
# 2. Calculates event count, length and mean intensity for each coast for both datasets;
# 3. Load all events and extract top three MHWs/ MCSs;
# 4. Extracts top three MHW/ MCS events for each coastal section and type of data;
# 5. Extracts top 1 MHW/ MCS events for each coastal section and type of data;
# 6. Calculate beginning and end dates for largest events;
# 7. Calcuate statistical significance between coastal sections;
# 8. Calculate co-occurrence between coastal sections;
# 9. Calcuate stats and statistical significance between coastal sections etc. for co-occurrence;
# 10. Rate of increase in MHWs/ MCSs;
# 11. R2 between in situ and OISST time series;
# 12. Co-occurrence within datasets and coastal sections

#############################################################################
## DEPENDS ON:
require(zoo); require(plyr); require(stringr); require(lubridate); require(xtable); library(magrittr)
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
    x <- x[,c(25:26,1,2,3,5)]
    colnames(x)[3:6] <- c("year","frequency","duration","intensity")
    x[is.na(x$frequency)] <- 0
    dat <- rbind(dat, x)
  }
  dat$coast <- factor(dat$coast, levels = c("wc", "sc", "ec"))
  return(dat)
}

mhwAnnual <- annualLoad(dir1)
mcsAnnual <- annualLoad(dir2)
mhwAnnualSST <- annualLoad(dir3)
mcsAnnualSST <- annualLoad(dir4)

# Load metadata
load("data/metaData2.Rdata")

# Clean Metadata for inclusion as a table in the paper
metaData3 <- metaData2
metaData3$site <- as.character(metaData3$site)
metaData3$start.date <- as.character(metaData3$start.date)
metaData3$end.date <- as.character(metaData3$end.date)
row.names(metaData3) <- NULL
metaData3 <- metaData3[c(1:4,22,5:17,23,18:21,24),c(1:5,7:9,11:15)]
metaData3[c(5,19,24),1] <- "coast"
metaData3[c(5,19,24),2] <- "mean"
metaData3[c(5,19,24),c(4:7)] <- NA
metaData3$length <- round(metaData3$length/365,1) # Conert to year
colnames(metaData3)[c(6:9)] <- c("start date", "end date", "duration (years)", "NA %")
xtable(metaData3, auto = TRUE)

#############################################################################
## 2. Calculates event count, length and mean intensity for each coast for both datasets

# Function used for calculations
resultsAnnualCoastal <- function(mhw1, mcs1){ # To be used with "annual" data frames only
  results <- data.frame(coast = as.factor("All"),
                        mhw_freq = round(mean(mhw1[, 4], na.rm = T), 1),
                        mhw_freq_sd = round(sd(mhw1[, 4], na.rm = T), 1),
                        mhw_dur = round(mean(mhw1[, 5], na.rm = T), 1),
                        mhw_dur_sd = round(sd(mhw1[, 5], na.rm = T), 1),
                        mhw_intens = round(mean(mhw1[, 6], na.rm = T), 2),
                        mhw_intens_sd = round(sd(mhw1[, 6], na.rm = T), 2),
                        mcs_freq = round(mean(mcs1[, 4], na.rm = T), 1),
                        mcs_freq_sd = round(sd(mcs1[, 4], na.rm = T), 1),
                        mcs_dur = round(mean(mcs1[, 5], na.rm = T), 1),
                        mcs_dur_sd = round(sd(mcs1[, 5], na.rm = T), 1),
                        mcs_intens = round(mean(mcs1[, 6], na.rm = T), 2),
                        mcs_intens_sd = round(sd(mcs1[, 6], na.rm = T), 2))
  for(i in 1:length(levels(mhw1$coast))){
    mhw2 <- droplevels(subset(mhw1, coast == levels(mhw1$coast)[i]))
    mcs2 <- droplevels(subset(mcs1, coast == levels(mcs1$coast)[i]))
    z <- data.frame(coast = levels(mhw1$coast)[i],
                    mhw_freq = round(mean(mhw2[, 4], na.rm = T), 1),
                    mhw_freq_sd = round(sd(mhw2[, 4], na.rm = T), 1),
                    mhw_dur = round(mean(mhw2[, 5], na.rm = T), 1),
                    mhw_dur_sd = round(sd(mhw2[, 5], na.rm = T), 1),
                    mhw_intens = round(mean(mhw2[, 6], na.rm = T), 2),
                    mhw_intens_sd = round(sd(mhw2[, 6], na.rm = T), 2),
                    mcs_freq = round(mean(mcs2[, 4], na.rm = T), 1),
                    mcs_freq_sd = round(sd(mcs2[, 4], na.rm = T), 1),
                    mcs_dur = round(mean(mcs2[, 5], na.rm = T), 1),
                    mcs_dur_sd = round(sd(mcs2[, 5], na.rm = T), 1),
                    mcs_intens = round(mean(mcs2[, 6], na.rm = T), 2),
                    mcs_intens_sd = round(sd(mcs2[, 6], na.rm = T), 2))
    results <- rbind(results, z)
  }
  return(results)
}

allAnnualResults <- resultsAnnualCoastal(mhwAnnual, mcsAnnual)
xtable(allAnnualResults)
write.csv(allAnnualResults, "data/allAnnualResults.csv")
allAnnualSSTResults <- resultsAnnualCoastal(mhwAnnualSST, mcsAnnualSST)
xtable(allAnnualSSTResults)
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
  l1 = llply(l1, plyr::arrange, -abs(intCum)) # sort by intCum in descending order
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
mhwn$index <- rep(1:3, length(levels(as.factor(mhwn$site)))) # Make sure to correct the rep() to match nCum
mhwn$type <- "insitu"
mcsn <- eventLoadn(dir6, nCum = 3) # produce a data frame with mcs data...
mcsn$index <- rep(1:3, length(levels(as.factor(mcsn$site))))
mcsn$type <- "insitu"
# SST
mhwnSST <- eventLoadn(dir7, nCum = 3)
mhwnSST$index <- rep(1:3, length(levels(as.factor(mhwnSST$site))))
mhwnSST$type <- "OISST"
mcsnSST <- eventLoadn(dir8, nCum = 3)
mcsnSST$index <- rep(1:3, length(levels(as.factor(mcsnSST$site))))
mcsnSST$type <- "OISST"

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
  result <- result[,c(32,28:29,4,11,13:14)]
  result$month <- format(result$month, "%Y-%m")
  result$'start date' <- paste(result$month, result$dayStrt, sep = "-")
  result <- result[,c(1,2,8,5:7)]
  return(result)
}

mhw3 <- topCoastn(mhwn, 3)
xtable(mhw3)
write.csv(mhw3, "data/mhw3.csv")
mcs3 <- topCoastn(mcsn, 3)
xtable(mcs3)
write.csv(mcs3, "data/mcs3.csv")
mhwSST3 <- topCoastn(mhwnSST, 3)
xtable(mhwSST3)
write.csv(mhwSST3, "data/mhwSST3.csv")
mcsSST3 <- topCoastn(mcsnSST, 3)
xtable(mcsSST3)
write.csv(mcsSST3, "data/mcsSST3.csv")

# Combine into one table for publication
  # Currently keeping them seperate as I think it looks less cluttered

#############################################################################
## 5. Extracts top 1 MHW/ MCS events for each coastal section and type of data

mhw1 <- topCoastn(mhwn, 1)
mcs1 <- topCoastn(mcsn, 1)
mhwSST1 <- topCoastn(mhwnSST, 1)
mcsSST1 <- topCoastn(mcsnSST, 1)

#############################################################################
## 6. Calculate beginning and end dates for largest events

## NB: The code from which this is calculated was changed so this no longer runs, but it is unneccessary anyway...

# eventDates <- function(dat){
#   dat$start.date <- as.Date(paste(dat$yearStrt, dat$monthStrt, dat$dayStrt, sep = "/"), "%Y/%m/%d")
#   dat$end.date <- as.Date(paste(dat$yearEnd, dat$monthEnd, dat$dayEnd, sep = "/"), "%Y/%m/%d")
#   dat <- dat[,c(28,31:34,11,13:14)]
#   return(dat)
# }
#
# mhw1dates <- eventDates(mhw1)
# write.csv(mhw1dates, "data/mhw1dates.csv")
# mcs1dates <- eventDates(mcs1)
# write.csv(mcs1dates, "data/mcs1dates.csv")
# mhwSST1dates <- eventDates(mhwSST1)
# mcsSST1dates <- eventDates(mcsSST1)

#############################################################################
## 7. Calcuate statistical significance between coastal sections

# A function that takes a column of data and tests for homoscedasticity and normality
# x <- mhwAnnual[6]
# assumptions <- function(x) {
#   var <- round(var(x, na.rm = T), 2)
#   nor <- round(shapiro.test(x)$p.value, 2)
#   result <- data.frame(var, nor)
#   return(result)
# }

# Check assumptions for all data
# test <- assumptions(x)
# test <- daply(mhwAnnual[mhwAnnual$coast == "wc",][,4:5],
#               .variables = colnames(mhwAnnual[,4:5]), .fun = mean())

# Combine annual data
allAnnual <- rbind(mhwAnnual, mcsAnnual)
allAnnual$event <-rep(c("mhw", "mcs"), each = nrow(mhwAnnual))
allAnnual$type <- "insitu"
allAnnualSST <- rbind(mhwAnnualSST, mcsAnnualSST)
allAnnualSST$event <-rep(c("mhw", "mcs"), each = nrow(mhwAnnualSST))
allAnnualSST$type <- "OISST"
allAllAnnual <- rbind(allAnnual, allAnnualSST)

# Combine event data
mhwEvent$event <- "mhw"
mcsEvent$event <- "mcs"
allEvent <- rbind(mhwEvent, mcsEvent)
allEvent$type <- "insitu"
mhwEventSST$event <- "mhw"
mcsEventSST$event <- "mcs"
allEventSST <- rbind(mhwEventSST, mcsEventSST)
allEventSST$type <- "OISST"
allAllEvent <- rbind(allEvent, allEventSST)

## Mean(sd) of cummulative intensity for the entire coastline
# in situ + MHW
mean(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mhw"], na.rm = T)
sd(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mhw"], na.rm = T)
# in situ + MCS
mean(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mcs"], na.rm = T)
sd(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mcs"], na.rm = T)
# OISST + MHW
mean(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mhw"], na.rm = T)
sd(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mhw"], na.rm = T)
# OISST + MCS
mean(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mcs"], na.rm = T)
sd(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mcs"], na.rm = T)

## ANOVA for in situ
aovFrequency <- aov(frequency ~ coast * event, data = allAllAnnual[allAllAnnual$type == "insitu",])
#summary(aovFrequency)
tukeyFrequency <- TukeyHSD(aovFrequency)
#tukeyFrequency
aovDuration <- aov(duration ~ coast * event, data = allAllAnnual[allAllAnnual$type == "insitu",])
#summary(aovDuration)
tukeyDuration <- TukeyHSD(aovDuration)
#tukeyDuration
aovIntensity <- aov(intensity ~ coast * event, data = allAllAnnual[allAllAnnual$type == "insitu",])
#summary(aovIntensity)
tukeyIntensity <- TukeyHSD(aovIntensity)
#tukeyIntensity

## ANOVA for OISST
aovFrequency <- aov(frequency ~ coast * event, data = allAllAnnual[allAllAnnual$type == "OISST",])
#summary(aovFrequency)
tukeyFrequency <- TukeyHSD(aovFrequency)
#tukeyFrequency
aovDuration <- aov(duration ~ coast * event, data = allAllAnnual[allAllAnnual$type == "OISST",])
#summary(aovDuration)
tukeyDuration <- TukeyHSD(aovDuration)
#tukeyDuration
aovIntensity <- aov(intensity ~ coast * event, data = allAllAnnual[allAllAnnual$type == "OISST",])
#summary(aovIntensity)
tukeyIntensity <- TukeyHSD(aovIntensity)
#tukeyIntensity

## ANOVA for MHWs
aovFrequency <- aov(frequency ~ coast * type, data = allAllAnnual[allAllAnnual$event == "mhw",])
#summary(aovFrequency)
tukeyFrequency <- TukeyHSD(aovFrequency)
#tukeyFrequency
aovDuration <- aov(duration ~ coast * type, data = allAllAnnual[allAllAnnual$event == "mhw",])
#summary(aovDuration)
tukeyDuration <- TukeyHSD(aovDuration)
#tukeyDuration
aovIntensity <- aov(intensity ~ coast * type, data = allAllAnnual[allAllAnnual$event == "mhw",])
#summary(aovIntensity)
tukeyIntensity <- TukeyHSD(aovIntensity)
#tukeyIntensity

## ANOVA for everything
aovFrequency <- aov(frequency ~ coast * event * type, data = allAllAnnual)
#summary(aovFrequency)
tukeyFrequency <- TukeyHSD(aovFrequency)
#tukeyFrequency
aovDuration <- aov(duration ~ coast * event * type, data = allAllAnnual)
#summary(aovDuration)
tukeyDuration <- TukeyHSD(aovDuration)
#tukeyDuration
aovIntensity <- aov(intensity ~ coast * event * type, data = allAllAnnual)
#summary(aovIntensity)
tukeyIntensity <- TukeyHSD(aovIntensity)
#tukeyIntensity

### START AJS...
# I used a series of planned comparisons (general linear hypotheses) as specific contrasts. It is unnecessary to test each and every possible thing that is testable. The analyses below addresses all the comparisons in Table 2 (also the tests between columns in the table, which is not reported there).
# first for the count data; the factors are...
# type (in situ and OISST)
# event (MHW and MCS)
# coast (wc, sc and ec)
# We want to test...
# 1. is there is diff. in the number of MHWs btw the in situ and OISST data?
# 2. is there is diff. in the number of MCSs btw the in situ and OISST data?
# 3. do in situ and OISST data yield the same number of events?
# 4. within the in situ data, is there a diff in the number of MHWs and MCSs?
# 5. within the OISST data, is there a diff in the number of MHWs and MCSs?
# these are encoded by the following contrasts and analysed as general linear hypotheses:
allAllAnnual$c1 <- interaction(allAllAnnual$type, allAllAnnual$event)
levels(allAllAnnual$c1) # gives the order of the four combinations
allAllAnnual$c2 <- interaction(allAllAnnual$type, allAllAnnual$event, allAllAnnual$coast)
levels(allAllAnnual$c2)
(mod1 <- aov(frequency ~ c1, data = allAllAnnual))
summary(mod1)
summary.lm(mod1) # a more useful output
contrasts(allAllAnnual$c1)
library(multcomp) # for general linear hypotheses (glht)
cntrMat1 <- rbind("in situ-OISST (MCS)"=c(1, -1, 0, 0), # 1.
                 "in situ-OISST (MHW)"=c(0, 0, 1, -1), # 2.
                 "in situ-OISST (overall)"=c(1, -1, 1,-1), # 3.
                 "MHW-MCS (in situ)"=c(-1, 0, 1, 0), # 4.
                 "MHW-MCS (OISST)"=c(0, -1, 0, 1))# 5.
mod1.glht <- glht(mod1, linfct = mcp(c1 = cntrMat1), alternative = "two.sided")
summary(mod1.glht, test = adjusted("none"))
# 6. within the in situ data, is there a difference in MCSs between coasts (wc vs sc)?
# 7. within the in situ data, is there a difference in MCSs between coasts (wc vs ec)?
# 8. within the in situ data, is there a difference in MHWs between coasts (wc vs sc)?
# 9. within the in situ data, is there a difference in MHWs between coasts (wc vs ec)?
# 10. within the OISST data, is there a difference in MCSs between coasts (wc vs sc)?
# 11. within the OISST data, is there a difference in MCSs between coasts (wc vs ec)?
# 12. within the OISST data, is there a difference in MHWs between coasts (wc vs sc)?
# 13. within the OISST data, is there a difference in MHWs between coasts (wc vs ec)?
(mod2 <- aov(frequency ~ c2, data = allAllAnnual))
summary(mod2)
summary.lm(mod2) # a different output
cntrMat2 <- rbind("wc-sc (in situ, MCS)"=c(1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0), # 6.
                  "wc-ec (in situ, MCS)"=c(1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0), # 7.
                  "wc-sc (in situ, MHW)"=c(0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0), # 8.
                  "wc-ec (in situ, MHW)"=c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0), # 9.
                  "wc-sc (OISST, MCS)"=c(0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0), # 10.
                  "wc-ec (OISST, MCS)"=c(0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0), # 11.
                  "wc-sc (OISST, MHW)"=c(0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0), # 12.
                  "wc-ec (OISST, MHW)"=c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1), # 13.
                  "MHW-MCS (in situ, wc)"=c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 14.
                  "MHW-MCS (in situ, sc)"=c(0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0), # 15.
                  "MHW-MCS (in situ, ec)"=c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0), # 16.
                  "MHW-MCS (OISST, wc)"=c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), # 17.
                  "MHW-MCS (OISST, sc)"=c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0), # 18.
                  "MHW-MCS (OISST, ec)"=c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1))# 19.
mod2.glht <- glht(mod2, linfct = mcp(c2 = cntrMat2), alternative = "two.sided")
summary(mod2.glht, test = adjusted("none")) # a bonferroni adjustment can be applied if needed, but I specified orthogonal contrasts as far as possible.

# Below is a plain vanilla approach with unplanned ad-hoc comparisons, i.e. every combination of factor levels is compared to every other combination. Since the contrasts are still orthogonal there should not be problems with inflated type I errors, but the above is more efficient as it demonstrates more thought in the selection of hypotheses. Just for fun I have also been playing a bit with the 'magrittr' package and pipes...
mod1 <- aov(frequency ~ event * type * coast, data = allAllAnnual)
summary(mod1)
mod1.Tukey <- TukeyHSD(mod1)

mod1.Tukey$`type:event` %>% # event only
  data.frame(comp = row.names(mod1.Tukey$`type:event`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod1.Tukey$`type:event:coast` %>% # event only
  data.frame(comp = row.names(mod1.Tukey$`type:event:coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

# reveal only the significant differences:
mod1.Tukey$`event` %>% # event only
  data.frame(comp = row.names(mod1.Tukey$`event`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod1.Tukey$`type` %>% # type only
  data.frame(comp = row.names(mod1.Tukey$`type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod1.Tukey$`coast` %>% # coast only
  data.frame(comp = row.names(mod1.Tukey$`coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod1.Tukey$`event:type` %>% # OISST has more MCS and MHW
  data.frame(comp = row.names(mod1.Tukey$`event:type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod1.Tukey$`event:type:coast` %>% # comparisons of coasts within events and type
  data.frame(comp = row.names(mod1.Tukey$`event:type:coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

print(model.tables(mod1, "means"),digits = 3)
boxplot(frequency ~ event * type, data = allAllAnnual)

# now duration:
mod2 <- aov(duration ~ event * type * coast, data = allAllAnnual)
summary(mod2)
mod2.Tukey <- TukeyHSD(mod2)

# reveal only the significant differences:
mod2.Tukey$`event` %>% # event only
  data.frame(comp = row.names(mod2.Tukey$`event`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod2.Tukey$`type` %>% # type only
  data.frame(comp = row.names(mod2.Tukey$`type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod2.Tukey$`coast` %>% # coast only
  data.frame(comp = row.names(mod2.Tukey$`coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod2.Tukey$`event:type` %>% # OISST has more MCSs
  data.frame(comp = row.names(mod2.Tukey$`event:type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod2.Tukey$`event:type:coast` %>% # comparisons of coasts within events and type
  data.frame(comp = row.names(mod2.Tukey$`event:type:coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

# now duration:
mod3 <- aov(intensity ~ event * type * coast, data = allAllAnnual)
summary(mod3)
mod3.Tukey <- TukeyHSD(mod3)

# reveal only the significant differences:
mod3.Tukey$`event` %>% # event only
  data.frame(comp = row.names(mod3.Tukey$`event`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod3.Tukey$`type` %>% # type only
  data.frame(comp = row.names(mod3.Tukey$`type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod3.Tukey$`coast` %>% # coast only
  data.frame(comp = row.names(mod3.Tukey$`coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod3.Tukey$`event:type` %>% # OISST has more MCSs
  data.frame(comp = row.names(mod3.Tukey$`event:type`)) %>%
  dplyr::filter(p.adj <= 0.05)

mod3.Tukey$`event:type:coast` %>% # comparisons of coasts within events and type
  data.frame(comp = row.names(mod3.Tukey$`event:type:coast`)) %>%
  dplyr::filter(p.adj <= 0.05)

### END AJS...

## ANOVA for in situ cummulative intensity
aovIntensCum <- aov(intCum ~ coast * event, data = allAllEvent[allAllEvent$type == "insitu",])
#summary(aovIntensCum)
tukeyIntensCum <- TukeyHSD(aovIntensCum)
#tukeyIntensCum

## ANOVA for in situ cummulative intensity
aovIntensCum <- aov(intCum ~ coast * event, data = allAllEvent[allAllEvent$type == "OISST",])
#summary(aovIntensCum)
tukeyIntensCum <- TukeyHSD(aovIntensCum)
#tukeyIntensCum

## ANOVA for all cummulative intensity
aovIntensCum <- aov(intCum ~ coast * event * type, data = allAllEvent)
#summary(aovIntensCum)
tukeyIntensCum <- TukeyHSD(aovIntensCum)
#tukeyIntensCum

## Check that data are normal. These ones are not.
## First pull out your Residuals
# res <- residuals(aovFrequency )
# ## Now pull out your fitted values
# yfit <- fitted.values(aovFrequency )
# #window(width=8, height=3)
# par(mfrow=c(1,3),mex=0.8)
# ## Looking at the spread of the residuals to check Variances
# plot(yfit,res,xlab="Predicted values",ylab="Residuals")
# abline(0,0,lty=3)
# ### Looking at Normality
# qqnorm(res,main="Model Validation Graphs")
# qqline(res)
# hist(res,main= "",xlab="Residuals")

#

# GLM for everything # Not currently using this
# glmFrequency <- lm(frequency ~ coast * event * type, data = allAllAnnual)
# summary(glmFrequency)

#############################################################################
## 8. Calculate co-occurrence between coastal sections

## The code to calculate co-occurrence for individual sites may be found in "graph/eventsPlots2.R"

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

# mhwCoastCO <- ddply(mhwCoastCO, .(coast, lag, quantile, direction), summarize,
#                     proportion = mean(proportion, na.rm = TRUE))
# write.csv(mhwCoastCO, "data/mhwCoastCO.csv")

# mcsCoastCO <- ddply(mcsCoastCO, .(coast, lag, quantile, direction), summarize,
#                     proportion = mean(proportion, na.rm = TRUE))
# write.csv(mcsCoastCO, "data/mcsCoastCO.csv")

# Some exploratory stats
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="b"])
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="x"])
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="a"])

mean(mcsCoastCO$proportion[mcsCoastCO$direction =="b"])
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="x"])
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="a"])

#############################################################################
## 9. Calcuate stats and statistical significance between coastal sections etc. for co-occurrence

## Mean(sd) of co-occurrence for the different types
# MHW + all coasts + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
# MHW + all coasts + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MCS + all coasts + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
# MCS + all coasts + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])


## Mean values considering lag windows
# MHW + all coasts + before + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MHW + all coasts + before + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MCS + all coasts + before + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
# MCS + all coasts + before + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])


# MHW + all coasts + before + 50th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14])
# MHW + all coasts + before + 50th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14])
# MCS + all coasts + before + 50th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14])
# MCS + all coasts + before + 50th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14])


## Mean values considering lag length
# MHW + south coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
# MHW + south coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MCS + south coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
# MCS + south coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])

# MHW + west coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
# MHW + west coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MCS + west coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
# MCS + west coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])

# MHW + east coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2])
# MHW + east coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14])
# MCS + east coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2])
# MCS + east coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14])

## ANOVA for different coasts co-occurrence
# MHW + both lag + 2 day lag
aovCo <- aov(proportion ~ coast, data = mhwCoastCO[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2,])
#summary(aovCo)
tukeyFrequency <- TukeyHSD(aovCo)
#tukeyFrequency
# MHW + both lag + 14 day lag
aovCo <- aov(proportion ~ coast, data = mhwCoastCO[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14,])
#summary(aovCo)
tukeyFrequency <- TukeyHSD(aovCo)
#tukeyFrequency

# MCS + both lag + 2 day lag
aovCo <- aov(proportion ~ coast, data = mcsCoastCO[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2,])
#summary(aovCo)
tukeyFrequency <- TukeyHSD(aovCo)
#tukeyFrequency
# MCS + both lag + 14 day lag
aovCo <- aov(proportion ~ coast, data = mcsCoastCO[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14,])
#summary(aovCo)
tukeyFrequency <- TukeyHSD(aovCo)
#tukeyFrequency

#############################################################################
## 10. Rate of increase in MHWs/ MCSs

## Decadal trends
# First  subset out in situ time series over 30 years
is30 <- droplevels(metaData2[metaData2$length >= (30*365),])
isLong <- allAnnual[allAnnual$site %in% is30$site,]
isShort <- allAnnual[!(allAnnual$site %in% is30$site),]

# Combine all long time series
allLong <- rbind(isLong, allAnnualSST)

# Calculate trends
trends <- data.frame()
for(i in 1:length(levels(as.factor(allLong$site)))){
  dat1 <- subset(allLong, site == levels(as.factor(allLong$site))[i])
  for(j in 1:length(levels(as.factor(dat1$type)))){
    dat2 <- subset(dat1, type == levels(as.factor(dat1$type))[j])
    for(k in 1:length(levels(as.factor(dat2$event)))){
      dat3 <- subset(dat2, event == levels(as.factor(dat2$event))[k])
      lmodel <- lm(dat3$frequency ~ seq(1:length(dat3$frequency)))
      dat4 <- data.frame(site = dat3$site[1], coast = dat3$coast[1], type = dat3$type[1],
                         event = dat3$event[1], trend = round(as.numeric(coef(lmodel)[2]*10),1))
      trends <- rbind(trends, dat4)
    }
  }
}

## OISST stats
# OISST MHW all
mean(trends$trend[trends$type == "OISST" & trends$event == "mhw"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mhw"])
# OISST MHW wc
mean(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "wc"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "wc"])
# OISST MHW sc
mean(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "sc"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "sc"])
# OISST MHW ec
mean(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "ec"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mhw" & trends$coast == "ec"])

# OISST MCS
mean(trends$trend[trends$type == "OISST" & trends$event == "mcs"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mcs"])
# OISST MCS wc
mean(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "wc"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "wc"])
# OISST MCS sc
mean(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "sc"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "sc"])
# OISST MCS ec
mean(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "ec"])
sd(trends$trend[trends$type == "OISST" & trends$event == "mcs" & trends$coast == "ec"])


## in situ stats
# insitu MHW
mean(trends$trend[trends$type == "insitu" & trends$event == "mhw"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mhw"])
# insitu MHW wc
mean(trends$trend[trends$type == "insitu" & trends$event == "mhw" & trends$coast == "wc"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mhw" & trends$coast == "wc"])
# insitu MHW sc
mean(trends$trend[trends$type == "insitu" & trends$event == "mhw" & trends$coast == "sc"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mhw" & trends$coast == "sc"])

# insitu MCS
mean(trends$trend[trends$type == "insitu" & trends$event == "mcs"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mcs"])
# insitu MCS wc
mean(trends$trend[trends$type == "insitu" & trends$event == "mcs" & trends$coast == "wc"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mcs" & trends$coast == "wc"])
# insitu MCS sc
mean(trends$trend[trends$type == "insitu" & trends$event == "mcs" & trends$coast == "sc"])
sd(trends$trend[trends$type == "insitu" & trends$event == "mcs" & trends$coast == "sc"])

## Short time series
shorts <- data.frame()
for(i in 1:length(levels(as.factor(isShort$site)))){
  dat1 <- subset(isShort, site == levels(as.factor(isShort$site))[i])
  dat1.1 <- subset(metaData2, site == levels(as.factor(isShort$site))[i])
  dat1.date <- seq(year(dat1.1$start.date), year(dat1.1$end.date), by = 1)
  dat1.date.1 <- dat1.date[1:(length(dat1.date)/2)]
  dat1.date.2 <- dat1.date[(ceiling(length(dat1.date)/2)+1):length(dat1.date)]
  for(j in 1:length(levels(as.factor(dat1$event)))){
    dat2 <- subset(dat1, event == levels(as.factor(dat1$event))[j])
    dat2.1 <- dat2[dat2$year %in% dat1.date.1,]
    dat2.2 <- dat2[dat2$year %in% dat1.date.2,]
    prop <- sum(dat2.2$frequency)/sum(dat2.1$frequency)
    dat4 <- data.frame(site = dat2$site[1], coast = dat2$coast[1], type = dat2$type[1],
                       event = dat2$event[1], prop = round(prop,2))
    shorts <- rbind(shorts, dat4)
  }
}

## Short stats
# MHW all
mean(shorts$prop[shorts$event == "mhw"])
sd(shorts$prop[shorts$event == "mhw"])
# MHW wc
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "wc"])
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "wc"])
# MHW sc
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "sc"])
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "sc"])
# MHW ec
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "ec"])
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "ec"])

# MCS all
mean(shorts$prop[shorts$event == "mcs"])
sd(shorts$prop[shorts$event == "mcs"])
# MCS wc
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "wc"])
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "wc"])
# MCS sc
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "sc"])
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "sc"])
# MCS ec
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "ec"])
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "ec"])


#############################################################################
## 11. R2 between in situ and OISST time series

# Load time series for both datasets
load("prep/SA_coastal_temps.RData")
SA_coastal_temps$date <- as.Date(SA_coastal_temps$date)
load("data/OISSTdaily.Rdata")

resultsR2 <- data.frame()
for(i in 1:length(levels(SA_coastal_temps$site))) {
  x <- subset(SA_coastal_temps, site == levels(SA_coastal_temps$site)[i])
  y <- subset(OISSTdaily, site == levels(SA_coastal_temps$site)[i])
  x <- x[x$date %in% y$date,]
  y <- y[y$date %in% x$date,]
  z <- data.frame(site = as.character(x$site[1]), R2 = round(coef(lm(y$temp~x$temp))[2],2))
  resultsR2 <- rbind(resultsR2, z)
}

# Add column for plotting
row.names(resultsR2) <- NULL
resultsR2$R22 <- paste0("R^2 ==", format(resultsR2$R2, digits=2))

#############################################################################
## 12. Co-occurrence within datasets and coastal sections



