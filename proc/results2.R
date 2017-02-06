#############################################################################
## This script does:
# 1. Load annual MHW/ MCS files and clean meta-data for publication;
# 2. Calculates event count, duration and mean intensity for each coast for both datasets;
# 3. Load all events and extract top three MHWs/ MCSs;
# 4. Extracts top three MHWs/ MCSs for each coastal section and type of data;
# 5. Extracts top 1 MHW/ MCS for each coastal section and type of data;
# 6. Calcuate statistical significance between coastal sections;
# 7. Calculate co-occurrence values;
# 8. Calcuate stats and statistical significance between coastal sections etc. for co-occurrence;
# 9. Decadal trends in MHWs/ MCSs
# 10. R2 between in situ and OISST time series
#############################################################################

#############################################################################
## DEPENDS ON:
# setwd("/Users/ajsmit/Dropbox/̧repos/MHW/proc")
library(zoo)
library(plyr)
library(dplyr)
library(stringr)
library(lubridate)
library(xtable)
library(magrittr)
library(multcomp)
library(doMC); doMC::registerDoMC(4)
source("setupParams/theme.R")
source("func/seqSites.R")
# "graph/eventsPlots2.R" # This script calculates the co-occurrence rates for sites
# "data/metaData2.csv"
# Several files from Eric CJ Oliver generated via Python
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


# 1. Load annual and event MHW/ MCS files and metadata2 -------------------

# First specify coastal sections
wc <- c("Hout Bay", "Kommetjie", "Port Nolloth", "Sea Point")
sc <- c("Fish Hoek", "Gordons Bay", "Hamburg", "Hermanus", "Humewood", "Knysna",
        "Mossel Bay", "Muizenberg", "Pollock Beach", "Tsitsikamma West",
        "Storms River Mouth", "Tsitsikamma East", "Ystervarkpunt")
ec <- c("Eastern Beach", "Nahoon Beach", "Orient Beach", "Sodwana")

# Then specify directories for loading
# Annual values
dir1 <- paste(getwd(), "/data/MHW/annual", sep = "")
dir2 <- paste(getwd(), "/data/MCS/annual", sep = "")
dir3 <- paste(getwd(), "/data/MHW/SST annual", sep = "")
dir4 <- paste(getwd(), "/data/MCS/SST annual", sep = "")
# Individual events
dir5 <- paste(getwd(), "/data/MHW/events", sep = "")
dir6 <- paste(getwd(), "/data/MCS/events", sep = "")
dir7 <- paste(getwd(), "/data/MHW/SST events", sep = "")
dir8 <- paste(getwd(), "/data/MCS/SST events", sep = "")

# Function that loads individual annual files
annual.load <- function(file) {
  x <- as.character(file)
  y <- read.csv(x, header = TRUE, skip = 2)
  y$site <- sapply(strsplit(as.character(x), "/"), "[[", 8)
  y$site <- sapply(strsplit(y$site, "_MHW_data.annual.csv"), "[[", 1)
  y$site <- sapply(strsplit(y$site, "_MCS_data.annual.csv"), "[[", 1)
  y$site <- str_replace_all(y$site, "_", " ")
  if(y$site[1] %in% wc) {
    y$coast <- "wc"
  } else if(y$site[1] %in% sc) {
    y$coast <- "sc"
  } else if(y$site[1] %in% ec) {
    y$coast <- "ec"
  } else {
    stop(paste(y$site, " not detected in any coastal region.", sep = ""))
  }
  y <- y[,c(25:26,1,2,3,5)]
  colnames(y)[3:6] <- c("year","frequency","duration","intensity")
  y[is.na(y$frequency)] <- 0
  return(y)
}

# Function that loads all annual files in a folder
annual.load.all <- function(directory) {
  fname <- dir(directory, full.names = TRUE)
  dat <- ldply(fname, annual.load, .parallel = T)
  dat$coast <- factor(dat$coast, levels = c("wc", "sc", "ec"))
  return(dat)
}

# Load all annual files
mhwAnnual <- annual.load.all(dir1)
mcsAnnual <- annual.load.all(dir2)
mhwAnnualSST <- annual.load.all(dir3)
mcsAnnualSST <- annual.load.all(dir4)

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
metaData3$length <- round(metaData3$length/365,1) # Convert to year
colnames(metaData3)[c(6:9)] <- c("start date", "end date", "duration (years)", "NA %")
xtable(metaData3, auto = TRUE)


# 2. Calculates event count, duration and mean intensity for each  --------

# Function used for calculations
results.annual.coastal <- function(mhw1, mcs1){ # To be used with "annual" data frames only
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

allAnnualResults <- results.annual.coastal(mhwAnnual, mcsAnnual)
xtable(allAnnualResults)
write.csv(allAnnualResults, "data/allAnnualResults.csv")
allAnnualSSTResults <- results.annual.coastal(mhwAnnualSST, mcsAnnualSST)
xtable(allAnnualSSTResults)
write.csv(allAnnualSSTResults, "data/allAnnualSSTResults.csv")


# 3. Load all events and extract top three MHWs/ MCSs ---------------------

# Pretty column names
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

# Function that loads individual event files
event.load <- function(file) {
  x <- as.character(file)
  y <- read.csv(x, header = FALSE, skip = 2, sep = ",",
                col.names = colNames)
  y$site <- sapply(strsplit(as.character(x), "/"), "[[", 8)
  y$site <- sapply(strsplit(y$site, "_MHW_data.events.csv"), "[[", 1)
  y$site <- sapply(strsplit(y$site, "_MCS_data.events.csv"), "[[", 1)
  y$site <- str_replace_all(y$site, "_", " ")
  return(y)
}

# Function that loads all event files in a folder
event.load.all <- function(directory) {
  fname <- dir(directory, full.names = TRUE)
  dat <- ldply(fname, event.load, .parallel = T)
  dat <- dat %>%
    group_by(site) %>%
    mutate(coast = metaData2$coast[metaData2$site == site[1]])
  dat$date <-  as.Date(paste(dat$yearStrt, dat$monthStrt, dat$dayStrt, sep = "-"))
  dat$month <- floor_date(as.Date(paste(dat$yearStrt, dat$monthStrt, dat$dayStrt, sep = "-")), "month")
  dat <- dat %>%
    group_by(site) %>%
    mutate(lon = metaData2$lon[metaData2$site == site[1]]) %>% 
    mutate(lat = metaData2$lat[metaData2$site == site[1]])
  return(dat)
}

# Load the files
mhwEvent <- event.load.all(dir5)
mcsEvent <- event.load.all(dir6)
mhwEventSST <- event.load.all(dir7)
mcsEventSST <- event.load.all(dir8)

# The top three events per site
event.load.n <- function(dir, nCum = 5) {
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
mhwn <- event.load.n(dir5, nCum = 3) # produce a data frame with mhw data...
mhwn$index <- rep(1:3, length(levels(as.factor(mhwn$site)))) # Make sure to correct the rep() to match nCum
mhwn$type <- "insitu"
mcsn <- event.load.n(dir6, nCum = 3) # produce a data frame with mcs data...
mcsn$index <- rep(1:3, length(levels(as.factor(mcsn$site))))
mcsn$type <- "insitu"
# SST
mhwnSST <- event.load.n(dir7, nCum = 3)
mhwnSST$index <- rep(1:3, length(levels(as.factor(mhwnSST$site))))
mhwnSST$type <- "OISST"
mcsnSST <- event.load.n(dir8, nCum = 3)
mcsnSST$index <- rep(1:3, length(levels(as.factor(mcsnSST$site))))
mcsnSST$type <- "OISST"


# 4. Extracts top three MHW/ MCS events for each coastal section a --------

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


# 5. Extracts top 1 MHW/ MCS events for each coastal section and t --------

mhw1 <- topCoastn(mhwn, 1)
mcs1 <- topCoastn(mcsn, 1)
mhwSST1 <- topCoastn(mhwnSST, 1)
mcsSST1 <- topCoastn(mcsnSST, 1)


# 6. Calcuate statistical significance between coastal sections -----------

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
mean(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mhw"], na.rm = T) # 26.10681
sd(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mhw"], na.rm = T) # 24.36918
# in situ + MCS
mean(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mcs"], na.rm = T) # 26.44975
sd(abs(allAllEvent$intCum)[allAllEvent$type == "insitu" & allAllEvent$event == "mcs"], na.rm = T) # 24.24552
# OISST + MHW
mean(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mhw"], na.rm = T) # 18.64982
sd(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mhw"], na.rm = T) # 15.10139
# OISST + MCS
mean(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mcs"], na.rm = T) # 23.16813
sd(abs(allAllEvent$intCum)[allAllEvent$type == "OISST" & allAllEvent$event == "mcs"], na.rm = T) # 23.4894

### START AJS...
# I used a series of planned comparisons (generalized linear hypothesis test) as specific contrasts. 
# It is unnecessary to test each and every possible thing that is testable. 
# The analyses below address all the comparisons in Table 2 (also the tests between columns in the table, which is not reported there).
# First for the count data; the factors are...
# type (in situ and OISST)
# event (MHW and MCS)
# coast (wc, sc and ec)
# We want to test...
# 1. is there a diff. in the number of MHWs btw the in situ and OISST data?
# 2. is there a diff. in the number of MCSs btw the in situ and OISST data?
# 3. do in situ and OISST data yield the same number of events?
# 4. within the in situ data, is there a diff in the number of MHWs and MCSs?
# 5. within the OISST data, is there a diff in the number of MHWs and MCSs?
# these are encoded by the following contrasts and analysed as general linear hypotheses:
allAllAnnual$c1 <- interaction(allAllAnnual$type, allAllAnnual$event)
levels(allAllAnnual$c1) # gives the order of the four combinations; used for specifying contrasts
allAllAnnual$c2 <- interaction(allAllAnnual$type, allAllAnnual$event, allAllAnnual$coast)
levels(allAllAnnual$c2)
(mod1 <- aov(frequency ~ c1, data = allAllAnnual))
print(model.tables(mod1, "means"),digits = 3)
boxplot(frequency ~ c1, data = allAllAnnual)
summary(mod1)
summary.lm(mod1) # a more useful output
contrasts(allAllAnnual$c1)
library(multcomp) # for generalized linear hypothesis test (glht)
cntrMat1 <- rbind("in situ-OISST (MCS)"=c(1, -1, 0, 0), # 1.
                 "in situ-OISST (MHW)"=c(0, 0, 1, -1), # 2.
                 "in situ-OISST (overall)"=c(1, -1, 1,-1), # 3.
                 "MHW-MCS (in situ)"=c(-1, 0, 1, 0), # 4.
                 "MHW-MCS (OISST)"=c(0, -1, 0, 1))# 5.
mod1.glht <- glht(mod1, linfct = mcp(c1 = cntrMat1), alternative = "two.sided")
library(sandwich)
# with sandwich estimator to account for the small differences in heterogeneity between the c1 factors (see plot above)
mod1a.glht <- glht(mod1, linfct = mcp(c1 = cntrMat1), alternative = "two.sided", vcov = sandwich)
summary(mod1a.glht)
summary(mod1a.glht, test = adjusted("none"))
coef(mod1a.glht)
confint(mod1a.glht) # confidence intervals of the difference between the contrasts
plot(mod1.glht)
plot(mod1a.glht)
## 6. within the in situ data, is there a difference in MCSs between coasts (wc vs sc)?
## 7. within the in situ data, is there a difference in MCSs between coasts (wc vs ec)?
## 8. within the in situ data, is there a difference in MCSs between coasts (sc vs ec)?
# 9. within the in situ data, is there a difference in MHWs between coasts (wc vs sc)?
# 10. within the in situ data, is there a difference in MHWs between coasts (wc vs ec)?
# 11. within the in situ data, is there a difference in MHWs between coasts (sc vs ec)?
## 12. within the OISST data, is there a difference in MCSs between coasts (wc vs sc)?
## 13. within the OISST data, is there a difference in MCSs between coasts (wc vs ec)?
## 14. within the OISST data, is there a difference in MCSs between coasts (sc vs ec)?
# 15. within the OISST data, is there a difference in MHWs between coasts (wc vs sc)?
# 16. within the OISST data, is there a difference in MHWs between coasts (wc vs ec)?
# 17. within the OISST data, is there a difference in MHWs between coasts (sc vs ec)?
## 18. within in situ and wc, are the number of cold and warm events the same?
## 19. within in situ and sc, are the number of cold and warm events the same?
## 20. within in situ and ec, are the number of cold and warm events the same?
# 21. within OISST and wc, are the number of cold and warm events the same?
# 22. within OISST and sc, are the number of cold and warm events the same?
# 23. within OISST and ec, are the number of cold and warm events the same?
(mod2 <- aov(frequency ~ c2, data = allAllAnnual))
summary(mod2)
summary.lm(mod2) # a different output
cntrMat2 <- rbind("wc-sc (in situ, MCS)"=c(1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0), # 6.
                  "wc-ec (in situ, MCS)"=c(1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0), # 7.
                  "sc-ec (in situ, MCS)"=c(0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0), # 8.

                  "wc-sc (in situ, MHW)"=c(0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0), # 9.
                  "wc-ec (in situ, MHW)"=c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0), # 10.
                  "sc-ec (in situ, MHW)"=c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0), # 11.

                  "wc-sc (OISST, MCS)"=c(0, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0), # 12.
                  "wc-ec (OISST, MCS)"=c(0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0), # 13.
                  "sc-ec (OISST, MCS)"=c(0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0), # 14.

                  "wc-sc (OISST, MHW)"=c(0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0), # 15.
                  "wc-ec (OISST, MHW)"=c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, -1), # 16.
                  "sc-ec (OISST, MHW)"=c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1), # 17.

                  "MHW-MCS (in situ, wc)"=c(-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0), # 18.
                  "MHW-MCS (in situ, sc)"=c(0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0), # 19.
                  "MHW-MCS (in situ, ec)"=c(0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0), # 20.

                  "MHW-MCS (OISST, wc)"=c(0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0), # 21.
                  "MHW-MCS (OISST, sc)"=c(0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0), # 22.
                  "MHW-MCS (OISST, ec)"=c(0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1))# 23.
mod2.glht <- glht(mod2, linfct = mcp(c2 = cntrMat2), alternative = "two.sided", vcov = sandwich)
summary(mod2.glht, test = adjusted("none")) # a bonferroni adjustment can be applied if needed, but I specified orthogonal contrasts as far as possible.
coef(mod2.glht)
confint(mod2.glht)

# now for duration:
(mod3 <- aov(duration ~ c1, data = allAllAnnual))
mod3a.glht <- glht(mod3, linfct = mcp(c1 = cntrMat1), alternative = "two.sided", vcov = sandwich)
summary(mod3a.glht, test = adjusted("none"))

(mod4 <- aov(duration ~ c2, data = allAllAnnual))
mod4a.glht <- glht(mod4, linfct = mcp(c2 = cntrMat2), alternative = "two.sided", vcov = sandwich)
summary(mod4a.glht, test = adjusted("none"))

# ... and intensity:
(mod5 <- aov(abs(intensity) ~ c1, data = allAllAnnual))
boxplot(abs(intensity) ~ c1, data = allAllAnnual)
mod5a.glht <- glht(mod5, linfct = mcp(c1 = cntrMat1), alternative = "two.sided", vcov = sandwich)
summary(mod5a.glht, test = adjusted("none")) # --> Note that it is important to use abs(intensity)

(mod6 <- aov(abs(intensity) ~ c2, data = allAllAnnual))
mod6a.glht <- glht(mod6, linfct = mcp(c2 = cntrMat2), alternative = "two.sided", vcov = sandwich)
summary(mod6a.glht, test = adjusted("none"))
### END AJS...


# 7. Calculate co-occurrence values ---------------------------------------

# Function used for calculating co-occurrence
cooccurrence <- function(dat1, dat2, lag = seq(2,14,2)){
  dat3 <- data.frame()
  direction <- c("b","x","a")
  for(i in 1:length(levels(as.factor(dat1$site)))) {
    x1 <- droplevels(subset(dat1, site == levels(as.factor(dat1$site))[i]))
    x2 <- droplevels(subset(dat2, site == levels(as.factor(dat1$site))[i]))
    x1 <- x1[x1$yearStrt >= min(x2$yearStrt), ] # Subset x so that dates match up
    x1 <- x1[x1$yearStrt <= max(x2$yearStrt), ]
    x2 <- x2[x2$yearStrt >= min(x1$yearStrt), ]
    x2 <- x2[x2$yearStrt <= max(x1$yearStrt), ]
    for(j in 1:length(lag)){
      for(k in 1:length(seq(0.0,1,0.1))){
        for(l in 1:length(direction)){
          x1.1 <- x1[x1$intCum >= quantile(x1$intCum, probs = seq(0.0,1,0.1)[k]),]
          x2.1 <- x2[x2$intCum >= quantile(x2$intCum, probs = seq(0.0,1,0.1)[k]),]
          y <- 0
          #x3 <- data.frame() # For test purposes to see which events match up
          for(m in 1:nrow(x1.1)) {
            x1.2 <- x1.1$date[m]
            if(direction[l] == "b"){
              x1.3 <- seq((x1.2 - days(lag[j])), x1.2, 1)
            } else if(direction[l] == "x"){
              x1.3 <- seq((x1.2 - days(lag[j])), (x1.2 + days(lag[j])), 1)
            } else if (direction[l] == "a") {
              x1.3 <- seq(x1.2, (x1.2 + days(lag[j])), 1)
            }
            x2.2 <- droplevels(subset(x2.1, date %in% x1.3))
            y <- y + nrow(x2.2)
          }
          z <- data.frame(site = x1$site[1], lon = x1$lon[1], lat = x1$lat[1],
                          lag = lag[j], quantile = seq(0.0,1,0.1)[k], direction = direction[l],
                          insitu = nrow(x1.1), OISST = nrow(x2.1),
                          cooccurrence= y, proportion = y/nrow(x1.1))
          dat3 <- rbind(dat3, z)
        }
      }
    }
  }
  return(dat3)
}

#mhwCO0 <- cooccurrence(mhw, mhwSST, lag = 0) # Test to see which happen on exact same day
#MHW
# mhwCoastCO <- cooccurrence(mhwEvent, mhwEventSST) # This takes several minutes to run... cursed for loops...
# write.csv(mhwCoastCO, "data/mhwCO.csv", row.names = F)
# mhwCoastCO <- read.csv("data/mhwCO.csv")
#MCS
# mcsCoastCO <- cooccurrence(mcsEvent, mcsEventSST)
# write.csv(mcsCoastCO, "data/mcsCO.csv", row.names = F)
# mcsCoastCO <- read.csv("data/mcsCO.csv")

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
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="b"]) # 0.1067157
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="x"]) # 0.2111034
mean(mhwCoastCO$proportion[mhwCoastCO$direction =="a"]) # 0.1171394

mean(mcsCoastCO$proportion[mcsCoastCO$direction =="b"]) # 0.04615112
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="x"]) # 0.0938059
mean(mcsCoastCO$proportion[mcsCoastCO$direction =="a"]) # 0.05763613


# 8. Calcuate stats and statistical significance between coastal s --------

## Mean(sd) of co-occurrence for the different types
# MHW + all coasts + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0894
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0683
# MHW + all coasts + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.3842
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.2004
# MCS + all coasts + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.1007
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0542
# MCS + all coasts + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.3030
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.1264


## Mean values considering lag windows
# MHW + all coasts + before + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.2216
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.1315
# MHW + all coasts + before + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.1762
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.0978
# MCS + all coasts + before + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.1628
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.0958
# MCS + all coasts + before + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.1682
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.0816


# MHW + all coasts + before + 50th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14]) # 0.1693
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "b" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14]) # 0.1126
# MHW + all coasts + before + 50th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14]) # 0.1514
sd(mhwCoastCO$proportion[mhwCoastCO$direction == "a" & mhwCoastCO$quantile == 0.5 & mhwCoastCO$lag == 14]) # 0.1271
# MCS + all coasts + before + 50th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14]) # 0.0503
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "b" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14]) # 0.0773
# MCS + all coasts + before + 50th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14]) # 0.0807
sd(mcsCoastCO$proportion[mcsCoastCO$direction == "a" & mcsCoastCO$quantile == 0.5 & mcsCoastCO$lag == 14]) # 0.0761


## Mean values considering lag length
# MHW + south coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.1040
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0709
# MHW + south coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.4540
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "south" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.1754
# MCS + south coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.1117
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0635
# MCS + south coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.3391
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "south" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.1421

# MHW + west coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0719
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0334
# MHW + west coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.2834
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "west" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.0567
# MCS + west coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0824
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0352
# MCS + west coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.1920
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "west" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.0223

# MHW + east coast + both directions + 0th quantile + 2 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0593
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 2]) # 0.0857
# MHW + east coast + both directions + 0th quantile + 14 day lag
mean(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.2580
sd(mhwCoastCO$proportion[mhwCoastCO$coast == "east" & mhwCoastCO$direction == "x" & mhwCoastCO$quantile == 0.0 & mhwCoastCO$lag == 14]) # 0.2926
# MCS + east coast + both directions + 0th quantile + 2 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0836
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 2]) # 0.0296
# MCS + east coast + both directions + 0th quantile + 14 day lag
mean(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.2967
sd(mcsCoastCO$proportion[mcsCoastCO$coast == "east" & mcsCoastCO$direction == "x" & mcsCoastCO$quantile == 0.0 & mcsCoastCO$lag == 14]) # 0.0556

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


# 9. Decadal trends in MHWs/ MCSs -----------------------------------------

## Decadal trends
# First  subset out in situ time series over 30 years
is30 <- droplevels(metaData2[metaData2$length >= (30*365),])
isLong <- allAnnual[allAnnual$site %in% is30$site,]
isShort <- allAnnual[!(allAnnual$site %in% is30$site),]

# Combine all long time series
allLong <- rbind(isLong, allAnnualSST)

# Calculate trends using a GLM (Poisson with log-link)
trends <- data.frame()
for(i in 1:length(levels(as.factor(allLong$site)))){
  dat1 <- subset(allLong, site == levels(as.factor(allLong$site))[i])
  for(j in 1:length(levels(as.factor(dat1$type)))){
    dat2 <- subset(dat1, type == levels(as.factor(dat1$type))[j])
    mhwlmodel <- glm(dat2$frequency[dat2$event == "mhw"] ~ seq(1:length(dat2$frequency[dat2$event == "mhw"])), 
                     family = poisson(link = "log"))
    mhwlmodel0 <- glm(dat2$frequency[dat2$event == "mhw"] ~ 1, family = poisson(link = "log")) # intercept only
    mcslmodel <- glm(dat2$frequency[dat2$event == "mcs"] ~ seq(1:length(dat2$frequency[dat2$event == "mcs"])), 
                     family = poisson(link = "log"))
    mcslmodel0 <- glm(dat2$frequency[dat2$event == "mcs"] ~ 1, family = poisson(link = "log")) # intercept only
    dat3 <- data.frame(site = dat2$site[1], coast = dat2$coast[1], type = dat2$type[1],
                       mhwtrend = round(as.numeric(coef(mhwlmodel)[2]*10),1),
                       mhwpR2 = round(1-logLik(mhwlmodel)/logLik(mhwlmodel0),2), # McFadden's pseudo-R2
                       mhwp.val = round(coef(summary(mhwlmodel))[,4][2],2),
                       mcstrend = round(as.numeric(coef(mcslmodel)[2]*10),1),
                       mcspR2 = round(1-logLik(mcslmodel)/logLik(mcslmodel0),2), # McFadden's pseudo-R2
                       mcsp.val = round(coef(summary(mcslmodel))[,4][2],2))
    trends <- rbind(trends, dat3)
  }
  row.names(trends) <- NULL
}


## Order data.frame for use in the paper
trendsTable <- trends[order(trends$type, trends$coast),]
rownames(trendsTable) <- NULL
trendsTable <- trendsTable[c(22:25,1:21),]
rownames(trendsTable) <- NULL
trendsTable <- trendsTable[c(1,2,4,3,7,8,5,6,9,16,10,12,21,15,14,20,18,19,17,13,11,22,24,23,25),]
trendsTable$ID <- as.character(c(1,2,6,7,1:21))
rownames(trendsTable) <- NULL
trendsTable <- trendsTable[,c(10,2,4:9)]
trendsTable$coast <- as.character(trendsTable$coast)
trendsTable$coast[trendsTable$coast == "wc"] <- "west"
trendsTable$coast[trendsTable$coast == "sc"] <- "south"
trendsTable$coast[trendsTable$coast == "ec"] <- "east"
print(xtable(trendsTable[,1:8], auto = TRUE), include.rownames=FALSE)

## OISST stats
# OISST MHW all
mean(trends$mhwtrend[trends$type == "OISST"]) # 0.22
sd(trends$mhwtrend[trends$type == "OISST"]) # 0.15
# OISST MHW wc
mean(trends$mhwtrend[trends$type == "OISST" & trends$coast == "wc"]) # 0.13
sd(trends$mhwtrend[trends$type == "OISST" & trends$coast == "wc"]) # 0.13
# OISST MHW sc
mean(trends$mhwtrend[trends$type == "OISST" & trends$coast == "sc"]) # 0.24
sd(trends$mhwtrend[trends$type == "OISST" & trends$coast == "sc"]) # 0.17
# OISST MHW ec
mean(trends$mhwtrend[trends$type == "OISST" & trends$coast == "ec"]) # 0.25
sd(trends$mhwtrend[trends$type == "OISST" & trends$coast == "ec"]) # 0.06

# OISST MCS all
mean(trends$mcstrend[trends$type == "OISST"]) # -0.34
sd(trends$mcstrend[trends$type == "OISST"]) # 0.27
# OISST MCS wc
mean(trends$mcstrend[trends$type == "OISST" & trends$coast == "wc"]) # -0.03
sd(trends$mcstrend[trends$type == "OISST" & trends$coast == "wc"]) # 0.21
# OISST MCS sc
mean(trends$mcstrend[trends$type == "OISST" & trends$coast == "sc"]) # -0.41
sd(trends$mcstrend[trends$type == "OISST" & trends$coast == "sc"]) # 0.25
# OISST MCS ec
mean(trends$mcstrend[trends$type == "OISST" & trends$coast == "ec"]) # -0.45
sd(trends$mcstrend[trends$type == "OISST" & trends$coast == "ec"]) # 0.17


## in situ stats
# in situ MHW all
mean(trends$mhwtrend[trends$type == "insitu"]) # 0.20
sd(trends$mhwtrend[trends$type == "insitu"]) # 0.37
# in situ MHW wc
mean(trends$mhwtrend[trends$type == "insitu" & trends$coast == "wc"]) # 0.00
sd(trends$mhwtrend[trends$type == "insitu" & trends$coast == "wc"]) # 0.28
# in situ MHW sc
mean(trends$mhwtrend[trends$type == "insitu" & trends$coast == "sc"]) # 0.40
sd(trends$mhwtrend[trends$type == "insitu" & trends$coast == "sc"]) # 0.42


# in situ MCS all
mean(trends$mcstrend[trends$type == "insitu"]) # -0.08
sd(trends$mcstrend[trends$type == "insitu"]) # 0.34
# in situ MCS wc
mean(trends$mcstrend[trends$type == "insitu" & trends$coast == "wc"]) # -0.15
sd(trends$mcstrend[trends$type == "insitu" & trends$coast == "wc"]) # 0.49
# in situ MCS sc
mean(trends$mcstrend[trends$type == "insitu" & trends$coast == "sc"]) # 0.00
sd(trends$mcstrend[trends$type == "insitu" & trends$coast == "sc"]) # 0.29

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
    dat3 <- data.frame(site = dat2$site[1], coast = dat2$coast[1], type = dat2$type[1],
                       event = dat2$event[1], prop = round(prop,2))
    shorts <- rbind(shorts, dat3)
  }
}

## ANOVA for event count proportions for short time series
aovShorts <- aov(prop ~ event * coast, data = shorts)
#summary(aovShorts)
tukeyIntensCum <- TukeyHSD(aovShorts)
#tukeyIntensCum

## Short stats
# MHW all
mean(shorts$prop[shorts$event == "mhw"]) # 1.72
sd(shorts$prop[shorts$event == "mhw"]) # 1.31
# MHW wc
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "wc"]) # 1.50
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "wc"]) # 0.60
# MHW sc
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "sc"]) # 2.14
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "sc"]) # 1.42
# MHW ec
mean(shorts$prop[shorts$event == "mhw" & shorts$coast == "ec"]) # 0.68
sd(shorts$prop[shorts$event == "mhw" & shorts$coast == "ec"]) # 0.36

# MCS all
mean(shorts$prop[shorts$event == "mcs"]) # 0.79
sd(shorts$prop[shorts$event == "mcs"]) # 0.65
# MCS wc
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "wc"]) # 1.81
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "wc"]) # 0.58
# MCS sc
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "sc"]) # 0.51
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "sc"]) # 0.34
# MCS ec
mean(shorts$prop[shorts$event == "mcs" & shorts$coast == "ec"]) # 1.05
sd(shorts$prop[shorts$event == "mcs" & shorts$coast == "ec"]) # 0.82

## Short counts for ANOVA
shortCount <- data.frame()
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
    dat3 <- data.frame(site = dat2$site[1], coast = dat2$coast[1], type = dat2$type[1],
                       event = dat2$event[1], count = sum(dat2.1$frequency), half = "first")
    dat4 <- data.frame(site = dat2$site[1], coast = dat2$coast[1], type = dat2$type[1],
                       event = dat2$event[1], count = sum(dat2.2$frequency), half = "second")
    shortCount <- rbind(shortCount, dat3, dat4)
  }
}

shortCount$c1 <- interaction(shortCount$coast, shortCount$event, shortCount$half)
levels(shortCount$c1)
(mod7 <- aov(count ~ c1, data = shortCount))
print(model.tables(mod7, "means"),digits = 3)
boxplot(count ~ c1, data = shortCount)
summary(mod7)
summary.lm(mod7)
contrasts(shortCount$c1)
# the factors are...
# coast (wc, sc and ec)
# event (MHW and MCS)
# half (first or second)
# We want to test...
# 1. is there a diff. in the number of MHWs btw the first and second halfs of the data?
# 2. is there a diff. in the number of MCSs btw the first and second halfs of the data?
# 3. within the first half of the data, is there a diff in the number of MHWs and MCSs?
# 4. within the second half of the data, is there a diff in the number of MHWs and MCSs?
# 5. within the wc, is there a difference in MHWs between halves?
# 6. within the wc, is there a difference in MCSs between halves?
# 7. within the sc, is there a difference in MHWs between halves?
# 8. within the sc, is there a difference in MCSs between halves?
# 9. within the ec, is there a difference in MHWs between halves?
# 10. within the ec, is there a difference in MCSs between halves?
cntrMat7 <- rbind("first-second (MHW)"=c(0, 0, 0, 1, 1, 1, 0, 0, 0, -1, -1, -1), # 1.
                  "first-second (MCS)"=c(1, 1, 1, 0, 0, 0, -1, -1, -1, 0, 0, 0), # 2.
                  
                  "MHW-MCS (first half)"=c(-1, -1, -1, 1, 1, 1, 0, 0, 0, 0, 0, 0), # 3.
                  "MHW-MCS (secnd half)"=c(0, 0, 0, 0, 0, 0, -1, -1, -1, 1, 1, 1), # 4.
                  
                  "first-second (MHW, wc)"=c(0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0), # 5.
                  "first-second (MCS, wc)"=c(1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0), # 6.
                  
                  "first-second (MHW, sc)"=c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0), # 7.
                  "first-second (MCS, sc)"=c(0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0), # 8.
                  
                  "first-second (MHW, ec)"=c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1), # 9.
                  "first-second (MCS, ec)"=c(0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0)) # 10.
mod7.glht <- glht(mod7, linfct = mcp(c1 = cntrMat7), alternative = "two.sided", vcov = sandwich)
summary(mod7.glht, test = adjusted("none"))


# 10. R2 between in situ and OISST time series ----------------------------

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