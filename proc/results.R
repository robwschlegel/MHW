## Extract results of heat wave numbers etc by coast and by time during time series
require(zoo)

#############################################################################
# Need two different load loops because .csv files have different lines to skip
test.annual <- read.csv("data/MHW/annual/Eastern_Beach_MHW_data.annual.csv", skip = 2)
test.events <- read.csv("data/MHW/events/Eastern_Beach_MHW_data.events.csv", skip = 1)

#############################################################################
## Load all meta-data
# Annual MHW
files_list <- dir("data/MHW/annual", full.names = FALSE)
annualMHW <- data.frame()
for(i in 1:length(files_list)) {
  x <- read.csv(paste(getwd(), "/data/MHW/annual/", files_list[i], sep = ""), header = TRUE, skip = 2, sep = ",")
  fileBase1 <- strsplit(files_list[i], "_MHW") # Split character vectors wherever there is a "_"
  fileBase1.1 <- sapply(fileBase1, "[[", 1) # Selects the first split, site name
  x$site <- fileBase1.1
  x <- x[, c(25, 1:24)]
  annualMHW <- rbind(annualMHW, x)
}
# Events MHW
files_list <- dir("data/MHW/events/", full.names = FALSE)
eventMHW <- data.frame()
for(i in 1:length(files_list)) {
  x <- read.csv(paste(getwd(), "/data/MHW/events/", files_list[i], sep = ""), header = TRUE, skip = 1, sep = ",")
  fileBase1 <- strsplit(files_list[i], "_MHW") # Split character vectors wherever there is a "_"
  fileBase1.1 <- sapply(fileBase1, "[[", 1) # Selects the first split, site name
  x$site <- fileBase1.1
  x <- x[, c(28, 1:27)]
  eventMHW <- rbind(eventMHW, x)
}
# Annual MCS
files_list <- dir("data/MCS/annual", full.names = FALSE)
annualMCS <- data.frame()
for(i in 1:length(files_list)) {
  x <- read.csv(paste(getwd(), "/data/MCS/annual/", files_list[i], sep = ""), header = TRUE, skip = 2, sep = ",")
  fileBase1 <- strsplit(files_list[i], "_MCS") # Split character vectors wherever there is a "_"
  fileBase1.1 <- sapply(fileBase1, "[[", 1) # Selects the first split, site name
  x$site <- fileBase1.1
  x <- x[, c(25, 1:24)]
  annualMCS <- rbind(annualMCS, x)
}
# Events MCS
files_list <- dir("data/MCS/events/", full.names = FALSE)
eventMCS <- data.frame()
for(i in 1:length(files_list)) {
  x <- read.csv(paste(getwd(), "/data/MCS/events/", files_list[i], sep = ""), header = TRUE, skip = 1, sep = ",")
  fileBase1 <- strsplit(files_list[i], "_MCS") # Split character vectors wherever there is a "_"
  fileBase1.1 <- sapply(fileBase1, "[[", 1) # Selects the first split, site name
  x$site <- fileBase1.1
  x <- x[, c(28, 1:27)]
  eventMCS <- rbind(eventMCS, x)
}

#############################################################################
# Note broad patterns along entire coast
  # This is done visually by looking at Eric's maps

#############################################################################
## Create dataframes for time series by coastal sections
wc <- c("Hout_Bay", "Kommetjie", "Port_Nolloth", "Sea_Point")
sc <- c("Fish_Hoek", "Gordons_Bay", "Hamburg", "Hermanus", "Humewood", "Knysna", "Mossel_Bay", "Muizenberg", "Pollock_Beach", "Storms_River_Mouth", "Tsitsikamma", "Ystervarkpunt")
ec <- c("Eastern_Beach", "Nahoon_Beach", "Orient_Beach", "Sodwana")

# Annual MHW
annualMHWwc <- droplevels(subset(annualMHW, site %in% wc))
annualMHWsc <- droplevels(subset(annualMHW, site %in% sc))
annualMHWec <- droplevels(subset(annualMHW, site %in% ec))
# Event MHW
eventMHWwc <- droplevels(subset(eventMHW, site %in% wc))
eventMHWsc <- droplevels(subset(eventMHW, site %in% sc))
eventMHWec <- droplevels(subset(eventMHW, site %in% ec))
# Annual MCS
annualMCSwc <- droplevels(subset(annualMCS, site %in% wc))
annualMCSsc <- droplevels(subset(annualMCS, site %in% sc))
annualMCSec <- droplevels(subset(annualMCS, site %in% ec))
# Event MCS
eventMCSwc <- droplevels(subset(eventMCS, site %in% wc))
eventMCSsc <- droplevels(subset(eventMCS, site %in% sc))
eventMCSec <- droplevels(subset(eventMCS, site %in% ec))

#############################################################################
## Analyse relevant statistiics per coast
# Create functions for cleaner calculating... better to not use for loops...
results <- function(x){ # To be used with "annual" data frames only
  event_counts <- data.frame()
  event_lengths <- data.frame()
  event_cum_ins <- data.frame()
  trend_counts <- data.frame()
  trend_lengths <- data.frame()
  trend_cum_ins <- data.frame()
  for(i in 1:length(levels(as.factor(x$site)))){
    y <- droplevels(subset(x, site == levels(as.factor(x$site))[i]))
    event_count <- mean(y[, 3], na.rm = T) # Subset by column number so it can apply to MHW and MCS
    event_counts <- rbind(event_counts, event_count)
    event_length <- mean(y[, 8], na.rm = T)
    event_lengths <- rbind(event_lengths, event_length)
    event_cum_in <- mean(y[, 9])
    event_cum_ins <- rbind(event_cum_ins, event_cum_in)
    trend_count <- as.numeric(coef(lm(y[, 3] ~ index(y)))[2])*10
    trend_counts <- rbind(trend_counts, trend_count)
    trend_length <- as.numeric(coef(lm(y[, 8] ~ index(y)))[2])*10
    trend_lengths <- rbind(trend_lengths, trend_length)
    trend_cum_in <- as.numeric(coef(lm(y[, 9] ~ index(y)))[2])*10
    trend_cum_ins <- rbind(trend_cum_ins, trend_cum_in)
  }
  # The mean count of trends for the coastal section
  z <- data.frame(mean_event_count = mean(event_counts[,1], na.rm = T), 
                  mean_event_length = mean(event_lengths[,1], na.rm = T),
                  mean_event_cum_in = mean(event_cum_ins[,1], na.rm = T),
                  mean_trend_count = mean(trend_counts[,1], na.rm = T),
                  mean_trend_length = mean(trend_lengths[,1], na.rm = T),
                  mean_trend_cum_in = mean(trend_cum_ins[,1], na.rm = T))
  return(z)
}

# Stats for MHW and MCS for all coastal sections
results(annualMHW)
results(annualMHWwc)
results(annualMHWsc)
results(annualMHWec)
results(annualMCS)
results(annualMCSwc)
results(annualMCSsc)
results(annualMCSec)

#############################################################################
## Extract two largest MHW and MCS events per coast for further analysis
# Create function to be applied to all MHW and MCS ddataframes
big2 <- function(x){ # To be used with "annual" data frames only
  y1 <- x[which.max(abs(x$Cumulative.intensity..deg.C.x.days.)), 1:15]
  x2 <- x[-which.max(abs(x$Cumulative.intensity..deg.C.x.days.)), ] # Remove largest so second largest can be found
  y2 <- x2[which.max(abs(x2$Cumulative.intensity..deg.C.x.days.)), 1:15]
  z <- rbind(y1, y2)
  return(z)
}

# Big 2 for MHW and MCS for all coastal sections
big2(eventMHW)
big2(eventMHWwc)
big2(eventMHWsc)
big2(eventMHWec)
big2(eventMCS)
big2(eventMCSwc)
big2(eventMCSsc)
big2(eventMCSec)
