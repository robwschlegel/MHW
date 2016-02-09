#############################################################################
## This script does:
# 1. loads daily coastal and SST time series;
# 2. loads coastal and SST MHW/ MCS results;
# 3. Calculates co-occurrence of events between coastal and SST time series;
# 4. Extracts top three MHW/ MCS events for each time series and type of data;
# 5. Prepares data for plotting;
# 6. Creates map of co-occurrence values;
# 7. Creates massive line graph of all time series with top three events
#############################################################################

#############################################################################
## DEPENDS ON:
require(ggplot2); require(stringr); require(plyr); require(lubridate)
source("setupParams/theme.R")
# "prep/SA_coastal_temps.RData"
# "data/metaData2.Rdata"
# All of the files from Eric
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "data/mhwCO.csv"
# "data/mhwCOm.csv"
# "data/mcsCO.csv"
# "data/mcsCOm.csv"
# "graph/mhwCOmap.pdf"
# "graph/mhwCOmmap.pdf"
# "graph/mcsCOmap.pdf"
# "graph/mcsCOmmap.pdf"
#############################################################################

#############################################################################
## 1. loads daily coastal and SST time series
load("prep/SA_coastal_temps.RData")
levels(SA_coastal_temps$site) # Check that the correct sites are being used

# Load metadata to get the lon/ lat values later
load("data/metaData2.Rdata")


#############################################################################
## 2. loads coastal and SST MHW/ MCS results

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
    x$date <-  as.Date(paste(x$yearStrt, x$monthStrt, x$dayStrt, sep = "-"))
    x$lon <- metaData2$lon[metaData2$site == x$site[1]]
    x$lat <- metaData2$lat[metaData2$site == x$site[1]]
    dat <- rbind(dat, x)
  }
  return(dat)
}

dir1 <- paste(getwd(), "/data/MHW/events/", sep = "")
dir2 <- paste(getwd(), "/data/MCS/events/", sep = "")
dir3 <- paste(getwd(), "/data/MHW/SST events/", sep = "")
dir4 <- paste(getwd(), "/data/MCS/SST events/", sep = "")
mhw <- eventLoad(dir1) # produce a data frame with mhw data...
mcs <- eventLoad(dir2) # produce a data frame with mcs data...
mhwSST <- eventLoad(dir3) # produce a data frame with SST mhw data...
mcsSST <- eventLoad(dir4) # produce a data frame with SST mcs data...

#############################################################################
## 3. Calculates co-occurrence of events between coastal and SST time series

cooccurrence <- function(dat1, dat2, lag = 14) {
  dat3 <- data.frame()
  for(i in 1:length(levels(as.factor(dat1$site)))) {
    x1 <- droplevels(subset(dat1, site == levels(as.factor(dat1$site))[i]))
    x2 <- droplevels(subset(dat2, site == levels(as.factor(dat1$site))[i]))
    y <- 0
    for(j in 1:nrow(x1)) {
      x1.1 <- x1$date[j]
      x1.2 <- x1.1 - days(lag)
      if(lag > 0) {
        x1.3 <- seq(x1.2, x1.1, 1)
      } else if(lag < 0) {
        x1.3 <- seq(x1.1, x1.2, 1)
      } else if(lag == 0){
        x1.3 <- x1.1
      }
      x2.1 <- droplevels(subset(x2, date %in% x1.3))
      if(nrow(x2.1) > 0){
        y <- y + 1
      }
    }
    z <- data.frame(site = x1$site[1], events = nrow(x1), 
                    cooccurrence = y, proportion = y/nrow(x1),
                    lon = x1$lon[1], lat = x1$lat[1])
    dat3 <- rbind(dat3, z)
  }
  return(dat3)
}

mhwCO <- cooccurrence(mhw, mhwSST)
write.csv(mhwCO, "data/mhwCO.csv", row.names = F)
mcsCO <- cooccurrence(mcs, mcsSST)
write.csv(mcsCO, "data/mcsCO.csv", row.names = F)
mhwCOm <- cooccurrence(mhw, mhwSST, lag = -14)
write.csv(mhwCOm, "data/mhwCOm.csv", row.names = F)
mcsCOm <- cooccurrence(mcs, mcsSST, lag = -14)
write.csv(mcsCOm, "data/mcsCOm.csv", row.names = F)


#############################################################################
## 4. Extracts top three MHW/ MCS events for each time series and type of data

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
mhwn <- eventLoadn(dir1, nCum = 3) # produce a data frame with mhw data...
mhwn$event <- rep(1:3, length(levels(as.factor(mhwn$site)))) # Make sure to correct the rep() to match nCum
mcsn <- eventLoadn(dir2, nCum = 3) # produce a data frame with mcs data...
mcsn$event <- rep(1:3, length(levels(as.factor(mcsn$site))))
# SST
mhwnSST <- eventLoadn(dir3, nCum = 3)
mhwnSST$event <- rep(1:3, length(levels(as.factor(mhwnSST$site))))
mcsnSST <- eventLoadn(dir2, nCum = 3)
mcsnSST$event <- rep(1:3, length(levels(as.factor(mcsnSST$site)))) 

# mhwWC <- droplevels(mhw[mhw$site %in% sitesWC, ]) # Currently not focussing on coastal values
# mhwSC <- droplevels(mhw[mhw$site %in% sitesSC, ])
# mhwEC <- droplevels(mhw[mhw$site %in% sitesEC, ])
# 
# mcsWC <- droplevels(mcs[mcs$site %in% sitesWC, ])
# mcsSC <- droplevels(mcs[mcs$site %in% sitesSC, ])
# mcsEC <- droplevels(mcs[mcs$site %in% sitesEC, ])

#############################################################################
## 5. Prepares data for plotting

# Setup up environment for plotting
latSA <- c(-35.5, -26); lonSA <- c(14, 34)

# Load SA map data
load("graph/south_africa_coast.RData")

#############################################################################
## 6. Creates map of co-occurrence values

# Create graphing function
cooccurrenceMap <- function(dat){
  p <- ggplot() + coord_equal() + theme_bw() +
    geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), 
                 colour = "black", fill = "grey80") +
    geom_point(data = dat, aes(x = lon, y = lat, colour = proportion), size = 3.2) +
    scale_colour_continuous(limits = c(0.0, 0.5), low = "blue", high = "red", breaks = seq(0.0, 0.5, 0.1)) +
    coord_map(xlim = lonSA, ylim = latSA, projection = "mercator") +
    theme(legend.background = element_blank(),
          legend.justification = c(1,0), legend.position = c(0.18, 0.08))
  p
}

# Create all the graphs
mhwCOmap <- cooccurrenceMap(mhwCO)
ggsave("graph/mhwCOmap.pdf", width = 8, height = 6)
mhwCOmmap <- cooccurrenceMap(mhwCOm)
ggsave("graph/mhwCOmmap.pdf", width = 8, height = 6)
mcsCOmap <- cooccurrenceMap(mcsCO)
ggsave("graph/mcsCOmap.pdf", width = 8, height = 6)
mcsCOmmap <- cooccurrenceMap(mcsCOm)
ggsave("graph/mcsCOmmap.pdf", width = 8, height = 6)

#############################################################################
## 7. Creates massive line graph of all time series with top three events

eventPlot <- function(dat, xvar = "month", yvar = "temp", mhw, mcs, width = 6, height = 6) {
  fName <- paste("graphs/", deparse(substitute(dat)), "_plot.pdf", sep = "")
  pdf(fName, width = width, height = height)
  p1 <- ggplot(data = dat, aes_string(x = xvar, y = yvar, group = "site")) +
    geom_line() +
    geom_vline(data = mhw, aes(xintercept = as.numeric(month)), col = "red", size = 0.4) +
    geom_vline(data = mcs, aes(xintercept = as.numeric(month)), col = "blue", size = 0.4) +
    scale_x_date(breaks = "1 year", minor_breaks = "1 month", labels = date_format("%Y")) +
    facet_grid(site ~ ., scale = "free_x") +
    ylab(expression(paste("Temperature (", degree~C, ")"))) + xlab("Date") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  print(p1) # necessary for save to pdf to work...
  dev.off()
  return(p1)
}
eventPlot(dat = monthsWC, mhw = mhwWC, mcs = mcsWC) # plot of West Coast mhw
eventPlot(dat = monthsSC, mhw = mhwSC, mcs = mcsSC) # plot of South Coast mhw
eventPlot(dat = monthsEC, mhw = mhwEC, mcs = mcsEC) # plot of East Coast mhw

## Calculate mean from annual time series
  ## Use the .annual results to calculate other statistics
  ## Count, dur, mean intensity
    ## From this one can deduce cummulative intensity

## If things match up/ coupling this shows ability to exchange water with offshore bodies
## If no match up, this means a lack of ability to exchange water with offshore bodies

## Calculate probability of co-occurence between in situ SST time series 
  ## Include lag to help expand the range of possibility
## Create a map showing the fraction of co-occurence happening at each site