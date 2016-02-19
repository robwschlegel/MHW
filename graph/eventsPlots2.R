#############################################################################
## This script does:
# 1. loads daily coastal and SST time series;
# 2. loads coastal and SST MHW/ MCS results;
# 3. Calculates co-occurrence of events based on size of event and lag;
# 4. Extracts top three MHW/ MCS events for each time series and type of data;
# 5. Prepares data for plotting;
# 6. Creates map of co-occurrence values;
# 7. Prepares data for plotting;
# 8. Creates massive line graph of all time series with top three events;
# 9. Prepares data for plotting;
# 10. Creates dot and line graph of co-occurrence depending on quantiles of events and lag
#############################################################################

#############################################################################
## DEPENDS ON:
require(ggplot2); require(stringr); require(plyr); require(zoo); require(lubridate)
source("setupParams/theme.R")
# "prep/SA_coastal_temps.RData"
# "data/metaData2.Rdata"
# "graph/sa_shore.Rdata"
# All of the files from Eric
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "data/mhwCO.csv"
# "data/mhwCOb.csv"
# "data/mhwCOa.csv"
# "data/mcsCO.csv"
# "data/mcsCOb.csv"
# "data/mcsCOa.csv"
# "data/mhwnCO.csv"
# "data/mcsnCO.csv"
# "data/mhwnCOb.csv"
# "data/mcsnCOb.csv"
# "data/mhwnCOa.csv"
# "data/mcsnCOa.csv"
# "graph/mhwCOmap.pdf"
# "graph/mhwCObmap.pdf"
# "graph/mhwCOamap.pdf"
# "graph/mcsCOmap.pdf"
# "graph/mcsCObmap.pdf"
# "graph/mcsCOamap.pdf"
# "graph/eventsALL.pdf"
# "graph/mhwnCOfig.pdf"
# "graph/mcsnCOfig.pdf"
# "graph/mhwnCObfig.pdf"
# "graph/mcsnCObfig.pdf"
# "graph/mhwnCOafig.pdf"
# "graph/mcsnCOafig.pdf"
#############################################################################

#############################################################################
## 1. loads daily coastal and SST time series
load("prep/SA_coastal_temps.RData")
#levels(SA_coastal_temps$site) # Check that the correct sites are being used

# Load satellite time series
load("data/OISSTdaily.Rdata")
#levels(OISSTdaily$site)

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
    x$month <- floor_date(as.Date(paste(x$yearStrt, x$monthStrt, x$dayStrt, sep = "-")), "month")
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
mhw$type <- "insitu"
mcs <- eventLoad(dir2) # produce a data frame with mcs data...
mcs$type <- "insitu"
mhwSST <- eventLoad(dir3) # produce a data frame with SST mhw data...
mhwSST$type <- "sst"
mcsSST <- eventLoad(dir4) # produce a data frame with SST mcs data...
mcsSST$type <- "sst"

#############################################################################
## 3. Calculates co-occurrence of events based on quantile of event and lag

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
# mhwCO <- cooccurrence(mhw, mhwSST) # This takes several minutes to run... cursed for loops...
# write.csv(mhwCO, "data/mhwCO.csv", row.names = F)
mhwCO <- read.csv("data/mhwCO.csv")
#MCS
# mcsCO <- cooccurrence(mcs, mcsSST)
# write.csv(mcsCO, "data/mcsCO.csv", row.names = F)
mcsCO <- read.csv("data/mcsCO.csv")

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
mhwn$type <- "insitu"
mcsn <- eventLoadn(dir2, nCum = 3) # produce a data frame with mcs data...
mcsn$event <- rep(1:3, length(levels(as.factor(mcsn$site))))
mcsn$type <- "insitu"
# SST
mhwnSST <- eventLoadn(dir3, nCum = 3)
mhwnSST$event <- rep(1:3, length(levels(as.factor(mhwnSST$site))))
mhwnSST$type <- "sst"
mcsnSST <- eventLoadn(dir4, nCum = 3)
mcsnSST$event <- rep(1:3, length(levels(as.factor(mcsnSST$site))))
mcsnSST$type <- "sst"

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

### NB: The cooccurrence values have been calculated differently so this will no longer run...

# Create graphing function
# cooccurrenceMap <- function(dat){
#   p <- ggplot() + coord_equal() + theme_bw() +
#     geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), 
#                  colour = "black", fill = "grey80") +
#     geom_point(data = dat, aes(x = lon, y = lat, colour = proportion), size = 3.2) +
#     scale_colour_continuous(limits = c(0.0, 0.8), low = "blue", high = "red", breaks = seq(0.2, 0.6, 0.2)) +
#     coord_map(xlim = lonSA, ylim = latSA, projection = "mercator") +
#     theme(legend.background = element_blank(),
#           legend.justification = c(1,0), legend.position = c(0.18, 0.08))
#   p
# }
# 
# # Create all the graphs
# mhwCOmap <- cooccurrenceMap(mhwCO)
# ggsave("graph/mhwCOmap.pdf", width = 8, height = 6)
# mhwCObmap <- cooccurrenceMap(mhwCOb)
# ggsave("graph/mhwCObmap.pdf", width = 8, height = 6)
# mhwCOamap <- cooccurrenceMap(mhwCOa)
# ggsave("graph/mhwCOamap.pdf", width = 8, height = 6)
# mcsCOmap <- cooccurrenceMap(mcsCO)
# ggsave("graph/mcsCOmap.pdf", width = 8, height = 6)
# mcsCObmap <- cooccurrenceMap(mcsCOb)
# ggsave("graph/mcsCObmap.pdf", width = 8, height = 6)
# mcsCOamap <- cooccurrenceMap(mcsCOa)
# ggsave("graph/mcsCOamap.pdf", width = 8, height = 6)

#############################################################################
## 7. Prepares data for plotting

coastalTemp <- SA_coastal_temps[,c(1,3,4)]
coastalTemp$date <- as.Date(coastalTemp$date)
coastalTemp$type <- "insitu"
OISSTdaily$type <- "sst"

# Combine and reduce to monthly values
tsALL <- rbind(coastalTemp, OISSTdaily)
tsALL$date <- floor_date(tsALL$date, "month")
tsALL <- ddply(tsALL, .(site, type, date), summarize,
                         temp = mean(temp, na.rm = TRUE))

# Reorder sites for plotting
siteOrder <- c("Port Nolloth", "Sea Point", "Hout Bay", "Kommetjie", "Fish Hoek", "Muizenberg", "Gordons Bay", "Hermanus", "Ystervarkpunt", "Mossel Bay", "Knysna", "Tsitsikamma West", "Storms River Mouth", "Tsitsikamma East", "Pollock Beach", "Humewood", "Hamburg", "Eastern Beach", "Orient Beach", "Nahoon Beach", "Sodwana")

tsALL$site <- factor(tsALL$site, levels = siteOrder)
mhwn$site <- factor(mhwn$site, levels = siteOrder)
mcsn$site <- factor(mcsn$site, levels = siteOrder)
mhwnSST$site <- factor(mhwnSST$site, levels = siteOrder)
mcsnSST$site <- factor(mcsnSST$site, levels = siteOrder)

#############################################################################
## 8. Creates massive line graph of all time series with top three events

p1 <- ggplot(data = tsALL, aes(x = date, y = temp)) + bw_update +
  geom_line() +
  geom_vline(data = mhwn, aes(xintercept = as.numeric(month), linetype = as.factor(event)), col = "red", size = 0.4) +
  geom_vline(data = mcsn, aes(xintercept = as.numeric(month), linetype = as.factor(event)), col = "blue", size = 0.4) +
  geom_vline(data = mhwnSST, aes(xintercept = as.numeric(month), linetype = as.factor(event)), col = "red", size = 0.4) +
  geom_vline(data = mcsnSST, aes(xintercept = as.numeric(month), linetype = as.factor(event)), col = "blue", size = 0.4) +
  scale_linetype_manual(values = c("solid","dashed", "dotted"), guide = FALSE) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.015,0)) +
  facet_grid(site ~ type) +
  ylab(expression(paste("Temperature (", degree~C, ")"))) + xlab("Date") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p1

ggsave("graph/figure2.pdf", height = 24, width = 12)

#############################################################################
# 9. Prepares data for plotting

mhwCO$site <- factor(mhwCO$site, levels = siteOrder)
mhwCO$direction <- factor(mhwCO$direction, levels = c("b", "x", "a"))
mhwCO$index <- paste(mhwCO$lag, mhwCO$direction, sep = "_")
mcsCO$site <- factor(mcsCO$site, levels = siteOrder)
mcsCO$direction <- factor(mcsCO$direction, levels = c("b", "x", "a"))
mcsCO$index <- paste(mcsCO$lag, mcsCO$direction, sep = "_")

#############################################################################
# 11. Creates dot and line graph of co-occurrence depending on quantiles of events

## Test for issues
test <- mcsCO[mcsCO$site == "Gordons Bay",]
test$index <- paste(test$lag, test$direction, sep = "_")
test <- droplevels(test[test$direction == "x",])

p2 <- ggplot(data = test, aes(x = quantile, y = proportion)) + bw_update +
  geom_line(aes(colour = as.factor(index))) + 
  #geom_point(aes(colour = as.factor(test$lag))) +
  facet_grid(site ~ direction) +
  ylab("proportion") + xlab("quantile (%)")
p2


cooccurrenceQuantFigure <- function(dat){
  p2 <- ggplot(data = dat, aes(x = quantile, y = proportion)) + bw_update +
    geom_line(aes(colour = as.factor(index))) + 
    geom_point(aes(colour = as.factor(index))) +
    facet_grid(site ~ direction) +
    ylab("proportion") + xlab("quantile (%)") #+
    #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p2
}

mhwCOfig <- cooccurrenceQuantFigure(mhwCO)
ggsave("graph/mhwCOfig.pdf", width = 16, height = 18)
mcsCOfig <- cooccurrenceQuantFigure(mcsCO)
ggsave("graph/mcsCOfig.pdf", width = 16, height = 18)



## Note currently focussing on individual coastlines
# eventPlot(dat = monthsWC, mhw = mhwWC, mcs = mcsWC) # plot of West Coast mhw
# eventPlot(dat = monthsSC, mhw = mhwSC, mcs = mcsSC) # plot of South Coast mhw
# eventPlot(dat = monthsEC, mhw = mhwEC, mcs = mcsEC) # plot of East Coast mhw

## Calculate mean from annual time series
  ## Use the .annual results to calculate other statistics
  ## Count, dur, mean intensity
    ## From this one can deduce cummulative intensity

## If things match up/ coupling this shows ability to exchange water with offshore bodies
## If no match up, this means a lack of ability to exchange water with offshore bodies
