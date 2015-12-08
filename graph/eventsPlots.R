library(lubridate)
library(scales) # for 'date_format' function
library(plyr)

#############################################################################
## 1. Load and set stuff.
#############################################################################

source("setupParams/theme.R") # also loads ggplot2 etc.

# Decide on which sites should belong to which coasts
wc <- c("Hout Bay", "Kommetjie", "Port Nolloth", "Sea Point")
sc <- c("Fish Hoek", "Gordons Bay", "Hamburg", "Hermanus", "Humewood", "Knysna", "Mossel Bay", "Muizenberg", "Pollock Beach", "Storms River Mouth", "Tsitsikamma", "Ystervarkpunt")
ec <- c("Eastern Beach", "Nahoon Beach", "Orient Beach", "Sodwana")

load("prep/SA_coastal_temps.RData") # These data come directly from the SACTN so are already POSIXCt etc.

mhwWC <- droplevels(SA_coastal_temps[SA_coastal_temps$site %in% wc, ])
mhwSC <- droplevels(SA_coastal_temps[SA_coastal_temps$site %in% sc, ])
mhwEC <- droplevels(SA_coastal_temps[SA_coastal_temps$site %in% ec, ])

mhwWC$date <- as.Date(mhwWC$date) # set dates
mhwSC$date <- as.Date(mhwSC$date)
mhwEC$date <- as.Date(mhwEC$date)

sitesWC <- levels(mhwWC$site) # extract site lists
sitesSC <- levels(mhwSC$site)
sitesEC <- levels(mhwEC$site)

mhwWC$month <- floor_date(mhwWC$date, "month")
monthsWC <- ddply(mhwWC, .(site, month), summarize,
                 temp = mean(temp, na.rm = TRUE)) # make monthly means
mhwSC$month <- floor_date(mhwSC$date, "month")
monthsSC <- ddply(mhwSC, .(site, month), summarize,
                 temp = mean(temp, na.rm = TRUE))
mhwEC$month <- floor_date(mhwEC$date, "month")
monthsEC <- ddply(mhwEC, .(site, month), summarize,
                 temp = mean(temp, na.rm = TRUE))

#############################################################################
## 2. Define column names for event data.
#############################################################################

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

#############################################################################
## 3. Write a function to pull all event data together and produce a data frame;
## also make sure to filter out only the top n MHWs at each site with regards to
## cumulative intensity.
#############################################################################

# TODO: make function more general to allow sorting by any named variable...
#==============================================================================

eventLoad <- function(dir, nCum = 5) {
  fname1 = dir(dir, pattern = "data.events", full.names = TRUE)
  fname2 = dir(dir, pattern = "data.events", full.names = FALSE)
  pf <- str_sub(fname2, -20, -1)
  require(plyr)
  l1 = llply(fname1, read.csv, skip = 2, col.names = colNames) # loops avoided!
  l1 = llply(l1, arrange, -abs(intCum)) # sort by intCum in descending order
  top_n = function(x) head(x, nCum)
  df1 = ldply(l1, top_n) # pick the n highest ones
  require(stringr)
  siteNames1 = unlist(strsplit(dir(dir, pattern = "data.events", full.names = FALSE), pf[1]))
  siteNames1 = str_replace(siteNames1, "_", " ") # parse names for site column
  df1$site = rep(siteNames1, each = nCum)
  df1$month = floor_date(as.Date(paste(df1$yearStrt, df1$monthStrt, df1$dayStrt, sep = "-")), "month")
  df1$site[df1$site == "Storms River_Mouth"] <- "Storms River Mouth" # str_replace() misses the second "_"
  return(df1)
}

dir1 <- paste(getwd(), "/data/MHW/events/", sep = "")
dir2 <- paste(getwd(), "/data/MCS/events/", sep = "")
mhw <- eventLoad(dir1, nCum = 3) # produce a data frame with mhw data...
mhw$event <- rep(1:3, length(levels(as.factor(mhw$site)))) # Make sure to correct the rep() to match nCum
mcs <- eventLoad(dir2, nCum = 3) # produce a data frame with mcs data...
mcs$event <- rep(1:3, length(levels(as.factor(mcs$site)))) 

#############################################################################
## 4. Extract mhw per coast.
#############################################################################

mhwWC <- droplevels(mhw[mhw$site %in% sitesWC, ])
mhwSC <- droplevels(mhw[mhw$site %in% sitesSC, ])
mhwEC <- droplevels(mhw[mhw$site %in% sitesEC, ])

mcsWC <- droplevels(mcs[mcs$site %in% sitesWC, ])
mcsSC <- droplevels(mcs[mcs$site %in% sitesSC, ])
mcsEC <- droplevels(mcs[mcs$site %in% sitesEC, ])

# fullWC <- left_join(x = monthsWC, y = mhwWC, by = c("site", "month"))
# fullSC <- left_join(x = monthsSC, y = mhwSC, by = c("site", "month"))
# fullEC <- left_join(x = monthsEC, y = mhwEC, by = c("site", "month"))

#############################################################################
## 5. Make plots of monthly time series.
#############################################################################

# TODO: make plot function more general to allow auto-naming by any sort variable named under .3...
#==============================================================================

eventPlot <- function(dat, xvar = "month", yvar = "temp", mhw, mcs, width = 6, height = 6) {
    fName <- paste(getwd(),"/graph/", deparse(substitute(dat)), "_plot.pdf", sep = "")
    pdf(fName, width = width, height = height)
      p1 <- ggplot(data = dat, aes_string(x = xvar, y = yvar, group = "site")) +
      geom_line() + bw_update +
      geom_vline(data = mhw, aes(xintercept = as.numeric(month), 
                                 linetype = as.factor(event)), col = "red", size = 0.4) +
      geom_vline(data = mcs, aes(xintercept = as.numeric(month),
                                 linetype = as.factor(event)), col = "blue", size = 0.4) +
      scale_linetype_manual(values = c("solid","dashed", "dotted")) + # Reorder linetypes... not ideal
      scale_x_date(breaks = "1 year", minor_breaks = "1 month", labels = date_format("%Y")) +
      facet_grid(site ~ ., scale = "free_x") +
      ylab(expression(paste("Temperature (", degree~C, ")"))) + xlab("Date") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
      print(p1) # necessary for save to pdf to work...
      dev.off()
      return(p1)
}
eventPlot(dat = monthsWC, mhw = mhwWC, mcs = mcsWC) # plot of West Coast mhw
eventPlot(dat = monthsSC, mhw = mhwSC, mcs = mcsSC, height = 12) # plot of South Coast mhw
eventPlot(dat = monthsEC, mhw = mhwEC, mcs = mcsEC) # plot of East Coast mhw

#Plot all of the time series
monthsALL <- rbind(monthsWC, monthsSC, monthsEC)
mhwALL <- rbind(mhwWC, mhwSC, mhwEC)
mcsALL <- rbind(mcsWC, mcsSC, mcsEC)
eventPlot(dat = monthsALL, mhw = mhwALL, mcs = mcsALL, height = 20)

#############################################################################
## 6. Focus on Hermanus, Mossel Bay and Knysna using anomalies.
#############################################################################

# TODO: Arrange sites by longitude or latitude...
# TODO: Also add Tsitsikamma.
#==============================================================================

monthsSC2 <- droplevels(dplyr::filter(monthsSC, site == c("Hermanus", "Knysna", "Mossel Bay")))
mhwSC2 <- droplevels(mhw[mhw$site %in% c("Hermanus", "Knysna", "Mossel Bay"), ])
mcsSC2 <- droplevels(mcs[mcs$site %in% c("Hermanus", "Knysna", "Mossel Bay"), ])

# TODO: Find a more effective way to find dates in common amongst all four (incl. Tsitsikamma) sites.
#==============================================================================

numberMo <- ddply(monthsSC2,. (month), nrow)
selectMo <- dplyr::filter(numberMo, V1 == 3) # find months shared by all three sites
commonMo <- droplevels(monthsSC2[monthsSC2$month %in% selectMo$month, ])

# TODO: use ddply or something to avoid the loop...
#==============================================================================

meanLT <- ddply(commonMo, .(site), summarize, temp = mean(temp, na.rm = TRUE)) # long-tern mean to subtract
anom <- vector("list", 4)
names(anom) <- levels(monthsSC2$site)
for (n in levels(monthsSC2$site)) {
  anom[[n]] = unlist(monthsSC2$temp[monthsSC2$site == n]) - unlist(meanLT$temp[meanLT$site == n])
} # calculate the anomalies
monthsSC2$anom <- as.vector(unlist(anom))

eventPlot(dat = monthsSC2, mhw = mhwSC2, mcs = mcsSC2, yvar = "anom", width = 6, height = 4) # plot of West Coast mhw
