# This function is designed to:
# 1. Make calculations on the daily data ("insituDaily_v3.2.RData", produced by "merge2Daily_insitu_v3.2.R"):
# 2. Create basic statistics;

# How long does the temperatures (at a certain threshold) persist?
# How frequently does a temperatures (at a certain threshold) occur?

# EOFs for each month (etc.)

# Search for diurnal pattens in temperatures...


library(plyr)
library(lubridate)
library(zoo) # for the "na.trim()" function; I'd prefer to not use "zoo" but will continue using it for now.
metaTemp <- function(x) {
  a <- as.data.frame(start(na.trim(x))) # Start date of time series
  b <- as.data.frame(end(na.trim(x))) # End date of time series
  c <- length(na.trim(x)) # Length in days of time series
  d <- length(na.omit(na.trim(x))) # Number days with temps
  e <- length(na.trim(x)) - length(na.omit(na.trim(x))) # Number of NaNs
  f <- 100 - ((length(na.omit(na.trim(x))) / length(na.trim(x))) * 100) # % NaNs in time series
  g <- mean(x, na.rm = TRUE)
  #h <- median(x, na.rm = TRUE) # Include as boxplots, not here
  h <- sd(x, na.rm = TRUE)
  i <- min(x, na.rm = TRUE)
  j <- max(x, na.rm = TRUE)
  k <- length(x[x > 18]) # probably interesting for the west coast (kelp area) only... Also, probably better expressed for Feb only - i.e. the percentage of Feb days with a temperature above 18 degr. C.
  all <- cbind(a, b, c, d, e, round(f, 1), round(g, 3), round(h, 3), round(i, 3), round(j, 3), k) # Put them into one data frame
  names(all) <- c("date.start", "date.end", "ts.length", "temp.days", "NaN.days", "NaN.perc",
                  "mean", "sd", "min", "max", "t18")
  return(all)
}