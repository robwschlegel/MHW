#############################################################################
## This script does:
# 1. loads the raw OISST time series;
# 2. Polishes and saves them for use with other scripts
#############################################################################

#############################################################################
## DEPENDS ON:
require(plyr); require(stringr)
# "data/metaData2.Rdata"
# All of the OISST time series from Eric
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "data/OISST.Rdata"
# "data/OISST.csv"
#############################################################################

#############################################################################
## 1. loads the raw OISST time series

dir <- paste(getwd(), "/data/OISST/", sep = "")
colNames <- c("year", "month", "day", "temp")

fname1 <-  dir(dir, full.names = TRUE)
fname2 <- sapply(strsplit(dir(dir, full.names = FALSE), ".csv"), "[[", 1)
pf <- str_sub(fname2, end = 12)
siteNames1 <- sapply(strsplit(fname2, pf), "[", 2)
siteNames1 <-  str_replace_all(siteNames1, "_", " ")

dat <- data.frame()
for(i in 1:length(fname1)) { # A shameful for loop... in order to label sites correctly
  x <- read.csv(fname1[i], header = FALSE, sep = ",",
                col.names = colNames)
  x$site <-  as.factor(siteNames1[i])
  x$date <-  as.Date(paste(x$year, x$month, x$day, sep = "-"))
  dat <- rbind(dat, x)
}

#############################################################################
## 2. Polishes and saves them for use with other scripts
OISSTdaily <- dat[,c(5,6,4)]

# Check for missing time series from in situ mismatch
load("data/metaData2.Rdata")
metaData2 <- droplevels(metaData2[1:21,]) # Remove coastal aggreagte values

check <- droplevels(metaData2$site[!(metaData2$site %in% OISSTdaily$site)])
check # Currently missing Nahoon Beach and Sodwana

# Insert placeholder rows for missing time series for plotting purposes
placeholder <- data.frame(site = as.factor(c("Nahoon Beach", "Sodwana")), date = OISSTdaily$date[1:2], temp = c(NA,NA))
OISSTdaily <- rbind(OISSTdaily, placeholder)

save(OISSTdaily, file = "data/OISSTdaily.Rdata")
write.csv(OISSTdaily, "data/OISSTdaily.csv")
