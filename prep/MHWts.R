#############################################################################
## This script does:
# 1. loads daily temperature time series;
# 2. Removes time series with too much NA%;
# 3. Saves for Eric Oliver to run his Python MHW script on
#############################################################################

#############################################################################
## DEPENDS ON:
source("func/seqSites.R")
# "data/SACTNdaily_v4.0.Rdata"
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
#"prep/SA_coastal_temps.RData"
#"prep/SA_coastal_temps.csv"
#"prep/SA_coastal_meta.RData"
#"prep/SA_coastal_meta.csv"
#############################################################################

#############################################################################
## Read in daily temperatures created by SACTN scripts and meta-data for it
load("data/SACTNdaily_v4.0.Rdata")
metaData <- read.csv("data/metaData.csv")
siteList <- read.csv("setupParams/site_list_v3.4.csv")

#############################################################################
## Investigate the state of NA% in data
NArundown <- data.frame()
for(i in 1:10){
  data1 <- data.frame(NA. = i*10, n = length(metaData$NA.[metaData$NA. <= i*10]))
  NArundown <- rbind(NArundown, data1)
}

# Best time series as decided by NA%
best <- metaData[metaData$NA. < 10, ]
mean(best$length)
max(best$length)
min(best$length)

# Then refine by length
best <- best[best$length >= 3650, ]
mean(best$length)
max(best$length)
min(best$length)

# Add index for subsetting
metaData$index <- paste(metaData$site, metaData$src, sep = "/")
best$index <- paste(best$site, best$src, sep = "/")
SACTNdaily_v4.0$index <- paste(SACTNdaily_v4.0$site, SACTNdaily_v4.0$src, sep = "/")
siteList$index <- paste(siteList$site, siteList$src, sep = "/")

# Add lon/ lat and order along coast
# SA_coastal_temps <- data.frame()
# for(i in 1:length(levels(as.factor(insitu$index)))){
#   data1 <- droplevels(subset(insitu, index == levels(as.factor(insitu$index))[i]))
#   coords <- droplevels(subset(siteList, index == data1$index[1]))
#   data1 <- cbind(data1, coords[,4:5])
#   SA_coastal_temps <- rbind(SA_coastal_temps, data1)
# }
#SA_coastal_temps <- seqSites(insitu) # Ordering the sites no longer works when the names of Tsitsikamma are changed

# Subset data and meta-data
SA_coastal_temps <- SACTNdaily_v4.0[SACTNdaily_v4.0$index %in% best$index, ]
SA_coastal_meta <- metaData[metaData$index %in% best$index, ]
#SA_coastal_meta <- seqSites(SA_coastal_meta)

# Differentiate between the two different Tsitsikamma sites
SA_coastal_temps$site <- as.character(SA_coastal_temps$site)
SA_coastal_temps$site[SA_coastal_temps$site == "Tsitsikamma" & SA_coastal_temps$src == "SAWS"] <- "Tsitsikamma West"
SA_coastal_temps$site[SA_coastal_temps$site == "Tsitsikamma" & SA_coastal_temps$src == "DEA"] <- "Tsitsikamma East"
SA_coastal_temps$site <- as.factor(SA_coastal_temps$site)
levels(SA_coastal_temps$site)

SA_coastal_meta$site <- as.character(SA_coastal_meta$site)
SA_coastal_meta$site[SA_coastal_meta$site == "Tsitsikamma" & SA_coastal_meta$src == "SAWS"] <- "Tsitsikamma West"
SA_coastal_meta$site[SA_coastal_meta$site == "Tsitsikamma" & SA_coastal_meta$src == "DEA"] <- "Tsitsikamma East"
SA_coastal_meta$site <- as.factor(SA_coastal_meta$site)
levels(SA_coastal_meta$site)

# Save
save(SA_coastal_temps, file = "prep/SA_coastal_temps.RData")
write.csv(SA_coastal_temps, file = "prep/SA_coastal_temps.csv", row.names = F)
save(SA_coastal_meta, file = "prep/SA_coastal_meta.RData")
write.csv(SA_coastal_meta, file = "prep/SA_coastal_meta.csv", row.names = F)
