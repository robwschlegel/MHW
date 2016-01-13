#############################################################################
## This script does:
# 1. Creates a meta-data table of the daily time series data;
# 2. Creates map to show time series locations;
# 3. Creates graph for each time series showing NA% on Y axis
## Each section may be run independent of the others
#############################################################################

#############################################################################
## DEPENDS ON:
require(plyr); require(zoo); require(lubridate); require(reshape2); require(maptools); require(sp); require(geosphere); require(PBSmapping); require(mapproj); require(ggplot2); require(gridExtra); require(grid); require(viridis); require(scales)
source("func/seqSites.R"); source("func/metaTemp.R"); source("setupParams/theme.R")
# "data/insituDaily_v3.4.RData"
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "metaData.csv"
# "graph/sa_site_NA.pdf"
# "graph/hist_NA.pdf"
#############################################################################

#############################################################################
## Read in daily temperatures created by SACTN scripts
load("data/insituDaily_v3.4.RData")

#############################################################################
## Create meta-data table

# Prep data for calculations
wide <- dcast(insituDaily_v3.4, date ~ site+src, value.var = "temp", mean)
wide_zoo <- zoo(wide[,2:length(colnames(wide))], as.Date(as.POSIXct(wide$date)))

# Calculate meta-data
data.summary <- adply(wide_zoo, 2, metaTemp)

# Clean up for human readability
names(data.summary)[c(1:7)] <- c("index", "start date", "end date", "length", "temp days", "NA days", "NA%")
data.summary$site <- sapply(strsplit(as.character(data.summary$index), "[_]"), "[[", 1)
data.summary$src <- sapply(strsplit(as.character(data.summary$index), "[_]"), "[[", 2)
data.summary <- seqSites(data.summary, data.summary$site)

# Combine data into one meta-data table
siteList <- read.csv("setupParams/site_list_v3.4.csv")
data.summary2 <- cbind(siteList[,2:7], data.summary[,1:10])

# Correctly round the metadata ds based on type
metaData <- data.frame()
for(i in 1:length(levels(as.factor(data.summary2$src)))){
  data <- droplevels(subset(data.summary2, src == levels(as.factor(data.summary2$src))[i]))
  if(levels(as.factor(data$src)) == "SAWS" | levels(as.factor(data$src)) == "KZNSB"){
    data[, 13:16] <- sapply(data[, 13:16], function(x){round_any(x, 0.1)})
    metaData <- rbind(metaData, data)
  }else if(levels(as.factor(data$src)) == "DEA"){
    for(j in 1:length(levels(data$site))){
      data2 <- droplevels(subset(data, site == levels(data$site)[j]))
      if(levels(data2$site) == "Aliwal Shoal" | levels(data2$site) == "Sodwana"){
        metaData <- rbind(metaData, data2)
      }else{
        data[, 13:16] <- sapply(data[, 13:16], function(x){round_any(x, 0.1)})
        metaData <- rbind(metaData, data2)
      }
    }
  }else{
    metaData <- rbind(metaData, data)
  }
}

# Order correctly for use in table
metaData <- seqSites(metaData)
save(metaData, file = "data/metaData.Rdata")
write.csv(metaData, "data/metaData.csv", row.names = F)

#############################################################################
## Create the site map

# Load meta-data
metaData <- read.csv("data/metaData.csv")
colnames(metaData)[1] <- "order"
colnames(metaData)[13] <- "NAper"

# Load SA map data
load("graph/south_africa_coast.RData")
load("graph/sa_provinces_new.RData")
sa_provinces_new$index <- 1:12 # Reduce it by 92%
sa_provinces_new <- droplevels(subset(sa_provinces_new, index == 1))

# Define plotting parameters
sa_lats <- c(-35.5, -26); sa_lons <- c(14, 34)
breaks <- seq(0, 75, 15)

## The scale bar and north arrow
# Create parameters for a scale bar:
scales <-  c(0, 100000, 200000) # Sets initial distances for scale in metres
scaleBar <- data.frame("DISTANCE" = scales/1000, # Changes the distance to look like km
                       "b.lon" = 29.5, # Set beginning longitude point for scale bar
                       "b.lat" = -34.8, # Set beginning latitude point for scale bar
                       destPoint(p = c(29.5, -34.8), b = 90, d = scales)) # Set start point, direction and length of scale bar
scaleLength <- scaleBar[3,] # The ending point

# Create parameters for a North arrow:
nArrow <- data.frame("b.lon" = scaleBar[3,4],
                     "b.lat" = -34.4,
                     destPoint(p = c(scaleBar[3,4], -34.4), b = 0, d = 50000))

# Create map
sa_site <- ggplot() + theme_bw() +
  geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), show_guide = FALSE, fill = "grey80") +
  ### Insert province borders
  geom_path(data = sa_provinces_new, aes(x = long, y = lat, group = group), colour = "grey50") +
  ### Insert site locations and labels
  geom_point(data = metaData, aes(x = lon, y = lat, colour = metaData$NAper), size = 3.2) +
  geom_text(data = metaData, aes(x = lon, y = lat, label = order), size = 1.4, colour = "black", show_guide = FALSE) +
  # Colour management
  scale_colour_continuous(limits = c(0, 75), low = "darkred", high = "pink", breaks = breaks) + # Sets colour palette
  guides(colour = guide_legend(title = "NA%")) +
  ### Annotate specific things:
  annotate("text", label = "Western \n Cape", x = 20.5, y = -33.6, size = 5) +
  annotate("text", label = "Eastern \n Cape", x = 26.5, y = -32.2, size = 5) +
  annotate("text", label = "KwaZulu \n Natal", x = 31.0, y = -28.3, size = 5) +
  annotate("text", label = "False Bay", x = 18.66, y = -34.55, size = 2.2) +
  annotate("text", label = "Algoa Bay", x = 26.4, y = -33.9, size = 2.2) +
  ### Scale bar and N arrow
  # Insert scale bar bottom:
  geom_segment(data = scaleLength, 
               aes(x = b.lon, y = b.lat, xend = lon, yend = b.lat), color = "BLACK", size = 0.3) +
  # Insert scale bar tips:
  geom_segment(data = scaleBar, 
               aes(x = lon, y = b.lat, xend = lon, yend = b.lat + 0.1), color = "BLACK", size = 0.3) +
  # Label distances:
  annotate("text", label = c("0", "100", "200"), 
           x = scaleBar[,4], y = scaleBar[,3] + 0.2, color = "BLACK", size = 2) +
  annotate("text", label = "km", x = scaleBar[3,4] + .2, y = scaleBar[1,3],  color = "BLACK", size = 2) + 
  # Insert N arrow:
  geom_segment(data = nArrow, 
               aes(x = b.lon, y = b.lat, xend = lon, yend = lat), arrow = arrow(length = unit(0.2, "cm"))) +
  annotate("text", label = "N", x = nArrow$b.lon + .15, y = nArrow$b.lat + .1, size = 2) + # Label N arrow
  ### End insert scale bar and arrow.
  #scale_x_continuous(expand = c(0, 0)) +
  coord_map(xlim = sa_lons, ylim = sa_lats, projection = "mercator") +
  theme(panel.background = element_rect(fill = "aquamarine", colour = NA),
        panel.border = element_rect(colour = NA, size = 0.4),
        axis.title = element_blank(),
        #axis.text = element_blank(),
        #axis.ticks = element_blank(),
        panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = "black", size = 0.2, linetype = "dotted"),
        plot.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.key = element_blank(),
        legend.key.height = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.justification = c(1,0), legend.position = c(0.12, 0.07))
#sa_site
ggsave("graph/sa_site_NA.pdf", width = 7.5, height = 4.5, pointsize = 10)

#############################################################################
## NA% fiugure

# Calculate NAs per month for each time series # This is greyed out as it takes a long time to run
# metaMonthNA <- data.frame()
# for(i in 1:length(levels(insituDaily_v3.4$site))){
#   data1 <- droplevels(subset(insituDaily_v3.4, site == levels(insituDaily_v3.4$site)[i]))
#   for(j in 1:length(levels(as.factor(data1$src)))){
#     data2 <- subset(data1, src == levels(as.factor(data1$src))[j])
#     data2$date <- format(data2$date, "%Y-%m")
#     data2 <- na.trim(data2)
#     for(k in 1:length(levels(as.factor(data2$date)))){
#       data3 <- subset(data2, date == levels(as.factor(data2$date))[k])
#       data4 <- length(data3$temp[is.na(data3$temp)])
#       data5 <- data.frame(site = data2$site[1], src = data2$src[1], date = data3$date[1], NA.days = data4)
#       metaMonthNA <- rbind(metaMonthNA, data5)
#     }
#   }
# }
# metaMonthNA <- seqSites(metaMonthNA)
# save(metaMonthNA, file = "data/metaMonthNA.RData")
load("data/metaMonthNA.RData")
# Create ordered factor by site/ src for plotting
metaMonthNA$index <- as.factor(paste(metaMonthNA$site, metaMonthNA$src, sep = "/ "))
siteList <- read.csv("setupParams/site_list_v3.4.csv")
siteList$index <- 1:length(siteList$site)
siteList$index2 <- paste(siteList$site, siteList$src, sep = "/ ")
siteList$index2 <- reorder(siteList$index2, siteList$index)
metaMonthNA2 <- data.frame()
for(i in 1:length(levels(as.factor(metaMonthNA$index)))){
  data1 <- droplevels(subset(metaMonthNA, index == levels(as.factor(siteList$index2))[i]))
  data1$index2 <- i
  metaMonthNA2 <- rbind(metaMonthNA2, data1)
}

metaMonthNA2$index <- reorder(metaMonthNA2$index, metaMonthNA2$index2)
#metaMonthNA2$index <- reorder(metaMonthNA2$index, metaMonthNA$site)

metaMonthNA2$date2 <- as.character(metaMonthNA2$date)
metaMonthNA2$date2 <- parse_date_time(metaMonthNA2$date, "ym", tz = "Africa/Johannesburg")

#metaMonthNA2 <- droplevels(subset(metaMonthNA2, index %in% best$index)) # for use when subsetting from other scripts

# Calculate average missing days
metaMonthNAmean <- ddply(metaMonthNA2, .(site, src), summarize, 
                   sd = round(sd(NA.days, na.rm = TRUE), 2),
                   mean = round(mean(NA.days, na.rm = TRUE), 2))
metaMonthNAmean$index <- paste(metaMonthNAmean$site, metaMonthNAmean$src, sep = "/ ")

# figure
hist_NA <- ggplot(metaMonthNA2, aes(x = date2, y = NA.days, fill = src)) + bw_update +
  geom_bar(stat = "identity") + # Create histogram
  geom_hline(data = metaMonthNAmean, linetype = "solid", size = 0.4, 
             aes(yintercept = mean, colour = src, alpha = 0.6)) +
  facet_wrap(~ index, ncol = 11, scales = "free_x") + # Create figure for each site
  scale_y_discrete(expand = c(0,0), breaks = c(0, 10, 20, 30)) +
  labs(title = NULL, x = NULL, y = expression(paste("NA%"))) +
  scale_fill_discrete(name = "Source", l = 20) +
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "bottom")
hist_NA
ggsave("graph/hist_NA.pdf", width = 15, height = 16, pointsize = 10)
