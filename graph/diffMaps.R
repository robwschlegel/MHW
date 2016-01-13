#############################################################################
## This script does:
# 1. Loads and combines metadata and MHW and MCS data for all sites;
# 2. Creates a function that creates comparison figures;
# 3. Creates graph of each statistics
#############################################################################

#############################################################################
## DEPENDS ON:
require(plyr); require(zoo); require(lubridate); require(reshape2); require(maptools); require(sp); require(geosphere); require(PBSmapping); require(mapproj); require(ggplot2); require(gridExtra); require(grid); require(viridis); require(scales)
source("func/seqSites.R"); source("func/metaTemp.R"); source("setupParams/theme.R")
# "metaData2.csv"
# "data/allSitesMHWMCS.Rdata"
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "graph/diff.pdf"
#############################################################################

#############################################################################
## Read in metaData and MHW/ MCS stats
load("data/metaData2.Rdata")
load("data/allSitesMHWMCS.Rdata")

# Manually fix mismatch in site numbers for now... 
  ## This must be corrected for by deciding which sites to use and how...
metaData2 <- metaData2[c(1:13,15:21),]
metaData2$site <- as.character(metaData2$site)
metaData2[12,1] <- "Tsitsikamma"
metaData2$site <- as.factor(metaData2$site)

# Calculate differences between MHWs and MCSs
diff <- data.frame()
for(i in 1:20){
  x <- data.frame(site = allSitesMHWMCS$site[i], 
                  allSitesMHWMCS[i,2:7] - allSitesMHWMCS[i+20,2:7], event = "diff")
  diff <- rbind(diff, x)
}

allSitesMHWMCS <- rbind(allSitesMHWMCS, diff)

# Squish it all together
metaData3 <- rbind(metaData2, metaData2, metaData2)
allData <- cbind(allSitesMHWMCS, metaData3)

#############################################################################
## Create function that will create paneled figures

# Load SA map data
load("graph/south_africa_coast.RData")

# Define plotting parameters
sa_lats <- c(-35.5, -26); sa_lons <- c(14, 34)

diffMap <- function(x){ # x = the column containing the statisitic to be displayed
  allData2 <- allData[,c(1,x,8:20)]
  allData2$stat <- allData2[,2]
  m <- ggplot() + bw_update +
    geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), 
                 show_guide = FALSE, fill = "grey20") +
    ### Insert site locations
    geom_point(data = allData2, aes(x = lon, y = lat, colour = stat), alpha = 0.8, size = 2.2) +
    # Colour management
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + # Sets colour palette
    guides(colour = guide_legend(title = colnames(allData2[2]))) +
    ### Annotate specific things:
    #annotate("text", label = "Cape Agulhas", x = 20.0, y = -35.1, size = 2.2) +
    #annotate("text", label = "False Bay", x = 18.66, y = -34.55, size = 2.2) +
    #annotate("text", label = "Algoa Bay", x = 26.4, y = -33.9, size = 2.2) +
    #scale_x_continuous(expand = c(0, 0)) +
    labs(title = NULL, x = NULL, y = NULL) +
    coord_map(xlim = sa_lons, ylim = sa_lats, projection = "mercator") +
    theme(panel.background = element_rect(fill = "grey80")) +
    facet_wrap(~event, nrow = 1)
  #m
}

#############################################################################
## Create figures

diff_event_count <- diffMap(2)
ggsave("graph/dif_event_count.pdf", width = 10, height = 3, pointsize = 8)

diff_event_length <- diffMap(3)
ggsave("graph/dif_event_length.pdf", width = 10, height = 3, pointsize = 8)

diff_event_cum_in <- diffMap(4)
ggsave("graph/dif_event_cum_in.pdf", width = 10, height = 3, pointsize = 8)

diff_trend_count <- diffMap(5)
ggsave("graph/dif_trend_count.pdf", width = 10, height = 3, pointsize = 8)

diff_trend_length <- diffMap(6)
ggsave("graph/dif_trend_length.pdf", width = 10, height = 3, pointsize = 8)

diff_trend_cum_in <- diffMap(7)
ggsave("graph/dif_trend_cum_in.pdf", width = 10, height = 3, pointsize = 8)
