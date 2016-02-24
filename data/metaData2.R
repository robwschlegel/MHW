#############################################################################
## This script does:
# 1. Loads and extracts meta-data for used sites;
# 2. Creates a table of meta-data for each site and coast;
# 3. Creates graph of SA showing the sites and what they are doing
#############################################################################

#############################################################################
## DEPENDS ON:
require(plyr); require(zoo); require(lubridate); require(reshape2); require(maptools); require(sp); require(geosphere); require(PBSmapping); require(mapproj); require(ggplot2); require(gridExtra); require(grid); require(viridis); require(scales)
source("func/seqSites.R"); source("func/metaTemp.R"); source("setupParams/theme.R")
# "metaData.csv""
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "data/metaData2.Rdata"
# "data/metaData2.csv"
# "graph/figure01.pdf"
#############################################################################

#############################################################################
## Read in metaData and extract used sites
load("data/metaData.Rdata")
names(metaData)[2] <- "coast" # Remove src column as type already shows this
metaData2 <- metaData[metaData$`NA%` <= 10,]
metaData2 <- metaData2[metaData2$length >= 3650,]
metaData2$site <- as.character(metaData2$site)
metaData2[12,1] <- "Tsitsikamma West"
metaData2[14,1] <- "Tsitsikamma East"
metaData2$site <- as.factor(metaData2$site)

#############################################################################
## Calculate mean meta data for each coast and add to metaData2
# Dictate coastal groupings
wc <- c("Hout Bay", "Kommetjie", "Port Nolloth", "Sea Point")
sc <- c("Fish Hoek", "Gordons Bay", "Hamburg", "Hermanus", "Humewood", "Knysna", "Mossel Bay", "Muizenberg", "Pollock Beach", "Storms River Mouth", "Tsitsikamma East", "Tsitsikamma West", "Ystervarkpunt")
ec <- c("Eastern Beach", "Nahoon Beach", "Orient Beach", "Sodwana")

# Extract sites and create single line data.frame for each coast
west <- metaData2[metaData2$site %in% wc,]
west$coast <- "west"
west2 <- data.frame(site = "west coast", coast = "west", lon = round(mean(west$lon), 4),
                    lat = round(mean(west$lat), 4), depth = round(mean(west$depth), 1),
                   type = NA, 'start date' = min(west$`start date`),'end date' = max(west$`end date`), 
                   length = round(mean(west$length), 0), 'temp days' = round(mean(west$'temp days'), 0), 
                   'NA days' = round(mean(west$'NA days'), 0), 'NA%' = round(mean(west$`NA%`),1), 
                   mean = round(mean(west$mean),1), sd = round(mean(west$sd),1), 
                   min = round(mean(west$min),1), max = round(mean(west$max),1))

south <- metaData2[metaData2$site %in% sc,]
south$coast <- "south"
south2 <- data.frame(site = "south coast", coast = "south", lon = round(mean(south$lon), 4),
                    lat = round(mean(south$lat), 4), depth = round(mean(south$depth), 1),
                    type = NA, 'start date' = min(south$`start date`),'end date' = max(south$`end date`), 
                    length = round(mean(south$length), 0), 'temp days' = round(mean(south$'temp days'), 0), 
                    'NA days' = round(mean(south$'NA days'), 0), 'NA%' = round(mean(south$`NA%`),1), 
                    mean = round(mean(south$mean),1), sd = round(mean(south$sd),1), 
                    min = round(mean(south$min),1), max = round(mean(south$max),1))

east <- metaData2[metaData2$site %in% ec,]
east$coast <- "east"
east2 <- data.frame(site = "east coast", coast = "east", lon = round(mean(east$lon), 4),
                    lat = round(mean(east$lat), 4), depth = round(mean(east$depth), 1),
                    type = NA, 'start date' = min(east$`start date`),'end date' = max(east$`end date`), 
                    length = round(mean(east$length), 0), 'temp days' = round(mean(east$'temp days'), 0), 
                    'NA days' = round(mean(east$'NA days'), 0), 'NA%' = round(mean(east$`NA%`),1), 
                    mean = round(mean(east$mean),1), sd = round(mean(east$sd),1), 
                    min = round(mean(east$min),1), max = round(mean(east$max),1))

coasts <- rbind(west2, south2, east2)
names(coasts)[c(7:8,10:12)] <- c("start date", "end date", "temp days", "NA days", "NA%")
metaData2 <- rbind(west, south, east, coasts)
metaData2 <- data.frame(ID = seq(1, length(metaData2$site),1), metaData2)
metaData2$lon <- round(metaData2$lon,2); metaData2$lat <- round(metaData2$lat,2)
metaData2$depth <- NULL; metaData2$NA.days <- NULL

save(metaData2, file = "data/metaData2.Rdata")
write.csv(metaData2, "data/metaData2.csv", row.names = F)

#############################################################################
## Create figure showing position of sites

# Load SA map data
load("graph/south_africa_coast.RData")

# Define plotting parameters
sa_lats <- c(-35.5, -26); sa_lons <- c(14, 34)

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
sa_site <- ggplot() + bw_update +
  geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), 
               show_guide = FALSE, fill = "grey80") +
  ### Insert province borders
  #geom_path(data = sa_provinces_new, aes(x = long, y = lat, group = group), colour = "grey50") +
  ### Insert site locations and labels
  geom_point(data = metaData2, aes(x = lon, y = lat, colour = metaData2$coast, 
                                   shape = metaData2$type), alpha = 0.8, size = 2.2) +
  #geom_text(data = metaData, aes(x = lon, y = lat, label = order), size = 1.4, colour = "black", show_guide = FALSE) +
  # Colour management
  #scale_colour_continuous(limits = c(0, 75), low = "darkred", high = "pink", breaks = breaks) + # Sets colour palette
  guides(colour = guide_legend(title = "coast"),
         shape = guide_legend(title = "type")) +
  ### Annotate specific things:
  #annotate("text", label = "Western \n Cape", x = 20.5, y = -33.6, size = 5) +
  #annotate("text", label = "Eastern \n Cape", x = 26.5, y = -32.2, size = 5) +
  #annotate("text", label = "KwaZulu \n Natal", x = 31.0, y = -28.3, size = 5) +
  annotate("text", label = "Cape Agulhas", x = 20.0, y = -35.1, size = 2.2) +
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
  theme(panel.background = element_rect(fill = NA, colour = NA),
        panel.border = element_rect(colour = "black", size = 0.4),
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
ggsave("graph/figure01.pdf", width = 7.5, height = 4.5, pointsize = 10)