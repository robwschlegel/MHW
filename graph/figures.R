#############################################################################
## This script does:
# 1. Figure 1
# 2. Figure 2
# 3. Figure 3
# 4. Figure 4
# 5. Figure 5
#############################################################################

#############################################################################
## DEPENDS ON:
require(ggplot2); require(lubridate); require(plyr); require(ggrepel)
source("setupParams/theme.R"); source("proc/results2.R")
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
# "graph/figure1.pdf"
# "graph/figure2.pdf"
# "graph/figure4.pdf"
# "graph/figure5.pdf"
#############################################################################

#############################################################################
## 1. Figure 1
# The coastal figure with sites

# Load meatadata and site pixels
metaData2 <- metaData2[-((length(metaData2$site)-2):length(metaData2$site)),] # Remove coastal means
metaData2$coast <- factor(metaData2$coast, levels = c("west", "south", "east"))
load("data/site_pixels.RData")
site_pixels2 <- read.csv("data/site_pixels2.csv")
site_pixels[52:53,1:3] <- site_pixels2
# Add coastal values to 
site_pixels3 <- data.frame()
for(i in 1:nrow(site_pixels)){
  x <- site_pixels[i,]
  if(x$site[1] %in% wc) {
    x$coast <- "west"
  } else if(x$site[1] %in% sc) {
    x$coast <- "south"
  } else if(x$site[1] %in% ec) {
    x$coast <- "east"
  }
  site_pixels3 <- rbind(site_pixels3, x)
}
site_pixels3$coast <- factor(site_pixels3$coast, levels = c("west", "south", "east"))

# Setup up environment for plotting
sa_lats <- c(-37, -27); sa_lons <- c(14, 34)

# Load SA map data and bathymetry
load("graph/south_africa_coast.RData")
names(south_africa_coast)[1] <- "lon"
load("data/bathy/bathy.RData") # HiRes for final image
load("data/bathy/sa_bathy.RData") # LowRes for tweaking

# 1. Figure 1
f1 <- ggplot() + theme_bw() + #coord_equal() + 
  geom_contour(data = bathy[bathy$depth >= -399,], aes(x = lon, y = lat, z = depth), 
               colour = "black", size = 0.4, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  stat_contour(data = bathy, aes(x = lon, y = lat, z = depth, alpha = ..level..), 
               colour = "black", size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  geom_polygon(data = south_africa_coast, aes(x = lon, y = lat, group = group), 
               size = 0.0, colour = NA, fill = "grey70") +
  geom_point(data = metaData2, aes(x = lon, y = lat, fill = coast, colour = coast), 
             alpha = 0.9, size = 2.6, shape = 21) +
  geom_point(data = site_pixels3, aes(x = lon, y = lat), 
             shape = 0, alpha = 1.0, size = 2.0, show.legend = FALSE) +
  geom_text(data = metaData2[-c(2:6,13,15),], aes(x = lon, y = lat, label = ID), size = 1.8) +
  geom_text_repel(data = metaData2[c(2:6),], aes(x = lon, y = lat, label = ID), size = 1.8) +
  geom_text(data = metaData2[c(20:21),], aes(x = lon, y = lat, label = ID), size = 1.8) +
  labs(title = NULL, x = NULL, y = NULL) +
  scale_colour_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  scale_fill_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  scale_y_continuous(breaks=seq(-35.0, -27.5, 2.5)) +
  scale_x_continuous(breaks=seq(15, 30, 5)) +
  ### Annotate specific things:
  #annotate("text", label = "Cape \n Agulhas", x = 20.0, y = -35.1, size = 2., colour = "black") +
  #annotate("text", label = "False Bay", x = 18.66, y = -34.55, size = 2.2) +
  #annotate("text", label = "Algoa \n Bay", x = 26.4, y = -34.0, size = 2., colour = "black") +
  #annotate("text", label = "Cape \n Point", x = 18.0, y = -34.4, size = 2, colour = "black") +
  #annotate("text", label = "Hamburg", x = 28.3, y = -33.28611111, size = 2, colour = "black") +
  theme(legend.key = element_blank(),
    legend.key.height = unit(0.4, "cm"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1,0), legend.position = c(0.65, 0.65)) +
  coord_cartesian(xlim = sa_lons, ylim = sa_lats)
f1
ggsave("graph/figures/figure1.pdf", width = 8, height = 4)

#############################################################################
## 2. Figure 2
# The line graphs of each site with largest events

# Load time series
load("prep/SA_coastal_temps.RData")
coastalTemp <- SA_coastal_temps[,c(1,3,4)]
coastalTemp$date <- as.Date(coastalTemp$date)
coastalTemp$type <- "insitu"
load("data/OISSTdaily.Rdata")
OISSTdaily$type <- "OISST"

# Combine and reduce to monthly values
tsALL <- rbind(coastalTemp, OISSTdaily)
tsALL$date <- floor_date(tsALL$date, "month")
tsALL <- ddply(tsALL, .(site, type, date), summarize,
               temp = mean(temp, na.rm = TRUE))

# Reorder sites for plotting
siteOrder <- c("Port Nolloth", "Sea Point", "Hout Bay", "Kommetjie", "Fish Hoek", "Muizenberg", "Gordons Bay", "Hermanus", "Ystervarkpunt", "Mossel Bay", "Knysna", "Tsitsikamma West", "Storms River Mouth", "Tsitsikamma East", "Pollock Beach", "Humewood", "Hamburg", "Eastern Beach", "Orient Beach", "Nahoon Beach", "Sodwana")

# Change site names for plotting
newNames <- c("1. Port Nolloth", "2. Sea Point", "3. Hout Bay", "4. Kommetjie", "5. Fish Hoek", "6. Muizenberg", "7. Gordons Bay", "8. Hermanus", "9. Ystervarkpunt", "10. Mossel Bay", "11. Knysna", "12. Tsitsikamma West", "13. Storms River Mouth", "14. Tsitsikamma East", "15. Pollock Beach", "16. Humewood", "17. Hamburg", "18. Eastern Beach", "19. Orient Beach", "20. Nahoon Beach", "21. Sodwana")

#tsALL2 <- rbind(tsALL, tsALL)
tsALL2 <- tsALL
tsALL2$site <- factor(tsALL2$site, levels = siteOrder)
#tsALL2$event <- rep(c("MHW", "MCS"), each = nrow(tsALL))
#tsALL2$index <- paste(tsALL2$site, tsALL2$date, tsALL2$event, tsALL2$type, sep = "-")
tsALL2$index <- paste(tsALL2$site, tsALL2$date, tsALL2$type, sep = "-")

# Run "proc/results2.R" first to populate the environment with the necessary data...

mhwn$event <- "MHW"
#mhwn$index2 <- paste(mhwn$site, mhwn$month, mhwn$event, mhwn$type, sep = "-")
mhwn$index2 <- paste(mhwn$site, mhwn$month, mhwn$type, sep = "-")

mcsn$event <- "MCS"
#mcsn$index2 <- paste(mcsn$site, mcsn$month, mcsn$event, mcsn$type, sep = "-")
mcsn$index2 <- paste(mcsn$site, mcsn$month, mcsn$type, sep = "-")

mhwnSST$event <- "MHW"
#mhwnSST$index2 <- paste(mhwnSST$site, mhwnSST$month, mhwnSST$event, mhwnSST$type, sep = "-")
mhwnSST$index2 <- paste(mhwnSST$site, mhwnSST$month, mhwnSST$type, sep = "-")

mcsnSST$event <- "MCS"
#mcsnSST$index2 <- paste(mcsnSST$site, mcsnSST$month, mcsnSST$event, mcsnSST$type, sep = "-")
mcsnSST$index2 <- paste(mcsnSST$site, mcsnSST$month, mcsnSST$type, sep = "-")

# Function to extract temperature during events
eventTemp <- function(x){
  events <- data.frame()
  for(i in 1:nrow(x)){#nrow(x)){
    x.1 <- droplevels(x[i,])
    x.1$temp <- tsALL2$temp[tsALL2$index == x.1$index2]
    events <- rbind(events, x.1)
  }
  #events$site <- factor(events$site, levels = siteOrder)
  events <- events[,28:34]
  return(events)
}

mhwn2 <- eventTemp(mhwn)
mhwn2$site <- factor(mhwn2$site, levels = siteOrder)

mcsn2 <- eventTemp(mcsn)
mcsn2$site <- factor(mcsn2$site, levels = siteOrder)

mhwnSST2 <- eventTemp(mhwnSST)
mhwnSST2$site <- factor(mhwnSST2$site, levels = siteOrder)

mcsnSST2 <- eventTemp(mcsnSST)
mcsnSST2$site <- factor(mcsn2$site, levels = siteOrder)

# Rename sites
levels(tsALL2$site) <- newNames
levels(mhwn2$site) <- newNames
levels(mcsn2$site) <- newNames
levels(mhwnSST2$site) <- newNames
levels(mcsnSST2$site) <- newNames

# The figure
f2 <- ggplot(data = tsALL2, aes(x = date, y = temp)) + bw_update +
  geom_line(data = tsALL2[tsALL2$type == "insitu",], colour = "#41b6c4", alpha = 1.0, show.legend = F) +
  geom_line(data = tsALL2[tsALL2$type == "OISST",], #linetype = "dotted",
            colour = "#081d58", alpha = 0.8,  show.legend = F) +
  # insitu MHW
  geom_point(data = mhwn2, aes(x = month, y = temp), shape = 21, alpha = 0.9,
             colour = "black", fill = "#41b6c4", size = 2.7, show.legend = T) +
  geom_text(data = mhwn2, aes(x = month, y = temp, label = index), size = 2.6, colour = "black") +
  # insitu MCS
  geom_point(data = mcsn2, aes(x = month, y = temp), shape = 22, alpha = 0.9,
             colour = "black", fill = "#41b6c4", size = 2.7, show.legend = T) +
  geom_text(data = mcsn2, aes(x = month, y = temp, label = index), size = 2.6, colour = "black") +
  # OISST MHW
  geom_point(data = mhwnSST2, aes(x = month, y = temp), shape = 21, alpha = 0.9,
             colour = "black", fill = "#081d58", size = 2.7, show.legend = T) +
  geom_text(data = mhwnSST2, aes(x = month, y = temp, label = index), colour = "white", size = 2.6) +
  # OISST MCS
  geom_point(data = mcsnSST2, aes(x = month, y = temp), shape = 22, alpha = 0.9,
             colour = "black", fill = "#081d58", size = 2.7, show.legend = T) +
  geom_text(data = mcsnSST2, aes(x = month, y = temp, label = index), colour = "white", size = 2.6) +
  # Additional stuff
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0.015,0)) +
  ylab(expression(paste("Temperature (", degree~C, ")"))) + xlab("Date") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_grid(site ~ .)
f2
ggsave("graph/figures/figure2.pdf", height = 28, width = 12)

#############################################################################
## 3. Figure 3
# The line graph provided by Eric

# No one here but us chickens

#############################################################################
## 4. Figure 4
# OISST coastal images

## Still need OISST data from AJ

#############################################################################
## 5. Figure 5
# Co-occurrence figure
mhwCoastCO <- read.csv("data/mhwCoastCO.csv")
mhwCoastCO$event <- "MHW"
mcsCoastCO <- read.csv("data/mcsCoastCO.csv")
mcsCoastCO$event <- "MCS"
allCoastCO <- rbind(mhwCoastCO, mcsCoastCO)
allCoastCO$index <- paste(allCoastCO$lag, allCoastCO$direction, sep = "_")
allCoastCO$coast <- factor(allCoastCO$coast, levels = c("west", "south", "east"))
allCoastCO$index2 <- paste(allCoastCO$event, allCoastCO$coast, sep = ": ")
allCoastCO$index2 <- factor(allCoastCO$index2, levels = levels(as.factor(allCoastCO$index2))[c(3,2,1,6,5,4)])
allCoastCO$direction <- factor(allCoastCO$direction, levels = c("b", "x", "a"))
levels(allCoastCO$direction) <- c("before", "combined", "after")

cols <- c("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")

figure5 <-  ggplot(data = allCoastCO, aes(x = quantile, y = proportion, group = index)) + bw_update +
  geom_line(aes(colour = as.factor(lag))) + 
  # geom_point(aes(colour = as.factor(lag)), size = 0.1) +
  facet_grid(index2 ~ direction, space = "free_y", shrink = T) +
  scale_color_grey() +
  ylab("proportion") + xlab("percentile (%)") +
  scale_colour_manual(values = rev(cols)) +
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2)) +
  guides(colour = guide_legend("lag [days]", nrow = 1, byrow = T, override.aes = list(size = 1.5))) +
  theme(axis.text = element_text(size = 7),
        legend.position = "bottom")
figure5
ggsave("graph/figures/figure5.pdf", width = 6, height = 7)
