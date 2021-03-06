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
require(ggplot2); require(lubridate); require(plyr); require(ggrepel); require(RColorBrewer)
library(doMC); doMC::registerDoMC(4)
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
geom_contour(data = bathy[bathy$depth >= -250,], aes(x = lon, y = lat, z = depth),
             colour = "black", alpha = 0.85, size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
stat_contour(data = bathy[bathy$depth < -250,], aes(x = lon, y = lat, z = depth, alpha = ..level..),
             colour = "black", size = 0.1, binwidth = 1000, na.rm = TRUE, show.legend = FALSE) +
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
  theme(legend.key.height = unit(0.4, "cm"),
    legend.background = element_blank(),
    legend.direction = "horizontal",
    legend.justification = c(1,0), legend.position = c(0.65, 0.65)) +
  coord_cartesian(xlim = sa_lons, ylim = sa_lats)
f1
ggsave("LaTeX/figure1.pdf", width = 8, height = 4)

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
               temp = mean(temp, na.rm = TRUE), .parallel = T)

# Reorder sites for plotting
siteOrder <- c("Port Nolloth", "Sea Point", "Hout Bay", "Kommetjie", "Fish Hoek", "Muizenberg", "Gordons Bay", "Hermanus", "Ystervarkpunt", "Mossel Bay", "Knysna", "Tsitsikamma West", "Storms River Mouth", "Tsitsikamma East", "Pollock Beach", "Humewood", "Hamburg", "Eastern Beach", "Orient Beach", "Nahoon Beach", "Sodwana")

# Change site names for plotting on the side of the facets
newNames1 <- c("1. Port \n Nolloth", "2. Sea \n Point", "3. Hout \n Bay", "4. Kommetjie", "5. Fish \n Hoek", "6. Muizenberg", "7. Gordons \n Bay", "8. Hermanus", "9. Yster- \n varkpunt", "10. Mossel \n Bay", "11. Knysna", "12. Tsitsikamma \n West", "13. Storms River \n Mouth", "14. Tsitsikamma \n East", "15. Pollock \n Beach", "16. Humewood", "17. Hamburg", "18. Eastern \n Beach", "19. Orient \n Beach", "20. Nahoon \n Beach", "21. Sodwana")
# Change site names for plotting on the top of the facets
newNames2 <- c("1. PortNolloth", "2. Sea Point", "3. Hout Bay", "4. Kommetjie", "5. Fish Hoek", "6. Muizenberg", "7. Gordons Bay", "8. Hermanus", "9. Yster varkpunt", "10. Mossel Bay", "11. Knysna", "12. Tsitsikamma West", "13. Storms River Mouth", "14. Tsitsikamma East", "15. Pollock Beach", "16. Humewood", "17. Hamburg", "18. Eastern Beach", "19. Orient Beach", "20. Nahoon Beach", "21. Sodwana")
# Change site names for plotting only numbers
newNames3 <- c("1. ", "2. ", "3. ", "4. ", "5. ", "6. ", "7. ", "8. ", "9. ", "10. ", "11. ", "12. ", "13. ", "14. ", "15. ", "16. ", "17. ", "18. ", "19. ", "20. ", "21. ")

## Create second data.frame for plotting and first add coast labels
tsALL2 <- data.frame()
for(i in 1:length(levels(tsALL$site))) {
  x <- subset(tsALL, site == levels(tsALL$site)[i])
  if(x$site[1] %in% wc) {
    x$coast <- "WC"
  } else if(x$site[1] %in% sc) {
    x$coast <- "SC"
  } else if(x$site[1] %in% ec) {
    x$coast <- "EC"
  }
  tsALL2 <- rbind(tsALL2, x)
}

# Reorder sites
tsALL2$site <- factor(tsALL2$site, levels = siteOrder)
tsALL2$index <- paste(tsALL2$site, tsALL2$date, tsALL2$type, sep = "-")

### Run "proc/results2.R" first to populate the environment with the necessary data... ###

# Change coast labels for plotting and add x y coords
tsALL2$coast <- as.factor(tsALL2$coast)
tsALL2$x <- as.Date("1994-03-01"); tsALL2$y <- 18.5

mhwn$event <- "MHW"
mhwn$index2 <- paste(mhwn$site, mhwn$month, mhwn$type, sep = "-")

mcsn$event <- "MCS"
mcsn$index2 <- paste(mcsn$site, mcsn$month, mcsn$type, sep = "-")

mhwnSST$event <- "MHW"
mhwnSST$index2 <- paste(mhwnSST$site, mhwnSST$month, mhwnSST$type, sep = "-")

mcsnSST$event <- "MCS"
mcsnSST$index2 <- paste(mcsnSST$site, mcsnSST$month, mcsnSST$type, sep = "-")

# Function to extract temperature during events
eventTemp <- function(x){
  events <- data.frame()
  for(i in 1:nrow(x)){#nrow(x)){
    x.1 <- droplevels(x[i,])
    x.1$temp <- tsALL2$temp[tsALL2$index == x.1$index2]
    events <- rbind(events, x.1)
  }
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
levels(tsALL2$site) <- newNames2
levels(mhwn2$site) <- newNames2
levels(mcsn2$site) <- newNames2
levels(mhwnSST2$site) <- newNames2
levels(mcsnSST2$site) <- newNames2

# Add x y coords for plotting R2 values
resultsR2$site <- factor(resultsR2$site, levels = siteOrder)
levels(resultsR2$site) <- newNames2
resultsR2$x <- as.Date("1976-01-01"); resultsR2$y <- 25
# Manually change row order as other methods don't work...
resultsR2 <- resultsR2[c(15,16,6,9,2,11,3,5,21,10,8,20,18,19,14,7,4,1,13,12,17),]

## Set the colour palette to be used
# colour
plot_colours <- c("#41b6c4", "#081d58")
# greyscale
# plot_colours <- c("grey50", "black")

## The figure
f2 <- ggplot(data = tsALL2, aes(x = date, y = temp)) + bw_update +
  geom_text(data = tsALL2, aes(x = x, y = y, label = coast), size = 12, colour = "grey80") +
  geom_text(data = resultsR2, aes(x = x, y = y, label = as.character(resultsR2$R22), group = site), parse = TRUE, size = 2.6) +
  ## Colour version
  geom_line(data = tsALL2[tsALL2$type == "insitu",], colour = plot_colours[1], alpha = 1.0, show.legend = F) +
  geom_line(data = tsALL2[tsALL2$type == "OISST",], #linetype = "dotted",
            colour = plot_colours[2], alpha = 0.8,  show.legend = F) +
  # insitu MHW
  geom_point(data = mhwn2, aes(x = month, y = temp), shape = 21, alpha = 0.9,
             colour = "black", fill = plot_colours[1], size = 2.7, show.legend = T) +
  geom_text(data = mhwn2, aes(x = month, y = temp, label = index), size = 2.6, colour = "black") +
  # insitu MCS
  geom_point(data = mcsn2, aes(x = month, y = temp), shape = 22, alpha = 0.9,
             colour = "black", fill = plot_colours[1], size = 2.7, show.legend = T) +
  geom_text(data = mcsn2, aes(x = month, y = temp, label = index), size = 2.6, colour = "black") +
  # OISST MHW
  geom_point(data = mhwnSST2, aes(x = month, y = temp), shape = 21, alpha = 0.9,
             colour = "black", fill = plot_colours[2], size = 2.7, show.legend = T) +
  geom_text(data = mhwnSST2, aes(x = month, y = temp, label = index), colour = "white", size = 2.6) +
  # OISST MCS
  geom_point(data = mcsnSST2, aes(x = month, y = temp), shape = 22, alpha = 0.9,
             colour = "black", fill = plot_colours[2], size = 2.7, show.legend = T) +
  geom_text(data = mcsnSST2, aes(x = month, y = temp, label = index), colour = "white", size = 2.6) +
  # Additional stuff
  scale_y_continuous(breaks = c(15,25)) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y", expand = c(0.015,0)) +
  ylab(expression(paste("Temperature (", degree~C, ")"))) + xlab("Date") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  facet_wrap(~site, ncol = 3)
f2
ggsave("LaTeX/figure2.pdf", height = 6, width = 10) # Use this line when saving a colour version
# ggsave("LaTeX/figure2_bw.pdf", height = 6, width = 10) # And this line for greyscale

#############################################################################
## 3. Figure 3
# The line graph provided by Eric

# No one here but us chickens

#############################################################################
## 4. Figure 4
# Co-occurrence MHW figure
# Run "proc/results2.R" first to populate the environment with the necessary data...

# The red colours to use
colsr <- c("#ffce39ff", "#feaa18ff", "#fe8418ff", "#f85302ff", "#f72703ff", "#b61a1cff", "#880027ff")

# Add index for plotting
mhwCoastCO$index <- paste(mhwCoastCO$lag, mhwCoastCO$direction, sep = "_")

# Correct order of windows
mhwCoastCO$direction <- factor(mhwCoastCO$direction, levels = c("b", "x", "a"))
levels(mhwCoastCO$direction) <- c("before", "combined", "after")

# Correct site order
mhwCoastCO$site <- factor(mhwCoastCO$site, levels = siteOrder)

# Rename sites
levels(mhwCoastCO$site) <- newNames3

# Add label column
mhwCoastCO$label <- paste(mhwCoastCO$site, mhwCoastCO$direction, sep = "")

# Create correct order for labels for plotting... not glamorous
labelOrder <- c("1. before", "1. combined", "1. after", "2. before", "2. combined", "2. after", "3. before", "3. combined", "3. after", "4. before", "4. combined", "4. after", "5. before", "5. combined", "5. after", "6. before", "6. combined", "6. after", "7. before", "7. combined", "7. after", "8. before", "8. combined", "8. after", "9. before", "9. combined", "9. after", "10. before", "10. combined", "10. after", "11. before", "11. combined", "11. after", "12. before", "12. combined", "12. after", "13. before", "13. combined", "13. after", "14. before", "14. combined", "14. after", "15. before", "15. combined", "15. after", "16. before", "16. combined", "16. after", "17. before", "17. combined", "17. after", "18. before", "18. combined", "18. after", "19. before", "19. combined", "19. after", "20. before", "20. combined", "20. after", "21. before", "21. combined", "21. after")

# Correct label order
mhwCoastCO$label <- factor(mhwCoastCO$label, levels = labelOrder)

# Correct lag order
# mhwCoastCO$lag <- as.factor(mhwCoastCO$lag)
# mhwCoastCO$lag <- factor(mhwCoastCO$lag, levels = c("2","4","6","8","10","12","14"))

# Remove combined columns
mhwCoastCO2 <- mhwCoastCO[mhwCoastCO$direction != "combined",]

# Change coast labels for plotting and add x y coords
mhwCoastCO2$coast <- as.factor(mhwCoastCO2$coast)
levels(mhwCoastCO2$coast) <- c("EC", "SC", "WC")
mhwCoastCO2$x <- 0.5; mhwCoastCO2$y <- 0.5

figure4 <-  ggplot(data = mhwCoastCO2, aes(x = quantile, y = proportion, group = index)) + bw_update +
  geom_text(aes(x, y, label = coast), size = 12.0, colour = "grey80") +
  geom_line(aes(colour = as.factor(lag))) +
  # geom_point(aes(colour = as.factor(lag)), size = 0.1) +
  facet_wrap(~ label, ncol = 6) +
  #label_both(labels = mhwCoastCO$label) +
  #scale_color_grey() +
  ylab("proportion") + xlab("percentile (%)") +
  #scale_colour_brewer(palette = "YlOrRd", direction = -1) +
  scale_colour_manual(values = rev(colsr)) + # Colour
  # scale_color_grey() + # greyscale
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2)) +
  scale_y_continuous(breaks = c(0.25, 0.50, 0.75), limits = c(0,1)) +
  guides(colour = guide_legend("lag [days]", nrow = 1, byrow = T, override.aes = list(size = 1.5))) +
  theme(axis.text = element_text(size = 7),
        legend.position = "bottom")
figure4
ggsave("LaTeX/figure4.pdf", height = 8, width = 6)
# ggsave("LaTeX/figure4_bw.pdf", height = 8, width = 6) # Use this line when generating a greyscale figure

#############################################################################
## 5. Figure 5
# Co-occurrence MCS figure
# Run "proc/results2.R" first to populate the environment with the necessary data...

# The blue colours to use
colsb <- c("#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")

# Add index for plotting
mcsCoastCO$index <- paste(mcsCoastCO$lag, mcsCoastCO$direction, sep = "_")

# Correct order of windows
mcsCoastCO$direction <- factor(mcsCoastCO$direction, levels = c("b", "x", "a"))
levels(mcsCoastCO$direction) <- c("before", "combined", "after")

# Correct site order
mcsCoastCO$site <- factor(mcsCoastCO$site, levels = siteOrder)

# Rename sites
levels(mcsCoastCO$site) <- newNames3

# Add label column
mcsCoastCO$label <- paste(mcsCoastCO$site, mcsCoastCO$direction, sep = "")

# Correct label order
mcsCoastCO$label <- factor(mcsCoastCO$label, levels = labelOrder)

# Remove combined columns
mcsCoastCO2 <- mcsCoastCO[mcsCoastCO$direction != "combined",]

# Change coast labels for plotting and add x y coords
mcsCoastCO2$coast <- as.factor(mcsCoastCO2$coast)
levels(mcsCoastCO2$coast) <- c("EC", "SC", "WC")
mcsCoastCO2$x <- 0.5; mcsCoastCO2$y <- 0.5

figure5 <-  ggplot(data = mcsCoastCO2, aes(x = quantile, y = proportion, group = index)) + bw_update +
  geom_text(aes(x, y, label = coast), size = 12.0, colour = "grey80") +
  geom_line(aes(colour = as.factor(lag))) +
  # geom_point(aes(colour = as.factor(lag)), size = 0.1) +
  facet_wrap(~ label, ncol = 6) +
  #scale_color_grey() +
  ylab("proportion") + xlab("percentile (%)") +
  scale_colour_manual(values = rev(colsb)) + # Colour
  # scale_color_grey() + # greyscale
  scale_x_continuous(breaks = seq(0.0, 1.0, 0.2)) +
  scale_y_continuous(breaks = c(0.25, 0.50, 0.75), limits = c(0,1)) +
  guides(colour = guide_legend("lag [days]", nrow = 1, byrow = T, override.aes = list(size = 1.5))) +
  theme(axis.text = element_text(size = 7),
        legend.position = "bottom")
figure5
ggsave("LaTeX/figure5.pdf", height = 8, width = 6)
# ggsave("LaTeX/figure5_bw.pdf", height = 8, width = 6)

