#############################################################################
## This script does:
# 1. Reads in MUR netCDF SST data;
# 2. Creates plotting parameters;
# 3. Plots temperatures and creates multiple panels;
# 4. Combines panels and saves as "graph/figure1_1km_inset_map.pdf"
#############################################################################

#############################################################################
## DEPENDS ON:
require(ncdf); require(colorRamps); library(reshape2); library(ggplot2); require(grid)
require(rgeos); require(maptools) # " maptools" must be loaded after "rgeos"
require(lubridate); require(plyr); require(ggrepel); library(viridis)
source("setupParams/theme.R"); source("proc/results2.R")

#############################################################################
## USED BY:
# Nothing
#############################################################################

#############################################################################
## CREATES:
# "graph/figure1_1km_inset_map.pdf"
#############################################################################

setwd("/Users/ajsmit/Dropbox/repos/MHW")

################################################################################
# Open a netCDF file
# ncsst <- open.ncdf("/Volumes/AGULHAS/OceanData/misc_hires/jplG1SST_2f50_4900_a129.nc")
ncsst <- open.ncdf("/Users/ajsmit/spatial/G1SST/jplG1SST_2016-02-01--2016-02-14.nc")
str(ncsst$dim)
# Get the sst etc.
sst <- get.var.ncdf(ncsst, "SST", start = c(1,1,10), count = c(-1,-1,1))
sst[sst == -1] <- NA

# Some basic info of the lats and lons
lat_min <- min(ncsst$dim$lat$vals); lat_max <- max(ncsst$dim$lat$vals)
lon_min <- min(ncsst$dim$lon$vals); lon_max <- max(ncsst$dim$lon$vals)
lat <- ncsst$dim$lat$vals
lon <- ncsst$dim$lon$vals

################################################################################
## Prepare graphing variables and environment
grid <- expand.grid(lon, lat) # Matches each possible lat to each possible lon in the sequences
msst <- melt(sst) # Melts SST lat and lon and temp into columns

rlat <- c(-45, -24); rlon <- c(13, 37)

# Load SA map data
load("graph/south_africa_coast.RData") # Lowres
names(south_africa_coast)[1] <- "lon"
load("graph/sa_shore.Rdata") # Hires
names(sa_shore)[4:5] <- c("lon","lat")

# Load meatadata and site pixels
#==============================================================================
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
fb_lats <- c(-34.4, -33.85); fb_lons <- c(18.0, 18.9)
ham_lats <- c(-33.4, -32.8); ham_lons <- c(27.3, 28.1)

# Load SA bathymetry
load("data/bathy/bathy.RData") # HiRes for final image
load("data/bathy/sa_bathy.RData") # LowRes for tweaking

################################################################################
## Assemble the graph

theme_set(theme_bw())

limits <- c(12,28) # for colour bar
breaks <- seq(6, 30, 2) # Create breaks to be used for colour bar


# The main map
#==============================================================================
sa <- ggplot(data = south_africa_coast, aes(x = lon, y = lat)) + bw_update +
  geom_raster(data = grid, aes(x = Var1, y = Var2, fill = msst$value)) +
  geom_contour(data = bathy[bathy$depth >= -250,], aes(x = lon, y = lat, z = depth),
               colour = "grey20", alpha = 0.7, size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  stat_contour(data = sa_bathy[sa_bathy$depth < -250,], aes(x = lon, y = lat, z = depth, alpha = ..level..),
               colour = "grey20", size = 0.1, binwidth = 1000, na.rm = TRUE, show.legend = FALSE) +
  geom_polygon(data = south_africa_coast, aes(x = lon, y = lat, group = group),
               fill = "#929292", colour = "#929292", size = 0.1, show.legend = FALSE) +
  # geom_point(data = metaData2, aes(x = lon, y = lat, colour = coast),
  #            alpha = 0.9, size = 2.6, shape = 1) +
  geom_point(data = metaData2, aes(x = lon, y = lat),
             alpha = 0.9, size = 2.6, shape = 1, colour = "ivory") +
  geom_point(data = site_pixels3, aes(x = lon, y = lat),
             shape = 0, alpha = 1.0, size = 1.2, show.legend = FALSE) +
  geom_text(data = metaData2[-c(2:6,13,15,18,19),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "ivory1") +
  # geom_text_repel(data = metaData2[c(2:6),], aes(x = lon, y = lat, label = ID),
  #                 size = 1.8, colour = "ivory1", segment.color = "ivory1") +
  geom_text(data = metaData2[c(20:21),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "ivory1") +
  coord_equal() +
  scale_x_continuous(limits = sa_lons, expand = c(0, 0), breaks = seq(15, 35, 5)) +
  scale_y_continuous(limits = sa_lats, expand = c(0, 0), breaks = seq(-35, -30, 5)) +
  scale_fill_viridis(breaks = breaks, limits = limits, expression(paste("Temp. (",degree,"C)"))) +
  # scale_colour_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colourbar(barheight = 1.00, barwidth = 10)) +
  theme(panel.background = element_rect(fill = "ivory", colour = NA),
        panel.border = element_rect(colour = "black", size = 0.5),
        panel.grid.minor = element_line(colour = "NA"),
        panel.grid.major = element_line(colour = "ivory", size = 0.2, linetype = "dotted"),
        legend.direction = "horizontal",
        legend.justification = c(1,0),
        legend.position = c(0.65, 0.80),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.5))
sa
#ggsave("graph/figures/figure1_1km.pdf", width = 8, height = 4)

# False Bay inset
#==============================================================================
fb <- ggplot(data = sa_shore, aes(x = lon, y = lat)) + bw_update +
  #geom_raster(data = grid, aes(x = Var1, y = Var2, fill = msst$value)) +
  geom_contour(data = bathy[bathy$depth >= -250,], aes(x = lon, y = lat, z = depth),
               colour = "grey20", alpha = 0.7, size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  # stat_contour(data = sa_bathy[sa_bathy$depth < -250,], aes(x = lon, y = lat, z = depth, alpha = ..level..),
  #              colour = "grey20", size = 0.1, binwidth = 1000, na.rm = TRUE, show.legend = FALSE) +
  geom_polygon(data = sa_shore, aes(x = lon, y = lat, group = PID),
               fill = "#929292", colour = "#929292", size = 0.1, show.legend = FALSE) +
  # geom_point(data = metaData2, aes(x = lon, y = lat, colour = coast),
  #            alpha = 0.9, size = 2.6, shape = 1) +
  geom_point(data = metaData2[c(2:6),], aes(x = lon, y = lat),
             alpha = 0.9, size = 2.6, shape = 1, colour = "black") +
  # geom_point(data = site_pixels3, aes(x = lon, y = lat),
  #            shape = 0, alpha = 1.0, size = 1.2, show.legend = FALSE) +
  geom_text(data = metaData2[c(2:6),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "black") +
  coord_equal() +
  coord_map(xlim = fb_lons, ylim = fb_lats, projection = "mercator") + 
  scale_x_continuous(breaks = c(18.2, 18.7)) +
  scale_y_continuous(breaks = c(-34.0, -34.3)) +
  # scale_fill_viridis(breaks = breaks, limits = limits, expression(paste("Temp. (",degree,"C)"))) +
  # scale_colour_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colourbar(barheight = 1.00, barwidth = 10)) +
  theme(panel.background = element_rect(fill = "ivory", colour = NA),
        panel.border = element_rect(colour = "ivory", size = 0.5),
        panel.grid.minor = element_line(colour = "NA"),
        panel.grid.major = element_line(colour = "ivory", size = 0.2, linetype = "dotted"),
        legend.direction = "horizontal",
        legend.justification = c(1,0),
        legend.position = c(0.65, 0.65),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.text = element_text(colour = "ivory"),
        axis.ticks = element_line(size = 0.5, colour = "ivory"))
fb

# Hamburg inset
#==============================================================================
ham <- ggplot(data = sa_shore, aes(x = lon, y = lat)) + bw_update +
  #geom_raster(data = grid, aes(x = Var1, y = Var2, fill = msst$value)) +
  # geom_contour(data = bathy[bathy$depth >= -250,], aes(x = lon, y = lat, z = depth),
  #              colour = "grey20", alpha = 0.7, size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  stat_contour(data = sa_bathy[sa_bathy$depth < -250,], aes(x = lon, y = lat, z = depth, alpha = ..level..),
               colour = "grey20", size = 0.1, binwidth = 1000, na.rm = TRUE, show.legend = FALSE) +
  geom_polygon(data = sa_shore, aes(x = lon, y = lat, group = PID),
               fill = "#929292", colour = "#929292", size = 0.1, show.legend = FALSE) +
  # geom_point(data = metaData2, aes(x = lon, y = lat, colour = coast),
  #            alpha = 0.9, size = 2.6, shape = 1) +
  geom_point(data = metaData2[c(17:20),], aes(x = lon, y = lat),
             alpha = 0.9, size = 2.6, shape = 1, colour = "black") +
  # geom_point(data = site_pixels3, aes(x = lon, y = lat),
  #            shape = 0, alpha = 1.0, size = 1.2, show.legend = FALSE) +
  geom_text(data = metaData2[c(17:20),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "black") +
  coord_equal() +
  coord_map(xlim = ham_lons, ylim = ham_lats, projection = "mercator") + 
  scale_x_continuous(breaks = c(27.5, 27.9)) +
  scale_y_continuous(breaks = c(-33.3, -33.0)) +
  # scale_fill_viridis(breaks = breaks, limits = limits, expression(paste("Temp. (",degree,"C)"))) +
  # scale_colour_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colourbar(barheight = 1.00, barwidth = 10)) +
  theme(panel.background = element_rect(fill = "ivory", colour = NA),
        panel.border = element_rect(colour = "ivory", size = 0.5),
        panel.grid.minor = element_line(colour = "NA"),
        panel.grid.major = element_line(colour = "ivory", size = 0.2, linetype = "dotted"),
        legend.direction = "horizontal",
        legend.justification = c(1,0),
        legend.position = c(0.65, 0.65),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.text = element_text(colour = "ivory"),
        axis.ticks = element_line(size = 0.5, colour = "ivory"))
ham

# Compile and save
#==============================================================================
pdf("graph/figure1_1km_inset_map.pdf", width = 8, height = 5, pointsize = 10) # Set PDF dimensions
vp1 <- viewport(x = 1.0, y = 1.0, w = 1.0, h = 1.0, just = c("right", "top")) # South Africa
vp2 <- viewport(x = 0.47, y = 0.64, w = 0.25, h = 0.25, just = c("right", "top"))  # False Bay
vp3 <- viewport(x = 0.76, y = 0.7, w = 0.25, h = 0.25, just = c("right", "top"))  # Hamburg
print(sa, vp = vp1)
print(fb, vp = vp2)
print(ham, vp = vp3)
dev.off()
