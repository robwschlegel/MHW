#############################################################################
## This script does:
# 1. Reads in MUR netCDF SST data;
# 2. Creates plotting parameters;
# 3. Plots temperatures;
# 4. Saves as "graph/Figure_5_MUR_v3.3.pdf"
#############################################################################

#############################################################################
## DEPENDS ON:
rotate <- function(x) t(apply(x, 2, rev))
require(ncdf); require(colorRamps); library(reshape2); library(ggplot2)
require(rgeos); require(maptools) # " maptools" must be loaded after "rgeos"
require(lubridate); require(plyr); require(ggrepel); library(viridis)
source("setupParams/theme.R"); source("proc/results2.R")

#############################################################################
## USED BY:
# Nothing
#############################################################################

#############################################################################
## CREATES:
#
#############################################################################

setwd("/Users/ajsmit/Dropbox/repos/MHW")

################################################################################
# Open a netCDF file
str(ncsst$dim)
ncsst <- open.ncdf("/Volumes/AGULHAS/OceanData/misc_hires/jplG1SST_2f50_4900_a129.nc")

# Get the sst etc.
sst <- get.var.ncdf(ncsst, "SST", start = c(1,1,1), count = c(-1,-1,1))
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

if(file.exists("C:/Users/Robert Schlegel/Documents/R/gshhg-bin-2.3.0/gshhs_f.b") == TRUE){
  shore <- getRgshhsMap("C:/Users/Robert Schlegel/Documents/R/gshhg-bin-2.3.0/gshhs_f.b",
                        xlim = rlon, ylim = rlat, level = 1, no.clip = FALSE, checkPolygons = TRUE)} else {
                          shore <- getRgshhsMap("/Users/ajsmit/spatial/gshhs_2.3.4/gshhg-bin-2.3.4/gshhs_f.b",
                                                xlim = rlon,
                                                ylim = rlat,
                                                level = 1, no.clip = FALSE,
                                                checkPolygons = TRUE)}

shore2 <- fortify(shore) # Convert shore polygon into a data.frame

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

# Load SA bathymetry
load("data/bathy/bathy.RData") # HiRes for final image
load("data/bathy/sa_bathy.RData") # LowRes for tweaking

################################################################################
## Assemble the graph

theme_set(theme_bw())

limits <- c(6,28) # for colour bar
breaks <- seq(6, 30.2, 2) # Create breaks to be used for colour bar

p <- ggplot(data = shore2, aes(x = long, y = lat)) +
  geom_raster(data = grid, aes(x = Var1, y = Var2, fill = msst$value)) +
  geom_contour(data = bathy[bathy$depth >= -250,], aes(x = lon, y = lat, z = depth),
               colour = "grey90", alpha = 0.7, size = 0.2, binwidth = 200, na.rm = TRUE, show.legend = FALSE) +
  stat_contour(data = bathy[bathy$depth < -250,], aes(x = lon, y = lat, z = depth, alpha = ..level..),
               colour = "grey90", size = 0.1, binwidth = 1000, na.rm = TRUE, show.legend = FALSE) +
  geom_polygon(data = shore2, aes(x = long, y = lat, group = group),
               fill = "#929292", colour = "black", size = 0.1, show.legend = FALSE) +
  # geom_point(data = metaData2, aes(x = lon, y = lat, colour = coast),
  #            alpha = 0.9, size = 2.6, shape = 1) +
  geom_point(data = metaData2, aes(x = lon, y = lat),
             alpha = 0.9, size = 2.6, shape = 1, colour = "white") +
  geom_point(data = site_pixels3, aes(x = lon, y = lat),
             shape = 0, alpha = 1.0, size = 1.2, show.legend = FALSE) +
  geom_text(data = metaData2[-c(2:6,13,15),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "ivory1") +
  geom_text_repel(data = metaData2[c(2:6),], aes(x = lon, y = lat, label = ID),
                  size = 1.8, colour = "ivory1", segment.color = "ivory1") +
  geom_text(data = metaData2[c(20:21),], aes(x = lon, y = lat, label = ID),
            size = 1.8, colour = "ivory1") +
  coord_equal() +
  scale_x_continuous(limits = rlon, expand = c(0, 0), breaks = seq(15, 35, 5)) +
  scale_y_continuous(limits = rlat, expand = c(0, 0)) +
  scale_fill_viridis(breaks = breaks, limits = limits, expression(paste("Temp. (",degree,"C)"))) +
  # scale_colour_manual(breaks = c("west", "south", "east"), values = c("#8dd3c7", "#4daf4a", "#e41a1c")) +
  xlab("") +
  ylab("") +
  guides(fill = guide_colourbar(barheight = 8.35)) +
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(colour = "black", size = 0.5),
        panel.grid.minor = element_line(colour = "NA"),
        panel.grid.major = element_line(colour = "grey20", size = 0.2, linetype = "dotted"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.key = element_blank(),
        legend.background = element_blank(),
        axis.ticks = element_line(size = 0.5))
# p
ggsave("graph/figures/figure1_1km.pdf", width = 8, height = 4)