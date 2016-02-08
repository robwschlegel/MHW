#############################################################################
## This script does:
# 1. Loads metadata for used sites;
# 2. Calculates shore normal transects, pixel points and bounding boxes;
# 4. Creates coastal bounding boxes/ polygon based on 200m isobath
# 5. Creates graph of bounding boxes with in situ points and isobath ranges (aw yisss)
#############################################################################

#############################################################################
## DEPENDS ON:
require(plyr); require(zoo); require(lubridate); require(reshape2); require(maptools); require(sp); require(geosphere); require(marmap); require(PBSmapping); require(mapproj); require(ggplot2); require(gridExtra); require(grid); require(viridis); require(scales)
source("func/shoreNormalTransectFunc.R"); source("setupParams/theme.R")
# "metaData2.csv"
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
# "graph/boundingBoxesAll.pdf"
#############################################################################

#############################################################################
## Read in metaData
load("data/metaData2.Rdata")
# Remove coastal mean data
metaData2 <- metaData2[-((length(metaData2$site)-2):length(metaData2$site)),]

#############################################################################
## Calculate shore normal transects, pixel points and bounding boxes

# Transects
site_transects <- data.frame()
for(i in 1:length(metaData2$site)){
 site <- metaData2[i,]
 site_transect <- shore.normal.transect(site, 2)
 site_transects <- rbind(site_transects, site_transect)
}

# Manually correct some transects
site_transects$heading[2:4] <- 290
site_transects$heading[5:6] <- 178
site_transects$heading[13:14] <- 177.9855
save(site_transects, file = "data/site_transects.RData")
write.csv(site_transects, file = "data/site_transects.csv", row.names = F)
load("data/site_transects.RData")

# Pixel points
site_pixels <- data.frame()
for(i in 1:length(metaData2$site)){
  site <- site_transects[i,]
  site_pixel <- transect.pixel.isobath(site, 25000, -200)
  site_pixels <- rbind(site_pixels, site_pixel)
}
save(site_pixels, file = "data/site_pixels.RData")
write.csv(site_pixels, file = "data/site_pixels.csv", row.names = F)
load("data/site_pixels.RData")

# Bounding box
  # Only one is made in order to know how large the the geom_point() squares should be made to match
test <- data.frame(xmin = destPoint(p = site_pixels[1,2:3], b = 270, d = 12500)[1],
                   xmax = destPoint(p = site_pixels[1,2:3], b = 90, d = 12500)[1],
                   ymin = destPoint(p = site_pixels[1,2:3], b = 180, d = 12500)[2],
                   ymax = destPoint(p = site_pixels[1,2:3], b = 0, d = 12500)[2])

#############################################################################
## Create 200m isobath polygon for plotting

# SA coastline
load("graph/south_africa_coast.RData")
names(south_africa_coast)[1] <- "lon"
south_africa_coast$site <- "SA"

# Manually divide up coastline
wc <- south_africa_coast[291:410,]
sc <- south_africa_coast[132:291,]
ec <- south_africa_coast[23:132,]

# Function for calculating bounding boxes
boundingBoxIso <- function(dat, distance = 400000, isobath = -200){
  #df <- data.frame()
  for(i in rev(1:length(dat$lon))){
    distances <- seq(from = 0, to = distance, by = 1000)
    heading <- shore.normal.transect(dat[i,])
    coords <- coords <- data.frame(lon = heading$lon, lat = heading$lat)
    distances2 <- as.data.frame(destPoint(p = coords, b = heading$heading, d = distances))
    sitesIdx <- knnx.index(sa_bathy[,1:2], as.matrix(distances2), k = 1)
    bathy2 <- sa_bathy[sitesIdx,]
    bathy2 <- bathy2[complete.cases(bathy2[,3]),]
    bathy2 <- bathy2[bathy2$depth >= isobath,]
    bathy2 <- bathy2[length(bathy2$depth),]
    if(nrow(bathy2) < 1){
      dat2 <- data.frame(dat[i, ])
    }else{
      dat2 <- cbind(bathy2[,1:2], dat[i,3:8])
    }
    dat <- rbind(dat, dat2)
  }
  return(dat)
}

wcBox <- boundingBoxIso(wc)
wcBox <- wcBox[-c(140:143),] # Remove island
scBox <- boundingBoxIso(sc)
ecBox <- boundingBoxIso(ec)


# isobath200 <- sa_bathy[sa_bathy$depth >= -200,]
# isobath200$coast <- NA
# isobath200$coast[1:1500] <- "west"
# isobath200$coast[1501:3200] <- "south"
# isobath200$coast[3201:4448] <- "east"

#############################################################################
## Create figure showing all the stuffs

# Define plotting parameters
sa_lats <- c(-37, -27); sa_lons <- c(14, 34)

# The Figure
f1 <- ggplot() + theme_bw() + #coord_equal() + 
  geom_raster(data = sa_bathy, aes(x = lon, y = lat, fill = depth)) +
  #geom_raster(data = isobath200, aes(x = lon, y = lat, fill = coast)) +
  stat_contour(data = sa_bathy, aes(x = lon, y = lat, z = depth, alpha = ..level..), 
               colour = "black", size = 0.2, binwidth = 200, na.rm = TRUE, show_guide = FALSE) +
  geom_polygon(data = south_africa_coast, aes(x = lon, y = lat, group = group), 
               size = 0.1, colour = "black", fill = "grey80") +
  #geom_path(data = sa_provinces_new, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = scBox, aes(x = lon, y = lat, group = group), alpha = 0.20, 
               colour = "green", fill = "green") +
  geom_polygon(data = wcBox, aes(x = lon, y = lat, group = group), alpha = 0.20, 
               colour = "blue", fill = "blue") +
  geom_polygon(data = ecBox, aes(x = lon, y = lat, group = group), alpha = 0.20, 
               colour = "orange", fill = "orange") +
  geom_point(data = metaData2, aes(x = lon, y = lat, colour = coast), alpha = 0.8, size = 1.6) +
  geom_point(data = site_pixels, aes(x = lon, y = lat), colour = "red", shape = 0, alpha = 0.8, size = 2.1) +
  geom_rect(data = test, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            alpha = 1, colour = "red", size = 0.1, linetype = 1) +
  scale_alpha_continuous(breaks = c(-200, -1000, -2000, -3000, -4000, -5000),
                         guide_legend(title = "depth (m)")) +
  scale_fill_gradient(low = "steelblue4", high = "steelblue1", na.value = "steelblue4", 
                      breaks = c(-1000, -2000, -3000, -4000, -5000),
                      guide_legend(title = "depth (m)")) +
  labs(title = NULL, x = NULL, y = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  #guide_legend(ncol = 2) +
  ### Annotate specific things:
  annotate("text", label = "Cape \n Agulhas", x = 20.0, y = -35.1, size = 2., colour = "white") +
  #annotate("text", label = "False Bay", x = 18.66, y = -34.55, size = 2.2) +
  annotate("text", label = "Algoa \n Bay", x = 26.4, y = -34.0, size = 2., colour = "white") +
  annotate("text", label = "Cape \n Point", x = 18.0, y = -34.4, size = 2, colour = "white") +
  annotate("text", label = "Hamburg", x = 28.3, y = -33.28611111, size = 2, colour = "white") +
  theme(#panel.background = element_rect(fill = "steelblue4", colour = NA),
        #panel.border = element_blank(),
        #axis.title = element_blank(),
        #axis.text = element_blank(),
        #axis.ticks = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        #plot.title = element_blank(),
        #legend.text = element_text(size = 7, colour = "White"),
        #legend.title = element_text(size = 7, colour = "White"),
        legend.key = element_rect(colour = NA, size = 0.2),
        legend.key.height = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.justification = c(1,0), legend.position = c(0.5, 0.4)) +
  coord_cartesian(xlim = sa_lons, ylim = sa_lats)
  #coord_map(xlim = sa_lons, ylim = sa_lats, projection = "mercator")
#f1
ggsave("graph/boundingBoxesAll.pdf", width = 8, height = 4)
