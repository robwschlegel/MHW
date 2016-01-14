###############################################################################
## DESCRIPTION: Given coordinates, these functions finds transect data for further use
## USAGE: To discern the heading of shore normal transects and other applications
## ARGUMENTS:
## DETAILS:
## VALUE:
## AUTHORS(S): Robert Schlegel
## REFERENCE(S):
## EXAMPLE(S):
## DEPENDS: "graph/south_africa_coast.RData", "data/bathy/bathy.RData"
## USED BY:
##############################################################################

require(FNN)
require(fossil)
require(rgeos)
require(maptools)
require(geosphere)
source("func/earthdist.R")

# low-res Africa Coastline
load("graph/africa_coast.RData")

# Download mid-res bathymetry data
# sa_lat <- c(-38, -24.5); sa_lon <- c(11.5, 35.5)
# sa_bathy <- as.xyz(getNOAA.bathy(lon1 = sa_lon[1], lon2 = sa_lon[2], lat1 = sa_lat[1], lat2 = sa_lat[2], resolution = 4))
# colnames(sa_bathy) <- c("lon", "lat", "depth")
# sa_bathy <- sa_bathy[sa_bathy$depth <= 0,]
# save(sa_bathy, file = "data/bathy/sa_bathy.RData")
load("data/bathy/sa_bathy.RData")
#plot(as.bathy(sa_bathy), image = T)

##############################################################################
# This function takes one site (e.g. one set of lon/lats) and calculates a shore normal transect
shore.normal.transect <- function(site, width = 2){
  # Find the site on the coastline and it's nearest neighbour points
  coords <- data.frame(lon = site$lon, lat = site$lat)
  coords2 <- knnx.index(africa_coast[,1:2], as.matrix(coords), k = 1)
  coords3 <- data.frame(site = site$site, africa_coast[c(coords2-width, coords2+width),]) # This line controls how wide the spread of points used to determine the shore-normal transect is. A wider spread gives a more even, less sensitive transect
  coords3 <- coords3[2:1,1:3]
  # Define the shore normal transect bearing
  heading <- earth.bear(coords3[1,2], coords3[1,3], coords3[2,2], coords3[2,3]) + 90
  if(heading >= 360){
    heading <- heading-360
  } else {
    heading <- heading
  }
  heading2 <- data.frame(site = site$site, lon = site$lon, lat = site$lat, heading)
  return(heading2)
}

##############################################################################
# This function takes one site (e.g. one set of lon/lats) and calculates a shore normal transect
  # It then extracts a lat/ lon point every X kilometres until reaching a specified isobath
transect.pixel.isobath <- function(site, distance = 25000, isobath = -200){
  # Extract coordinates
  coords <- data.frame(lon = site$lon, lat = site$lat)
  # Find lon/ lats every X metres until the chosen isobath is reached
  pixels <- data.frame()
  deep <- 0
  while(deep > isobath){
    # Find depth at distance step
    coords2 <- as.data.frame(destPoint(p = coords, b = site$heading, d = distance))
    sitesIdx <- knnx.index(sa_bathy[,1:2], as.matrix(coords2), k = 1)
    bathy2 <- sa_bathy[sitesIdx,]
    bathy2 <- bathy2[complete.cases(bathy2$depth),]
    #if(nrow(bathy2) < 1){
    #  bathy2 <- data.frame(coords2, depth = 0)
    #}else{
    #  bathy2 <- bathy2
    #}
    bathy3 <- data.frame(site = site$site, lon = bathy2$lon, lat = bathy2$lat, 
                         heading = site$heading, depth = bathy2$depth)
    pixels <- rbind(pixels, bathy3)
    coords <- coords2
    deep <- bathy2$depth
  }
  pixels2 <- pixels[pixels$depth > -200,] # This may remove the only row of data
  if(nrow(pixels2) < 1){
    pixels2 <- data.frame(site, depth = NA)
  }else{
    pixels2 <- pixels2
  }
  return(pixels2)
}