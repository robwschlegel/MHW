# Create figures of the desired bounding boxes that will be used to extract SST that will be compared to in situ data

# Load required packages
require(ggplot2); require(geosphere)

# Load necessary files
load("graph/sa_provinces_new.RData")
sa_provinces_new$index <- 1:12 # Reduce it by 92%
sa_provinces_new <- droplevels(subset(sa_provinces_new, index == 1))
load("graph/sa_shore.Rdata") # HiRes for figure
load("graph/south_africa_coast.RData") # LowRes for bounding box
names(south_africa_coast)[1] <- "lon" # Change to differentiate from bounding box edge

# Subset coastline depending on desired bounding box range
  # I decided to do this manually for now...
wc <- south_africa_coast[296:410,]
sc <- south_africa_coast[132:295,]
ec <- south_africa_coast[23:132,]

# Calculate 50km distance line from coastline for appropriate headings per coast
  # Using logic loop as aaply efforts were not succesful...
# x <- aaply(wc, 1, destPoint(p = c(wc[,1], wc[,2]), b = 270, d = 50000))

boundingBox <- function(dat, bearing = 0, distance = 50000){
  #df <- data.frame()
  for(i in rev(1:length(dat$lon))){
    x <- destPoint(p = c(dat$lon[i], dat$lat[i]), b = bearing, d = distance)
    y <- cbind(x, dat[i,3:7])
    dat <- rbind(dat, y)
  }
  return(dat)
}

wcBox <- boundingBox(wc, 270)
scBox <- boundingBox(sc, 180)
ecBox <- boundingBox(ec, 90)

# Setup up environment for plotting
latSA <- c(-35.5, -26); lonSA <- c(14, 34)

# Create figure showing the different bounding boxes
f1 <- ggplot() + coord_equal() + theme_bw() +
  geom_polygon(data = south_africa_coast, aes(x = lon, y = lat, group = group), 
               colour = "black", fill = "grey80") +
  geom_path(data = sa_provinces_new, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = wcBox, aes(x = lon, y = lat, group = group, alpha = 20), 
               colour = "blue", fill = "blue") +
  geom_polygon(data = scBox, aes(x = lon, y = lat, group = group, alpha = 20), 
               colour = "green", fill = "green") +
  geom_polygon(data = ecBox, aes(x = lon, y = lat, group = group, alpha = 20), 
               colour = "orange", fill = "orange") +
  guides(alpha = FALSE) +
  annotate("text", label = "Cape \n Point", x = 18.0, y = -34.7, size = 2) +
  annotate("text", label = "Hamburg", x = 28.3, y = -33.28611111, size = 2) +
  coord_map(xlim = lonSA, ylim = latSA, projection = "mercator")
f1
ggsave("graph/boundingBoxes.pdf", width = 6, height = 4)
