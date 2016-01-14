### Distance functions found at:
# http://r.789695.n4.nabble.com/Geographic-distance-between-lat-long-points-in-R-td3442338.html --
## According to that blogpost, I am using the function gcd.hf using the haversine formula. I wrapped it up in a function called CalcDists so that I can get a distance matrix between N sites.

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian lat/lon using the haversine formula:
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

# Function to calculate matrix of distances between each two sites usng the haversine function, above:
CalcDists <- function(latlongs) {
  name <- list(rownames(latlongs), rownames(latlongs))
  n <- nrow(latlongs)
  z <- matrix(0, n, n, dimnames = name)
  for (i in 1:n) {
    for (j in 1:n) z[i, j] <- gcd.hf(long1 = latlongs[i, 1],
                                     lat1 = latlongs[i, 2], long2 = latlongs[j, 1], lat2 = latlongs[j,2])
  }
  z <- as.dist(z)
  return(z)
}

# Distances between consecutive pairs of sites in a list. This function requires a data.frame with site in column 1, lon in column 2 and lat in column 3. Beforehand lats and lons in degrees need to be convereted to radians using the haversine function above.
PairsDists <- function(latlongs) {
  n <- nrow(latlongs)
  z <- matrix(0, n, 1, dimnames = list(latlongs[,1]))
  for (i in 1:n) {
    z[i] <- gcd.hf(long1 = latlongs[i, 2], lat1 = latlongs[i, 3],
                   long2 = latlongs[i+1, 2], lat2 = latlongs[i+1,3])
  }
  return(z)
}

# Distance between multisites
## NB: Automatically uses meta data in "proc/multisite_meta_v3.3.csv", which was compiled manually
## You must alter this file if you want different results
MultiDists <- function(){
  multisite_meta <- read.csv("proc/multisite_meta_v3.3.csv")
  multisite_meta[,c(5,6)] <- deg2rad(multisite_meta[,c(5,6)]) # Converts lon/lat from degree to radian
  multisite_meta <- seqSites(multisite_meta, multisite_meta$site) # Reorder factors so that sorting function works correctly
  # Calculate distance between multisite points
  multisite_dists <- data.frame()
  for (i in 1:length(levels(multisite_meta$site))){
    data <- droplevels(subset(multisite_meta, site == levels(multisite_meta$site)[i]))
    if (length(levels(as.factor(data$src))) == 2){
      dist <- gcd.hf(long1 = data[1, 5], lat1 = data[1, 6],
                     long2 = data[2, 5], lat2 = data[2, 6])
      depth <- abs(data$depth[1] - data$depth[2]) # Show difference in depth
      if (data$type[1] == "UTR" & data$type[2] == "UTR") {tool <- "UTR"
      } else if (data$type[1] == "thermo" & data$type[2] == "thermo") {tool <- "thermo"
      } else {tool <- "diff"} # Show if measurement devices are same or different
      src <- paste(data$src[1], data$src[2], sep = "/") # Show sources
      site <- levels(data$site)[1] # Show site name
      data2 <- data.frame(site, src, dist, depth, tool) # Bind it together
      multisite_dists <- rbind(multisite_dists, data2) # Add it to the data.frame
    } else if (length(levels(as.factor(data$src))) == 3) {
      data <- data[order(data$src),]
      x1 <- droplevels(subset(data, src == data$src[1])) # Remove all other sources of data
      x2 <- droplevels(subset(data, src == data$src[1])) # Remove all other sources of data
      y1 <- droplevels(subset(data, src == data$src[2])) # Remove all other sources of data
      y2 <- droplevels(subset(data, src == data$src[2])) # Remove all other sources of data
      z1 <- droplevels(subset(data, src == data$src[3])) # Remove all other sources of data
      z2 <- droplevels(subset(data, src == data$src[3])) # Remove all other sources of data

      ### First round
      dist <- gcd.hf(long1 = x1[1, 5], lat1 = x1[1, 6],
                     long2 = y1[1, 5], lat2 = y1[1, 6])
      depth <- abs(x1$depth[1] - y1$depth[1]) # Show difference in depth
      if (x1$type[1] == "UTR" & y1$type[1] == "UTR") {tool <- "UTR"
      } else if (x1$type[1] == "thermo" & y1$type[1] == "thermo") {tool <- "thermo"
      } else {tool <- "diff"} # Show if measurement devices are same or different
      src <- paste(x1$src[1], y1$src[1], sep = "/") # Show sources
      site <- levels(x1$site)[1] # Show site name
      data2 <- data.frame(site, src, dist, depth, tool) # Bind it together
      multisite_dists <- rbind(multisite_dists, data2) # Add it to the data.frame

      ### Second round
      dist <- gcd.hf(long1 = x2[1, 5], lat1 = x2[1, 6],
                     long2 = z1[1, 5], lat2 = z1[1, 6])
      depth <- abs(x2$depth[1] - z1$depth[1]) # Show difference in depth
      if (x2$type[1] == "UTR" & z1$type[1] == "UTR") {tool <- "UTR"
      } else if (x2$type[1] == "thermo" & z1$type[1] == "thermo") {tool <- "thermo"
      } else {tool <- "diff"} # Show if measurement devices are same or different
      src <- paste(x2$src[1], z1$src[1], sep = "/") # Show sources
      site <- levels(x2$site)[1] # Show site name
      data2 <- data.frame(site, src, dist, depth, tool) # Bind it together
      multisite_dists <- rbind(multisite_dists, data2) # Add it to the data.frame

      ### Third round
      dist <- gcd.hf(long1 = y2[1, 5], lat1 = y2[1, 6],
                     long2 = z2[1, 5], lat2 = z2[1, 6])
      depth <- abs(y2$depth[1] - z2$depth[1]) # Show difference in depth
      if (y2$type[1] == "UTR" & z2$type[1] == "UTR") {tool <- "UTR"
      } else if (y2$type[1] == "thermo" & z2$type[1] == "thermo") {tool <- "thermo"
      } else {tool <- "diff"} # Show if measurement devices are same or different
      src <- paste(y2$src[1], z2$src[1], sep = "/") # Show sources
      site <- levels(y2$site)[1] # Show site name
      data2 <- data.frame(site, src, dist, depth, tool) # Bind it together
      multisite_dists <- rbind(multisite_dists, data2) # Add it to the data.frame
    }
  }
  multisite_dists$dist <- round(multisite_dists$dist, 2)
  colnames(multisite_dists) <- c("site", "src", "dist (km)", "depth diff (m)", "tool")
  return(multisite_dists)
}