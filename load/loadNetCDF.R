#############################################################################
## This script does:
# 1. loads netcdf file from eric and coastal metadata;
# 2. finds nearest neighbour in the ocean;
# 3. saves coords as .csv for eric
#############################################################################

#############################################################################
## DEPENDS ON:
require(ggplot2); require(ncdf); require(FNN); require(reshape2)
source("setupParams/theme.R")
# "setupParams/mhw_mcs_census_SA.nc"
# "data/metaData2.Rdata"
# "graph/sa_shore.Rdata"
#############################################################################

#############################################################################
## USED BY:
# 
#############################################################################

#############################################################################
## CREATES:
#"data/site_pixels2.csv"
#############################################################################

#############################################################################
## 1. loads netcdf file from eric and coastal metadata
ncsst <- open.ncdf("setupParams/mhw_mcs_census_SA.nc")
lat <- ncsst$dim$lat$vals
lon <- ncsst$dim$lon$vals
grid <- expand.grid(lon, lat)
z <- melt(get.var.ncdf(ncsst, "MHW_freq"))
grid <- cbind(grid,mhw = z[,3])

# Load metadata for comparison purposes
load("data/metaData2.Rdata")
metaData3 <- metaData2[20:21,]

#############################################################################
## 2. finds nearest neighbour in the ocean

coords1 <- knnx.index(as.matrix(grid[1:2]), as.matrix(metaData3[,3:4]), k = 4)
coords2 <- grid[coords1,]

# Setup up environment for plotting
latSA <- c(-35.5, -26); lonSA <- c(14, 34)

# Load SA map data
load("graph/south_africa_coast.RData")

p <- ggplot() + coord_equal() + theme_bw() +
  geom_polygon(data = south_africa_coast, aes(x = long, y = lat, group = group), 
               colour = "black", fill = "grey80") +
  geom_point(data = coords3, aes(x = Var1, y = Var2), colour = "red", size = 3.2) +
  geom_point(data = metaData3, aes(x = lon, y = lat), colour = "blue", size = 3.2) +
  coord_map(xlim = lonSA, ylim = latSA, projection = "mercator")
p

# After visually inspecting the match-ups, manually select and touch up
coords3 <- coords2[c(4,7),1:2]
names(coords3) <- c("lon", "lat")
coords3$site <- rev(metaData3$site)
coords3 <- coords3[2:1,c(3,1:2)]

#############################################################################
## 3. saves coords as .csv for eric
write.csv(coords3, "data/site_pixels2.csv", row.names = F)
