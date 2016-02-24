# Co-occurrence figure
mhwCoastCO <- read.csv("data/mhwCoastCO.csv")
mhwCoastCO$event <- "MHW"
mcsCoastCO <- read.csv("data/mcsCoastCO.csv")
mcsCoastCO$event <- "MCS"
allCoastCO <- rbind(mhwCoastCO, mcsCoastCO)
allCoastCO$index <- paste(allCoastCO$lag, allCoastCO$direction, sep = "_")
allCoastCO$coast <- factor(allCoastCO$coast, levels = c("west", "south", "east"))
allCoastCO$direction <- factor(allCoastCO$direction, levels = c("b", "x", "a"))

figure5 <-  ggplot(data = allCoastCO, aes(x = quantile, y = proportion)) + bw_update +
  geom_line(aes(colour = as.factor(index))) + 
  geom_point(aes(colour = as.factor(index))) +
  facet_grid(event + coast ~ direction) +
  ylab("proportion") + xlab("quantile (%)")
figure5
ggsave("graph/figure5.pdf")
