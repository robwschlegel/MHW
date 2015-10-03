###############################################################################
## DESCRIPTION: This function subsets a data.frame (data) with a premade list of sites (sites) extracting the sites of interest and rearranging them in the correct west to east order.
## USAGE: You must be in the project root directory to use this function. 
## ARGUMENTS: The (data) must have columns labeled "site" and "src"
## DETAILS:
## VALUE:
## AUTHORS(S):
## REFERENCE(S):
## EXAMPLE(S)
##############################################################################

seqSites <- function(data, sites = data$site, src = data$src, ...){
  seqSite <- read.csv("setupParams/site_list_v3.4.csv")
  seqSite$index <- 1:length(seqSite$site)
  seqSite$index2 <- paste(seqSite$site, seqSite$src, sep = "/ ")
  data$index2 <- paste(data$site, data$src, sep = "/ ")
  seqSite$index2 <- reorder(seqSite$index2, seqSite$index) # Convert character vector into ordered factor vector
  seqSite2 <- droplevels(subset(seqSite, (index2 %in% data$index2)))
  data2 <- droplevels(subset(data, (index2 %in% seqSite2$index2)))
  index <- data.frame()
  for(i in 1:length(levels(seqSite2$site))){
    data3 <- droplevels(subset(data2, site == levels(seqSite2$site)[i]))
    seqSite3 <- droplevels(subset(seqSite2, site == levels(seqSite2$site)[i]))
    for(j in 1:length(levels(as.factor(seqSite3$src)))){
      seqSite4 <- subset(seqSite3, src == seqSite3$src[j])
      data4 <- subset(data3, src == seqSite3$src[j])
      data4$index <- seqSite4$index
      index <- rbind(index, data4) 
    }
  }
  index <- index[order(index$index),]
  #index$index2 <- paste(index$site, index$src, sep = "/ ")
  #index$index <- order(index$index) # Reorder factors, used for plotting
  index$index <- NULL; index$index2 <- NULL # Remove both index columns
  rownames(index) <- NULL
  return(index)
}