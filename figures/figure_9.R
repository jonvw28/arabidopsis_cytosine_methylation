At_tiles.data <- readRDS("data_complete.RData")
#
#LIBRARIES
#
source('functions_read_in_first.R')
#

if(!require("dplyr")){
        install.packages("dplyr")
}
library(dplyr)
#
# Apply filtering to data
#
data <- dplyr::select(At_tiles.data,c(relative_meth,mappab_20mer_1msh,cytosinesCountCG,
                                      cytosinesCountCHG,cytosinesCountCHH)) %>%
        dplyr::mutate(totalC = cytosinesCountCG + cytosinesCountCHG + cytosinesCountCHH)
data <- col.ratio(data,c(3,6),6,"CpG_prop") %>%
        select.narm(c(1,2,3,7)) %>%
        quant.bin(index = 1,title = "rel_meth_quant", zero.exclude = F, split = 4,
                  class.number = 5)
rm(At_tiles.data)
#
# plot histograms and fitted Gaussian kernal distributions
#
cols <-  adjustcolor(rainbow(5), alpha.f=0.3)
par(mfrow = n2mfrow(5))
labels <- character(5)
#
for(i in 1:5){
        labels[i] <- paste(20*(i-1),"-",20*i,sep="")
}

for (i in 1:5){
        temp <- dplyr::filter(data,rel_meth_quant == i)
        hist(temp$mappab_20mer_1msh, ylim = c(0,12), xlim = c(0,1),freq = FALSE,
             breaks=50, col = cols[i],
             main = paste("Histogram of mappability ",labels[i]," percentile",sep=""),
             xlab = "mappability")
        lines(stats::density(temp$mappab_20mer_1msh))
        rm(temp)
        
}