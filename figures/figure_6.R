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
# Select relevant data and then apply quantile binning (dealing with zeros in
# bottom class)
#
data <- dplyr::select(At_tiles.data,cytosinesCountCG,mappab_20mer_1msh,
                      relative_meth) %>%
        quant.bin(index = 3, title = "rel_meth_quant", zero.exclude = T, split = 3,
                  class.number = 10) %>%
        select.narm(index = 1:3,select = FALSE)
rm(At_tiles.data,breakpoints)
#
# plot a grid of smoothed scatterplots
#
par(mfrow = c(5,2))
size <- numeric(10)
title <- character(10)
#
title[1] <- "Bottom 10th percentile\n and complete loss"
for(i in 2:10){
        title[i] <- paste((i-1)*10,"th - ",10*i,"th percentile",sep="")
}
rm(i)

for(i in 1:10){
        temp <- dplyr::filter(data, rel_meth_quant == i)
        smoothScatter(temp$mappab_20mer_1msh,temp$cytosinesCountCG,
                      xlab = 'mappability', ylab = 'CpG Count',
                      main = title[i])
        size[i] <- nrow(temp)
        rm(temp)
}
rm(i,title)