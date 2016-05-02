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
# Apply filters
#
GEL.data <- dplyr::filter(At_tiles.data, cg_class == 1 & chg_class == 0 & chh_class == 0) %>%
        dplyr::select(mappab_20mer_1msh,cytosinesCountCG)

TEL.data <- dplyr::filter(At_tiles.data, cg_class == 1 & chg_class == 1 & chh_class == 1) %>%
        dplyr::select(mappab_20mer_1msh,cytosinesCountCG)
#
# plot figure 5
#
par(mfcol = c(2,1))
cols <-  rainbow(7)
#
smoothScatter(x = GEL.data[,1], y = GEL.data[,2], nbin =512,
              xlab="Mappability Score", ylab="CpG Count",
              colramp = colorRampPalette(c("white",cols[2])),
              ylim = c(0,75), xlim = c(0,1),
              main="Smoothed scatterplot of CpG context cytosine count \nand mappability score for GEL tiles"
)
#
smoothScatter(x = TEL.data[,1], y = TEL.data[,2],nbin = 512,
              xlab="Mappability Score", ylab="CpG Count",
              colramp = colorRampPalette(c("white",cols[5])),
              ylim = c(0,75), xlim = c(0,1),
              main="Smoothed scatterplot of CpG context cytosine count \nand mappability score for TEL tiles"
)