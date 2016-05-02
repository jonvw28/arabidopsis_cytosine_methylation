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
GEL <- dplyr::filter(At_tiles.data,cg_class == T, chg_class == F, 
                     chh_class == F) %>%
        dplyr::filter(relative_meth <= 2)
TEL <- dplyr::filter(At_tiles.data,cg_class == T, chg_class == T,
                      chh_class == T) %>%
        dplyr::filter(relative_meth <= 2)
rm(At_tiles.data)
#
# plot histograms
#
par(mfcol = c(2,1))
cols <-  adjustcolor(rainbow(7), alpha.f=0.3)
#
hist(GEL[,16],freq = F,ylim = c(0,10),
     main = expression("Histogram of relative CpG context methylation in "~italic("met")~"1-1 mutant vs WT for GEL tiles"),
     ylab = "Relative Frequency", 
     xlab = expression("Relative methylation of cytosines in CpG context in"~italic("met")~"1-1 mutant as proportion of WT levels"),
     col = cols[2])
hist(TEL[,16],freq = F,ylim = c(0,2),
     main = expression("Histogram of relative CpG context methylation in "~italic("met")~"1-1 mutant vs WT for TEL tiles"),
     ylab = "Relative Frequency", 
     xlab = expression("Relative methylation of cytosines in CpG context in"~italic("met")~"1-1 mutant as proportion of WT levels"),
     col = cols[5])