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
# Apply filtering to the data
#
data <- dplyr::select(At_tiles.data,c(relative_meth,mappab_20mer_1msh,cytosinesCountCG,
                               cytosinesCountCHG,cytosinesCountCHH)) %>%
        dplyr::mutate(totalC = cytosinesCountCG + cytosinesCountCHG + cytosinesCountCHH)
data <- col.ratio(data,c(3,6),6,"CpG_prop") %>%
        select.narm(c(1,2,3,7))
#
# Regime below will plot ecdfs for CpG count, mappability and CpG proportion
# for 3 through 10 quantiles and save all of these as png files
#
for(num.bins in 3:10){

t.data <- quant.bin(data, index = 1, title = "rel_meth_bin", zero.exclude = T, 
                  split = 1,
                  class.number = num.bins)

CG.emp <- list()
map.emp <- list()
CGprop.emp <- list()

for(i in 1:num.bins){
        temp <- dplyr::filter(t.data, rel_meth_bin == i)
        CG.emp[[i]] <- stats::ecdf(temp[,4])
        map.emp[[i]] <- stats::ecdf(temp[,3])
        CGprop.emp[[i]] <- stats::ecdf(temp[,5])
}
rm(temp,i)

library(RColorBrewer)

tcols <- RColorBrewer::brewer.pal(9,"OrRd")
mycols <- colorRampPalette(c(tcols[3],tcols[9]))(num.bins)
rm(tcols)

temp.num <- vector("numeric", length = num.bins)
labels <- vector("character", length = num.bins)
for (i in 1:num.bins){
        temp.num[i] <- round(i*100/num.bins,1)
}
labels[1] <- paste("Bottom ",temp.num[1], "% and complete loss", sep = "")
for(i in 2:num.bins){
        labels[i] <- paste(temp.num[i-1],"%-",temp.num[i],"%",sep="")
}
rm(temp.num,i)
#
# CpG Count
#
png(paste("CpG_count_ecdf_",num.bins,"_bins.png",sep = ""),width = 960, 
    height = 960)

plot(CG.emp[[1]], verticals = TRUE, do.points = FALSE, col = mycols[1],
     xlab = "Cytosines in CpG Context - Count", ylab = "Cumulative Probability",
     main = paste("Empirical Cumulative Distribution Functions for CpG Count\nSeperated by ",
                  num.bins ," Quantiles of Relative Methylation (Met1-1 vs WT)",
                  sep = ""))
for(i in 2:num.bins){
        plot(CG.emp[[i]], verticals = TRUE, do.points = FALSE,
             col = mycols[i], add = TRUE) 
}

legend("bottomright", legend = labels, lty = rep(1,num.bins), col = mycols)
dev.off()
#
# Mappab
#
png(paste("Mappab_20mer_1msh_ecdf_",num.bins,"_bins.png",sep = ""),width = 960, 
    height = 960)

plot(map.emp[[1]], verticals = TRUE, do.points = FALSE, col = mycols[1],
     xlab = "Mappability Score", ylab = "Cumulative Probability",
     main = paste("Empirical Cumulative Distribution Functions for Mappability\nSeperated by ",
                  num.bins ," Quantiles of Relative Methylation (Met1-1 vs WT)",
                  sep = ""))
for(i in 2:num.bins){
        plot(map.emp[[i]], verticals = TRUE, do.points = FALSE,
             col = mycols[i], add = TRUE) 
}

legend("bottomright", legend = labels, lty = rep(1,num.bins), col = mycols)
dev.off()
#
# CpG proportion
#
png(paste("CpG_prop_ecdf_",num.bins,"_bins.png",sep = ""),width = 960, 
    height = 960)

plot(CGprop.emp[[1]], verticals = TRUE, do.points = FALSE, col = mycols[1],
     xlab = "Cytosines in CpG Context - Proportion of all Cytosines",
     main = paste("Empirical Cumulative Distribution Functions for CpG Proportion\nSeperated by ",
                  num.bins ," Quantiles of Reative Methylation (Met1-1 vs WT)",
                  sep = ""))
for(i in 2:num.bins){
        plot(CGprop.emp[[i]], verticals = TRUE, do.points = FALSE,
             col = mycols[i], add = TRUE) 
}

legend("bottomright", legend = labels, lty = rep(1,num.bins), col = mycols)
dev.off()
#
rm(mycols,CG.emp,CGprop.emp,map.emp,labels,i,t.data)
}