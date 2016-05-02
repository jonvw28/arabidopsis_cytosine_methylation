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
        select.narm(c(1,2,3,7))
#
KS.CG <- list()
KS.CGprop <- list()
KS.map <- list()
classes <- list()
#
# Code below will apply K-S tests to the data split into 3 through 10 quantiles
# These will be applied to CpG content, CpG proportion and mappability. The 
# results are contained in the lists outputted. 3 quantiles goes to element 1,
# 4 quantiles to element 2 etc.
#
# Classes then gives a list of which quantiles are being compared for each test 
# result. The structure of classes is exactly the same as that for the
# results.
#
#
for(num.bins in 3:10){
        
        t.data <- quant.bin(data, index = 1, title = "rel_meth_bin", 
                            zero.exclude = T, split = 1,class.number = num.bins)
        
        tmp.CG <- list()
        tmp.CGprop <- list()
        tmp.map <- list()
        tmp.classes <- list()
        #
        # CG count
        #
        k <- 1
        for(i in 1:(num.bins-1)){
                for(j in (i+1):num.bins){
                        tmp1 <- dplyr::filter(t.data, rel_meth_bin == i)
                        tmp2 <- dplyr::filter(t.data, rel_meth_bin == j)
                        tmp.CG[[k]] <- stats::ks.test(tmp1[,4],tmp2[,4])
                        tmp.classes[[k]] <- c(i,j)
                        k <- k + 1
                }
        }
        #
        rm(i,j,k,tmp1,tmp2)
        classes[[num.bins - 2]] <- tmp.classes
        KS.CG[[num.bins - 2]] <- tmp.CG
        rm(tmp.CG,tmp.classes)
        #
        # CG prop
        #
        k <- 1
        for(i in 1:(num.bins-1)){
                for(j in (i+1):num.bins){
                        tmp1 <- dplyr::filter(t.data, rel_meth_bin == i)
                        tmp2 <- dplyr::filter(t.data, rel_meth_bin == j)
                        tmp.CGprop[[k]] <- stats::ks.test(tmp1[,5],tmp2[,5])
                        k <- k + 1
                }
        }
        #
        rm(i,j,k,tmp1,tmp2)
        KS.CGprop[[num.bins - 2]] <- tmp.CGprop
        rm(tmp.CGprop)
        #
        # mappab
        #
        k <- 1
        for(i in 1:(num.bins-1)){
                for(j in (i+1):num.bins){
                        tmp1 <- dplyr::filter(t.data, rel_meth_bin == i)
                        tmp2 <- dplyr::filter(t.data, rel_meth_bin == j)
                        tmp.map[[k]] <- stats::ks.test(tmp1[,3],tmp2[,3])
                        k <- k + 1
                }
        }
        #
        rm(i,j,k,tmp1,tmp2)
        KS.map[[num.bins - 2]] <- tmp.map
        rm(tmp.map)
        #
}
rm(num.bins,t.data)