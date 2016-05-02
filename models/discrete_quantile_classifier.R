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
# Select Relevant Data
#
tags <- which.index(At_tiles.data,c('relative_meth','cytosinesCountCG', 
                                    'cytosinesCountCHG', 'cytosinesCountCHH',
                                    'mappab_20mer_1msh','Annotation'))
#
# measure CpG proportion
#
data <- dplyr::select(At_tiles.data,tags) %>%
        dplyr::mutate(prop_CG = cytosinesCountCG/(cytosinesCountCG + 
                                                          cytosinesCountCHG +
                                                          cytosinesCountCHH))
rm(tags,At_tiles.data,breakpoints)
#
# Remove missing data points
#
data <- dplyr::select(data,c(1,2,7,5,6))
clear.data <- select.narm(data,index = 1:5,select = FALSE)
#
# set up training, CV and test sets
#
set.seed(1234)
clear.data <- train.partition(clear.data,0.7)
#
# Set out number of classes in system
#
classes <- c('meth.class' = 5,'CpG.class' = 5,'map.class' = 5)
#
# Add actual results for rel meth quants
#
clear.data <- quant.bin(clear.data,index=1,split = ncol(clear.data),
                        title = "Rel_meth_quant_act",class.number = classes[1])
#
################################################################################
#
# Training Phase
#
################################################################################
#
temp.data <- dplyr::filter(clear.data, set =="training")
#
# Get table of counts for all possible bins
#
count.table <- simple.quant.count(temp.data,index = c(1,2,4), classes = classes,
                                  title = c('rel_meth_quant','CpG_count_quant',
                                            'Mappab_quant'))
#
# Total these for each relative meth quantile and normalise counts to give
# probabilities
#
totals <- dplyr::group_by(count.table,CpG_count_quant,Mappab_quant)%>%
        dplyr::summarise(total = sum(count))

for (i in 1:classes[1]){
        prod <- classes[2]*classes[3]
        count.table[((i-1)*prod + 1):(i*prod),4] <- count.table[((i-1)*prod + 1):(i*prod),4]/totals[,3]
}
rm(totals,prod,i)
#
# Select which class is top for each quantile bin
#
top.prob <- count.table %>%
        dplyr::group_by(CpG_count_quant,Mappab_quant)%>%
        dplyr::summarise(which.max(count))
#
# Calculate cut offs for training set to define classes for all data
#
meth.cuts <- stats::quantile(temp.data[,1],probs=seq(0,1,length.out=classes[1]+1))
CpG.cuts <- stats::quantile(temp.data[,2],probs=seq(0,1,length.out=classes[2]+1))
map.cuts <- stats::quantile(temp.data[,4],probs=seq(0,1,length.out=classes[3]+1))
#
rm(temp.data)
#
# Place all data into quantiles for CpG based on training data
#
temp <- matrix(nrow = nrow(clear.data),ncol = classes[2])
#
for (i in 1:classes[2]){
        temp[,i] <- clear.data[,2] <= rep(CpG.cuts[(i+1)],nrow(clear.data))       
}
clear.data <- clear.data %>% 
        dplyr::mutate(
                train_CpG_quant = rep(classes[2]+1,
                                      nrow(clear.data)) - rowSums(temp))
rm(temp,i)
#
# Place all data into quantiles for Mappab based on training data
#
temp <- matrix(nrow = nrow(clear.data),ncol = classes[3])
#
for (i in 1:classes[3]){
        temp[,i] <- clear.data[,4] <= rep(map.cuts[(i+1)],nrow(clear.data))       
}
clear.data <- clear.data %>% 
        dplyr::mutate(
                train_map_quant = rep(classes[3]+1, 
                                      nrow(clear.data)) - rowSums(temp))
rm(temp,i)
#
################################################################################
#
# Scores Step
#
################################################################################
#
# Pull up most likely class for each tile
#
temp <- matrix(nrow = nrow(clear.data), ncol = 1)
clear.data <- cbind(clear.data,temp)
rm(temp)
tmp.l <- ncol(clear.data)
#
for(i in 1:nrow(clear.data)){
        temp <- top.prob[which(top.prob[,1]==clear.data[i,9]),]
        clear.data[i,tmp.l] <- temp[which(temp[,2]==clear.data[i,10]),3]
        rm(temp)
}
names(clear.data)[tmp.l] <- 'predicted_rel_meth_class'
rm(tmp.l,i)
#
# Score the classifier
#
#
temp.data <- dplyr::filter(clear.data, set =="training")
results <- classifier.test(temp.data, index = c(8,11))
rm(temp.data)

results