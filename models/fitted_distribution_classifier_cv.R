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
clear.data <- select.narm(data,index = 1:5,select = FALSE) %>%
        dplyr::filter(relative_meth > 0)
rm(data)
#
# set up training, CV and test sets
#
set.seed(1234)
clear.data <- cv.partition(clear.data,0.6,0.2)
#
# Set out number of classes in system
#
CV.score <- NULL
meth.breaks <- NULL
meth.count <- NULL
#
################################################################################
#
# Cross Validation For Loop
#
################################################################################
#
for (classes in 3:10){
#
# Add actual results for rel meth quants
#
clear.data <- quant.bin(clear.data,index=1,split = ncol(clear.data),
                        title = "Rel_meth_quant_act",class.number = classes)
meth.breaks[[classes-2]] <- stats::quantile(clear.data[,1],
                                            probs = seq(0,1,
                                                        length.out = (classes +1)))
meth.count[[classes-2]] <- table(clear.data[,ncol(clear.data)])

#
#
#
# Training Phase
################################################################################
#
temp.data <- dplyr::filter(clear.data, set =="training")
train.set <- quant.bin(temp.data, index = 1, split = ncol(temp.data), 
                       title = "class_train",class.number = classes)
rm(temp.data)
#
# Train CpG COunt Model with negative Binomial
#
param_CpG <- NULL
for (i in 1:classes){
        param_CpG[[i]] <- dplyr::filter(train.set, class_train == i)[,2] %>%
                fitdistrplus::fitdist("nbinom",method = "mle")
}
#
#
#
# Scoring Phase
###############################################################################
#
# Complete table of probabilites
#
probs <- matrix(nrow = nrow(clear.data),ncol = classes*3+3)
for(i in 1:classes){
        probs[,i] <- dnbinom(clear.data[,2], size = param_CpG[[i]]$estimate[1],
                                mu = param_CpG[[i]]$estimate[2])
}
rm(param_CpG,i)
#
# Train Mappab Score with emprical ideas
#
for(i in 1:classes){
        probs[,i+classes] <- kern.density(clear.data[,4],
                                    dplyr::filter(train.set,class_train == i)[,4])
}
rm(train.set,i)
#
# Calculate scores for the combined distributions 
#
for(i in 1:classes){
        probs[,i+2*classes] <- probs[,i]*probs[,i+classes]
}
rm(i)
#
# Normalise Probabilities
#
probs[,3*classes + 1] <- apply(probs[,(2*classes +1):(3*classes)],1,sum)
probs[,(2*classes +1):(3*classes)] <- probs[,(2*classes +1):(3*classes)]/probs[,3*classes+1]
#
# Pick most likely Class and what the socre for this is
#
probs[,3*classes + 2] <- apply(probs[,(2*classes +1):(3*classes)],1,which.max)
probs[,3*classes + 3] <- apply(probs[,(2*classes +1):(3*classes)],1,max)
#
# Append raw data with the predicitons
#
clear.data <- cbind(clear.data,probs[,(3*classes+2):(3*classes+3)])
rm(probs)
tmp <- ncol(clear.data)
names(clear.data)[(tmp-1):tmp] <- c("Predicted_meth_quant","probability")
rm(tmp)
#
#
# Cross-Validation step
################################################################################
#
temp.data <- dplyr::filter(clear.data,set == "CV")
CV.score[[classes-2]] <- classifier.test(temp.data,index = c(8,9),top.num = classes)



clear.data <- clear.data[,1:7]
}
################################################################################
#
# Cross Validation For Loop Ends above
#
################################################################################

CV.score