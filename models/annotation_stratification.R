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
rm(data)
annots <- c('exon','intergenic','intron','non-unique','te')
results <- NULL
#
###############################################################################
#
# For loop across annotations
#
###############################################################################
#
for (j in 1:5){

# set up training and test sets
#
set.seed(1234)
ann.data <- dplyr::filter(clear.data, Annotation == annots[j])
ann.data <- train.partition(ann.data,0.7)
#
# Set out number of classes in system
#
classes <- 5
#
################################################################################
#
# Training Phase
#
################################################################################
#
# Add actual results for rel meth quants, adjust for introns as too big a 
# proportion have no relative methylation
#
if (j == 3){
        ann.data <- quant.bin(ann.data,index=1,split = ncol(ann.data),
                              title = "Rel_meth_quant_act",class.number = classes, zero.exclude = T)
} else{
        ann.data <- quant.bin(ann.data,index=1,split = ncol(ann.data),
                                title = "Rel_meth_quant_act",class.number = classes)
}
#
#
temp.data <- dplyr::filter(ann.data, set =="training")

if (j == 3){
        train.set <- quant.bin(temp.data, index = 1, split = ncol(temp.data), 
                               title = "class_train",class.number = classes, zero.exclude = T)
}else{
        train.set <- quant.bin(temp.data, index = 1, split = ncol(temp.data), 
                               title = "class_train",class.number = classes)      
}
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
# Prediction
###############################################################################
#
# Complete table of probabilites
#
probs <- matrix(nrow = nrow(ann.data),ncol = classes*3+3)
for(i in 1:classes){
        probs[,i] <- dnbinom(ann.data[,2], size = param_CpG[[i]]$estimate[1],
                             mu = param_CpG[[i]]$estimate[2])
}
rm(param_CpG,i)
#
# Train Mappab Score with gaussian kernal
#
for(i in 1:classes){
        probs[,i+classes] <- kern.density(ann.data[,4],
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
predictions <- as.numeric(apply(probs[,(2*classes +1):(3*classes)],1,which.max))

#
# Append raw data with the predicitons
#
ann.data <- cbind(ann.data,predictions)
rm(probs)
tmp <- ncol(ann.data)
names(ann.data)[tmp] <- c("Predicted_meth_quant")
rm(tmp)
#
################################################################################
#
# Score on test set
#
################################################################################
#
test.data <- dplyr::filter(ann.data, set =="test")
results[[j]] <- classifier.test(test.data, index = c(8,9))
rm(ann.data,test.data,predictions)

}

results
