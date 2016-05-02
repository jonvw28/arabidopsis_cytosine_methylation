At_tiles.data <- readRDS("data_complete.RData")

model1 <- lm(relative_meth ~ cytosinesCountCG,
             At_tiles.data, na.action = na.omit)

model2 <- lm(relative_meth ~ blast,
             At_tiles.data, na.action = na.omit)

model3 <- lm(relative_meth ~ (cytosinesCountCG + blast)^2,
             At_tiles.data,na.action=na.omit)
