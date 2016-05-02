
At_tiles.data <- readRDS("data_complete.RData")

model6 <- lm(relative_meth ~ (cytosinesCountCG + mappab_20mer_0msh)^2,
             At_tiles.data, na.action = na.omit)

model7 <- lm(relative_meth ~ (cytosinesCountCG + mappab_20mer_1msh)^2,
             At_tiles.data, na.action = na.omit)

model8 <- lm(relative_meth ~ (cytosinesCountCG + mappab_20mer_2msh)^2,
             At_tiles.data, na.action = na.omit)

model9 <- lm(relative_meth ~ (cytosinesCountCG + mappab_20mer_3msh)^2,
             At_tiles.data, na.action = na.omit)

model10 <- lm(relative_meth ~ (cytosinesCountCG + mappab_200mer)^2,
             At_tiles.data, na.action = na.omit)

