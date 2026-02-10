load("~/subset_list.RData")
spe <- subset_list$sham_b1
keep <- sample(dim(spe)[2],size=dim(spe)[2]/4)
spe_sub <- spe[,keep]

#repeat 10times


