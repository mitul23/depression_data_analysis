library(dplyr)
library(tidyr)
library(data.table)
library(dmetar)

dat <- fread("discovery_metav3.0.meta.wAF.wBETA.tsv") %>% as_tibble() 
head(dat)

SE <- se.from.p(effect.size =  dat$OR, 
                p = dat$P, 
                N = dat$N,  
                effect.size.type= "ratio")

data <- cbind(dat, SE)
head(data)
sum(is.na(data$logStandardError))
idx <- which(is.na(data$logStandardError))
dim(data[idx,])
summary(data$logStandardError)

dat_nans <- data %>% filter(OR == 0)
dim(dat_nans)
write.csv(dat_nans,"dat_withnas.csv")
