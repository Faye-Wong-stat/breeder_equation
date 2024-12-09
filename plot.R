setwd("~/breeder_equation_project/")
library(ggplot2)
library(cowplot)

data1 <- readRDS("estimate_accuracy/accuracy_table_long.rds")
data2 <- readRDS("estimate_accuracy2/accuracy_table_long.rds")
data5 <- readRDS("estimate_accuracy5/accuracy_table_long.rds")
data4 <- readRDS("estimate_accuracy4/accuracy_table_long.rds")

data1 <- cbind(data = "poplar_Rincent", 
               data1)
data2 <- cbind(data = "wheat_Rincent", 
               data2)
data5 <- cbind(data = "Winn", 
               trait = "GrainYield", 
               data5)
colnames(data4)[1] <- "jth_fold"
data4 <- cbind(data = "Krause", 
               trait = "GrainYield", 
               data4)

accuracy <- rbind(rbind(rbind(data1, data2), data5), data4)

accuracy_summary <- aggregate(. ~ data + trait + adjustment + `prediction method`, 
                              data=accuracy[, !colnames(accuracy) %in% c("type", "jth_fold")], 
                              FUN=function(x){c(mean=mean(x), se=sd(x)/sqrt(length(x)))})


accuracy1 



