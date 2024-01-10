setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)



file.names <- list.files("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", ".txt")
# file.names.poplar <- gsub(".txt", "", file.names[grepl("Poplar", file.names, fixed=T)], )
file.names.wheat <- gsub(".txt", "", file.names[grepl("Wheat", file.names, fixed=T)], )

for (i in 1:length(file.names.wheat)){
  eval(call("<-", as.name(file.names.wheat[i]),
            read.table(paste("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/",
                             file.names.wheat[i],
                             ".txt", sep=""),
                       header=T)))
}
ls()
# [1] "estimate_gcor"                 "file.names"                   
# [3] "file.names.wheat"              "Genotyping_Wheat"             
# [5] "i"                             "NIRS_NormDer_Grain_DRY_Wheat" 
# [7] "NIRS_NormDer_Grain_IRR_Wheat"  "NIRS_NormDer_Leaf_DRY_Wheat"  
# [9] "NIRS_NormDer_Leaf_IRR_Wheat"   "Phenotyping_GrainYield_Wheat" 
# [11] "Phenotyping_HeadingDate_Wheat"



rownames(Phenotyping_GrainYield_Wheat) <- Phenotyping_GrainYield_Wheat$Variety
Phenotyping_GrainYield_Wheat <- Phenotyping_GrainYield_Wheat[, -1]
rownames(Phenotyping_HeadingDate_Wheat) <- Phenotyping_HeadingDate_Wheat$Variety
Phenotyping_HeadingDate_Wheat <- Phenotyping_HeadingDate_Wheat[, -1]
Phenotyping_Wheat <- cbind(Phenotyping_GrainYield_Wheat[, 1:2], Phenotyping_HeadingDate_Wheat[, 1:2])
colnames(Phenotyping_Wheat) <- c("IRR_GrainYield", "DRY_GrainYield", "IRR_HeadingDate", "DRY_HeadingDate")

Genotyping_Wheat <- scale(Genotyping_Wheat)
NIRS_NormDer_Grain_DRY_Wheat <- scale(NIRS_NormDer_Grain_DRY_Wheat)
NIRS_NormDer_Grain_IRR_Wheat <- scale(NIRS_NormDer_Grain_IRR_Wheat)
NIRS_NormDer_Leaf_DRY_Wheat <- scale(NIRS_NormDer_Leaf_DRY_Wheat)
NIRS_NormDer_Leaf_IRR_Wheat <- scale(NIRS_NormDer_Leaf_IRR_Wheat)
str(Phenotyping_Wheat)

kinship <- Genotyping_Wheat %*% t(Genotyping_Wheat) / nrow(Genotyping_Wheat)

set.seed(1) 
five_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
five_folds <- split(five_folds, 1:5)
four_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
four_folds <- split(four_folds, 1:4)
three_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
three_folds <- split(three_folds, 1:3)



accuracy_table_five <- data.frame(fold="five folds", 
                                  trait=rep(colnames(Phenotyping_Wheat), each=5), 
                                  jth_fold=rep(1:5, length(colnames(Phenotyping_Wheat))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear_grain=NA, 
                                  p_acc_pear_leaf=NA, 
                                  p_acc_gcor_grain=NA, 
                                  p_acc_gcor_leaf=NA)
for (i in 1:nrow(accuracy_table_five)){
  tryCatch({
    trait = accuracy_table_five$trait[i]
    jth_fold = accuracy_table_five$jth_fold[i]
    selected = five_folds[[jth_fold]][!is.na(Phenotyping_Wheat[five_folds[[jth_fold]], 
                                                            trait])]
    
    Y <- data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
                    obs=Phenotyping_Wheat[selected, trait])
    kin = kinship[Y$ID, Y$ID]
    rownames(kin) = Y$ID
    colnames(kin) = Y$ID
    
    gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                         Genotyping_Wheat[-selected, ])
    Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
    
    accuracy_table_five$g_acc_pear[i] = cor(Y$obs, Y$pred)
    accuracy_table_five$g_acc_gcor[i] = 
      estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
    if (grepl("IRR", trait)){
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_five$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_five$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_five$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_five$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    } else {
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_five$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_five$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_five$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_five$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    }
  }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
}

accuracy_table_four <- data.frame(fold="four folds", 
                                  trait=rep(colnames(Phenotyping_Wheat), each=4), 
                                  jth_fold=rep(1:4, length(colnames(Phenotyping_Wheat))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear_grain=NA, 
                                  p_acc_pear_leaf=NA, 
                                  p_acc_gcor_grain=NA, 
                                  p_acc_gcor_leaf=NA)
for (i in 1:nrow(accuracy_table_four)){
  tryCatch({
    trait = accuracy_table_four$trait[i]
    jth_fold = accuracy_table_four$jth_fold[i]
    selected = four_folds[[jth_fold]][!is.na(Phenotyping_Wheat[four_folds[[jth_fold]], 
                                                               trait])]
    
    Y <- data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
                    obs=Phenotyping_Wheat[selected, trait])
    kin = kinship[Y$ID, Y$ID]
    rownames(kin) = Y$ID
    colnames(kin) = Y$ID
    
    gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                         Genotyping_Wheat[-selected, ])
    Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
    
    accuracy_table_four$g_acc_pear[i] = cor(Y$obs, Y$pred)
    accuracy_table_four$g_acc_gcor[i] = 
      estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
    if (grepl("IRR", trait)){
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_four$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_four$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_four$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_four$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    } else {
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_four$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_four$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_four$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_four$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    }
  }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
}

accuracy_table_three <- data.frame(fold="three folds", 
                                  trait=rep(colnames(Phenotyping_Wheat), each=3), 
                                  jth_fold=rep(1:3, length(colnames(Phenotyping_Wheat))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear_grain=NA, 
                                  p_acc_pear_leaf=NA, 
                                  p_acc_gcor_grain=NA, 
                                  p_acc_gcor_leaf=NA)
for (i in 1:nrow(accuracy_table_three)){
  tryCatch({
    trait = accuracy_table_three$trait[i]
    jth_fold = accuracy_table_three$jth_fold[i]
    selected = three_folds[[jth_fold]][!is.na(Phenotyping_Wheat[three_folds[[jth_fold]], 
                                                               trait])]
    
    Y <- data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
                    obs=Phenotyping_Wheat[selected, trait])
    kin = kinship[Y$ID, Y$ID]
    rownames(kin) = Y$ID
    colnames(kin) = Y$ID
    
    gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                         Genotyping_Wheat[-selected, ])
    Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
    
    accuracy_table_three$g_acc_pear[i] = cor(Y$obs, Y$pred)
    accuracy_table_three$g_acc_gcor[i] = 
      estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
    if (grepl("IRR", trait)){
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_three$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_three$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_three$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_three$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    } else {
      pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                 NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
      accuracy_table_three$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
      accuracy_table_three$p_acc_gcor_grain[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
                                NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
      Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
      accuracy_table_three$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
      accuracy_table_three$p_acc_gcor_leaf[i] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
    }
  }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
}

saveRDS(accuracy_table_five, "estimate_accuracy2/accuracy_table_five.rds")
saveRDS(accuracy_table_four, "estimate_accuracy2/accuracy_table_four.rds")
saveRDS(accuracy_table_three, "estimate_accuracy2/accuracy_table_three.rds")



accuracy_table <- rbind(rbind(accuracy_table_five, accuracy_table_four), accuracy_table_three)
accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=6), 1:3]
accuracy_table_long$type <- rep(colnames(accuracy_table)[4:9], nrow(accuracy_table))
accuracy_table_long$accuracy <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*6+1) : ((i-1)*6+6)] = 
    accuracy_table[i, 4:9]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
accuracy_table_long <- na.omit(accuracy_table_long)



traits <- unique(accuracy_table_long$trait)

pdf(paste("estimate_accuracy2/plots/", "accuracy.pdf", sep=""), width=14)
for (i in 1:length(traits)){
  print(ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i], ], 
               aes(x=type, y=accuracy, fill=type)) + 
          geom_boxplot() + 
          facet_wrap(~fold) + 
          scale_fill_hue(labels=c("g_acc_pear"="genomic accuracy Pearson", 
                                  "g_acc_gcor"="genomic accuracy adjusted", 
                                  "p_acc_pear_grain"="phenomic accuracy Pearson, Grain NIRS", 
                                  "p_acc_pear_leaf"="phenomic accuracy Pearson, Leaf NIRS", 
                                  "p_acc_gcor_grain"="phenomic accuracy adjusted, Grain NIRS", 
                                  "p_acc_gcor_leaf"="phenomic accuracy adjusted, Leaf NIRS")) + 
          theme_minimal_grid(font_size=10) + 
          theme(axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(), 
                legend.position="bottom") + 
          ggtitle(paste("trait: ", traits[i], sep=""))
  )
}
dev.off()





