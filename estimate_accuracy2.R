# Renaud Rincent on wheat
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

# Genotyping_Wheat <- scale(Genotyping_Wheat)
Genotyping_Wheat <- as.matrix(Genotyping_Wheat)
Genotyping_Wheat <- Genotyping_Wheat*2 - 1
NIRS_NormDer_Grain_DRY_Wheat <- scale(NIRS_NormDer_Grain_DRY_Wheat)
NIRS_NormDer_Grain_IRR_Wheat <- scale(NIRS_NormDer_Grain_IRR_Wheat)
NIRS_NormDer_Leaf_DRY_Wheat <- scale(NIRS_NormDer_Leaf_DRY_Wheat)
NIRS_NormDer_Leaf_IRR_Wheat <- scale(NIRS_NormDer_Leaf_IRR_Wheat)
str(Phenotyping_Wheat)

kinship <- A.mat(Genotyping_Wheat)

set.seed(1) 
eight_folds <- vector(mode="list", length=ncol(Phenotyping_Wheat))
eight_folds <- lapply(eight_folds, FUN=function(x){
  vector(mode="list", length=20)
})
for (i in 1:ncol(Phenotyping_Wheat)){
  eight_folds[[i]] = lapply(eight_folds[[i]], FUN=function(x){
    sample(1:dim(Genotyping_Wheat)[1], round(dim(Genotyping_Wheat)[1]/8))
  })
}


# eight_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
# eight_folds <- split(eight_folds, 1:5)
# four_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
# four_folds <- split(four_folds, 1:4)
# three_folds <- sample(1:nrow(Genotyping_Wheat), nrow(Genotyping_Wheat))
# three_folds <- split(three_folds, 1:3)



accuracy_table_eight <- data.frame(
                                  trait=rep(colnames(Phenotyping_Wheat), each=20), 
                                  jth_fold=rep(1:20, length(colnames(Phenotyping_Wheat))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear_grain=NA, 
                                  p_acc_pear_leaf=NA, 
                                  p_acc_gcor_grain=NA, 
                                  p_acc_gcor_leaf=NA, 
                                  h2=NA)
for (i in 1:ncol(Phenotyping_Wheat)){
  for (j in 1:20){
    tryCatch({
      selected = eight_folds[[i]][[j]][!is.na(Phenotyping_Wheat[eight_folds[[i]][[j]], i])]
      
      Y = data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
                     obs=Phenotyping_Wheat[selected, i])
      
      kin = kinship[Y$ID, Y$ID]
      scale_by = mean(diag(kin))
      K = kinship/scale_by
      
      K_train = K[rownames(Genotyping_Wheat[-selected, ]), rownames(Genotyping_Wheat[-selected, ])]
      K_test = K[Y$ID, Y$ID]
      
      gmodel = mixed.solve(Phenotyping_Wheat[-selected, i], 
                           Genotyping_Wheat[-selected, ])
      Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
      
      accuracy_table_eight$g_acc_pear[(i-1)*20 + j] = cor(Y$obs, Y$pred)
      
      gmodel2 = mixed.solve(Phenotyping_Wheat[selected, i], K=K_test)
      h2 = gmodel2$Vu / var(Y$obs)
      accuracy_table_eight$g_acc_gcor[(i-1)*20 + j] = cor(Y$obs, Y$pred) / sqrt(h2)
      accuracy_table_eight$h2[(i-1)*20 + j] = h2
      
      
      if (grepl("IRR", colnames(Phenotyping_Wheat)[i])){
        pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, i], 
                                   NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
        Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
        accuracy_table_eight$p_acc_pear_grain[(i-1)*20 + j] = cor(Y$obs, Y$pred)
        accuracy_table_eight$p_acc_gcor_grain[(i-1)*20 + j] = 
          estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
        
        pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, i], 
                                  NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
        Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
        accuracy_table_eight$p_acc_pear_leaf[(i-1)*20 + j] = cor(Y$obs, Y$pred)
        accuracy_table_eight$p_acc_gcor_leaf[(i-1)*20 + j] = 
          estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
      } else {
        pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, i], 
                                   NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
        Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
        accuracy_table_eight$p_acc_pear_grain[(i-1)*20 + j] = cor(Y$obs, Y$pred)
        accuracy_table_eight$p_acc_gcor_grain[(i-1)*20 + j] = 
          estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
        
        pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, i], 
                                  NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
        Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
        accuracy_table_eight$p_acc_pear_leaf[(i-1)*20 + j] = cor(Y$obs, Y$pred)
        accuracy_table_eight$p_acc_gcor_leaf[(i-1)*20 + j] = 
          estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
      }
      
    }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
  }
}
# 
# for (i in 1:nrow(accuracy_table_eight)){
#   tryCatch({
#     trait = accuracy_table_eight$trait[i]
#     jth_fold = accuracy_table_eight$jth_fold[i]
#     selected = eight_folds[[jth_fold]][!is.na(Phenotyping_Wheat[eight_folds[[jth_fold]], 
#                                                             trait])]
#     
#     Y = data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
#                     obs=Phenotyping_Wheat[selected, trait])
#     
#     kin = kinship[Y$ID, Y$ID]
#     scale_by = mean(diag(kin))
#     K = kinship/scale_by
#     
#     K_train = K[rownames(Genotyping_Wheat[-selected, ]), rownames(Genotyping_Wheat[-selected, ])]
#     K_test = K[Y$ID, Y$ID]
#     
#     gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                          Genotyping_Wheat[-selected, ])
#     Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
#     
#     accuracy_table_eight$g_acc_pear[i] = cor(Y$obs, Y$pred)
#     
#     # K = kinship[-selected, -selected]
#     gmodel2 = mixed.solve(Phenotyping_Wheat[selected, trait], K=K_test)
#     h2 = gmodel2$Vu / var(Y$obs)
#     accuracy_table_eight$g_acc_gcor[i] = cor(Y$obs, Y$pred) / sqrt(h2)
#     accuracy_table_eight$h2[i] = h2
#     
#     if (grepl("IRR", trait)){
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_eight$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_eight$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_eight$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_eight$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     } else {
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_eight$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_eight$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_eight$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_eight$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     }
#   }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
# }
# 
# accuracy_table_four <- data.frame(fold="four folds", 
#                                   trait=rep(colnames(Phenotyping_Wheat), each=4), 
#                                   jth_fold=rep(1:4, length(colnames(Phenotyping_Wheat))), 
#                                   g_acc_pear=NA, 
#                                   g_acc_gcor=NA, 
#                                   p_acc_pear_grain=NA, 
#                                   p_acc_pear_leaf=NA, 
#                                   p_acc_gcor_grain=NA, 
#                                   p_acc_gcor_leaf=NA)
# for (i in 1:nrow(accuracy_table_four)){
#   tryCatch({
#     trait = accuracy_table_four$trait[i]
#     jth_fold = accuracy_table_four$jth_fold[i]
#     selected = four_folds[[jth_fold]][!is.na(Phenotyping_Wheat[four_folds[[jth_fold]], 
#                                                                trait])]
#     
#     Y = data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
#                     obs=Phenotyping_Wheat[selected, trait])
#     
#     kin = kinship[Y$ID, Y$ID]
#     scale_by = mean(diag(kin))
#     K = kinship/scale_by
#     
#     K_train = K[rownames(Genotyping_Wheat[-selected, ]), rownames(Genotyping_Wheat[-selected, ])]
#     K_test = K[Y$ID, Y$ID]
#     
#     gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                          Genotyping_Wheat[-selected, ])
#     Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
#     
#     accuracy_table_four$g_acc_pear[i] = cor(Y$obs, Y$pred)
#     
#     # K = kinship[-selected, -selected]
#     gmodel2 = mixed.solve(Phenotyping_Wheat[selected, trait], K=K_test)
#     h2 = gmodel2$Vu / var(Y$obs)
#     accuracy_table_four$g_acc_gcor[i] = accuracy_table_four$g_acc_pear[i] / sqrt(h2)
#     accuracy_table_four$h2[i] = h2
#     
#     if (grepl("IRR", trait)){
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_four$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_four$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_four$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_four$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     } else {
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_four$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_four$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_four$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_four$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     }
#   }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
# }
# 
# accuracy_table_three <- data.frame(fold="three folds", 
#                                   trait=rep(colnames(Phenotyping_Wheat), each=3), 
#                                   jth_fold=rep(1:3, length(colnames(Phenotyping_Wheat))), 
#                                   g_acc_pear=NA, 
#                                   g_acc_gcor=NA, 
#                                   p_acc_pear_grain=NA, 
#                                   p_acc_pear_leaf=NA, 
#                                   p_acc_gcor_grain=NA, 
#                                   p_acc_gcor_leaf=NA)
# for (i in 1:nrow(accuracy_table_three)){
#   tryCatch({
#     trait = accuracy_table_three$trait[i]
#     jth_fold = accuracy_table_three$jth_fold[i]
#     selected = three_folds[[jth_fold]][!is.na(Phenotyping_Wheat[three_folds[[jth_fold]], 
#                                                                trait])]
#     
#     Y = data.frame(ID=rownames(Genotyping_Wheat[selected, ]), 
#                     obs=Phenotyping_Wheat[selected, trait])
#     
#     kin = kinship[Y$ID, Y$ID]
#     scale_by = mean(diag(kin))
#     K = kinship/scale_by
#     
#     K_train = K[rownames(Genotyping_Wheat[-selected, ]), rownames(Genotyping_Wheat[-selected, ])]
#     K_test = K[Y$ID, Y$ID]
#     
#     gmodel = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                          Genotyping_Wheat[-selected, ])
#     Y$pred = Genotyping_Wheat[Y$ID, ] %*% gmodel$u
#     
#     accuracy_table_three$g_acc_pear[i] = cor(Y$obs, Y$pred)
#     
#     # K = kinship[-selected, -selected]
#     gmodel2 = mixed.solve(Phenotyping_Wheat[selected, trait], K=K_test)
#     h2 = gmodel2$Vu / var(Y$obs)
#     accuracy_table_three$g_acc_gcor[i] = accuracy_table_three$g_acc_pear[i] / sqrt(h2)
#     accuracy_table_three$h2[i] = h2
#     
#     if (grepl("IRR", trait)){
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_IRR_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_three$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_three$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_IRR_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_IRR_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_three$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_three$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     } else {
#       pmodel_grain = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                  NIRS_NormDer_Grain_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Grain_DRY_Wheat[selected, ] %*% pmodel_grain$u
#       accuracy_table_three$p_acc_pear_grain[i] = cor(Y$obs, Y$pred)
#       accuracy_table_three$p_acc_gcor_grain[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#       pmodel_leaf = mixed.solve(Phenotyping_Wheat[-selected, trait], 
#                                 NIRS_NormDer_Leaf_DRY_Wheat[-selected, ])
#       Y$pred = NIRS_NormDer_Leaf_DRY_Wheat[selected, ] %*% pmodel_leaf$u
#       accuracy_table_three$p_acc_pear_leaf[i] = cor(Y$obs, Y$pred)
#       accuracy_table_three$p_acc_gcor_leaf[i] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#     }
#   }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
# }

saveRDS(accuracy_table_eight, "estimate_accuracy2/accuracy_table_eight.rds")
# saveRDS(accuracy_table_four, "estimate_accuracy2/accuracy_table_four.rds")
# saveRDS(accuracy_table_three, "estimate_accuracy2/accuracy_table_three.rds")
accuracy_table_eight <- readRDS("estimate_accuracy2/accuracy_table_eight.rds")
# accuracy_table_four <- readRDS("estimate_accuracy2/accuracy_table_four.rds")
# accuracy_table_three <- readRDS("estimate_accuracy2/accuracy_table_three.rds")



accuracy_table <- accuracy_table_eight
accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=6), 1:2]
accuracy_table_long$type <- rep(colnames(accuracy_table)[3:8], nrow(accuracy_table))
accuracy_table_long$accuracy <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*6+1) : ((i-1)*6+6)] = 
    accuracy_table[i, 3:8]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
h2_table <- aggregate(accuracy_table$h2, list(accuracy_table$trait), mean)
rownames(h2_table) <- h2_table$Group.1
accuracy_table_long$h2 <- h2_table[accuracy_table_long$trait, "x"]
accuracy_table_long <- na.omit(accuracy_table_long)
accuracy_table_long$h2 <- as.character(round(accuracy_table_long$h2, 5))
accuracy_table_long$type <- factor(accuracy_table_long$type, levels=c("g_acc_pear", "p_acc_pear_grain", 
                                                                      "p_acc_pear_leaf", "g_acc_gcor", 
                                                                      "p_acc_gcor_grain", "p_acc_gcor_leaf"))
accuracy_table_long$adjustment <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "p_acc_pear_grain", 
                                                                         "p_acc_pear_leaf"), 
                                         "before correction", "after correction")
accuracy_table_long$adjustment <- factor(accuracy_table_long$adjustment, 
                                         levels=c("before correction", "after correction"))
accuracy_table_long$`prediction method` <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "g_acc_gcor"), 
                                         "genomic", ifelse(accuracy_table_long$type %in% 
                                                             c("p_acc_pear_grain", "p_acc_gcor_grain"), 
                                                           "grain-phenomic", "leaf-phenomic"))
accuracy_table_long$`prediction method` <- factor(accuracy_table_long$`prediction method`, 
                                         levels = c("genomic", "grain-phenomic", "leaf-phenomic"))
accuracy_table_long$trait2 <- ifelse(accuracy_table_long$trait %in% c("IRR_GrainYield", "DRY_GrainYield"), 
                                     "GrainYield", "HeadingDate")

saveRDS(accuracy_table_long, "estimate_accuracy2/accuracy_table_long.rds")
accuracy_table_long <- readRDS("estimate_accuracy2/accuracy_table_long.rds")



traits <- unique(accuracy_table_long$trait)
h2 <- unique(accuracy_table_long$h2)


for (i in 1:length(traits)){
  p1 = ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i] & 
                                    accuracy_table_long$adjustment == "before correction", ],
              aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) +
    geom_boxplot() +
    facet_wrap(~adjustment) +
    ylim(0, 1) + 
    ylab("predictive ability") + 
    theme_minimal_grid(font_size=12) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="bottom") + 
    ggtitle(paste("trait: ", traits[i], ", h2: ", h2[i], sep=""))
  
  p2 = ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i] & 
                                    accuracy_table_long$adjustment == "after correction", ],
              aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) +
    geom_boxplot() +
    facet_wrap(~adjustment) +
    ylim(0, 1) + 
    ylab("prediction accuracy") + 
    theme_minimal_grid(font_size=12) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position="none") +
    ggtitle(paste("trait: ", traits[i], ", h2: ", h2[i], sep=""))
  
  save_plot(paste("estimate_accuracy2/plots/", traits[i], "_accuracy.pdf", sep=""), 
            plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))
}

pdf(paste("estimate_accuracy2/plots/", "accuracy.pdf", sep=""), width=14, height=6)
for (i in 1:length(traits)){
  print(
    ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i] , ],
           aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) +
      geom_boxplot() +
      facet_wrap(~adjustment) +
      ylim(0, 1) + 
      ylab("predictive ability/prediction accuracy") + 
      theme_minimal_grid(font_size=12) +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position="bottom") +
      ggtitle(paste("trait: ", traits[i], ", h2: ", h2[i], sep=""))
  )
  
  # print(ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i], ], 
  #              aes(x=type, y=accuracy, fill=type)) + 
  #         geom_boxplot() + 
  #         facet_wrap(~fold) + 
  #         scale_fill_hue(labels=c("g_acc_pear"="genomic accuracy Pearson", 
  #                                 "g_acc_gcor"="genomic accuracy adjusted", 
  #                                 "p_acc_pear_grain"="phenomic accuracy Pearson, Grain NIRS", 
  #                                 "p_acc_pear_leaf"="phenomic accuracy Pearson, Leaf NIRS", 
  #                                 "p_acc_gcor_grain"="phenomic accuracy adjusted, Grain NIRS", 
  #                                 "p_acc_gcor_leaf"="phenomic accuracy adjusted, Leaf NIRS")) + 
  #         theme_minimal_grid(font_size=10) + 
  #         theme(axis.title.x = element_blank(),
  #               axis.text.x = element_blank(),
  #               axis.ticks.x = element_blank(), 
  #               legend.position="bottom") + 
  #         ggtitle(paste("trait: ", traits[i], ", h2: ", h2[i], sep=""))
  # )
}
dev.off()





