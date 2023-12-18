setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)



file.names <- list.files("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", ".txt")
file.names.poplar <- gsub(".txt", "", file.names[grepl("Poplar", file.names, fixed=T)], )
file.names.wheat <- gsub(".txt", "", file.names[grepl("Wheat", file.names, fixed=T)], )

for (i in 1:length(file.names.poplar)){
  eval(call("<-", as.name(file.names.poplar[i]), 
            read.table(paste("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", 
                             file.names.poplar[i], 
                             ".txt", sep=""), 
                       header=T)))
}
# for (i in 1:length(file.names.wheat)){
#   eval(call("<-", as.name(file.names.wheat[i]), 
#             read.table(paste("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", 
#                              file.names.wheat[i], 
#                              ".txt", sep=""), 
#                        header=T)))
# }

ls()
# [1] "estimate_gcor"                "file.names"                  
# [3] "file.names.poplar"            "file.names.wheat"            
# [5] "Genotyping_Poplar"            "i"                           
# [7] "NIRS_NormDer_Wood_Poplar_ORL" "NIRS_NormDer_Wood_Poplar_SAV"
# [9] "Phenotyping_Poplar"   

rownames(Phenotyping_Poplar) <- Phenotyping_Poplar$Accession
Phenotyping_Poplar <- Phenotyping_Poplar[, -1]
rownames(Genotyping_Poplar) <- Genotyping_Poplar$Accession
Genotyping_Poplar <- Genotyping_Poplar[, -1]
rownames(NIRS_NormDer_Wood_Poplar_ORL) <- NIRS_NormDer_Wood_Poplar_ORL$Accession
NIRS_NormDer_Wood_Poplar_ORL <- NIRS_NormDer_Wood_Poplar_ORL[, -1]
rownames(NIRS_NormDer_Wood_Poplar_SAV) <- NIRS_NormDer_Wood_Poplar_SAV$Accession
NIRS_NormDer_Wood_Poplar_SAV <- NIRS_NormDer_Wood_Poplar_SAV[, -1]

Genotyping_Poplar <- scale(Genotyping_Poplar)
NIRS_NormDer_Wood_Poplar_ORL <- scale(NIRS_NormDer_Wood_Poplar_ORL)
NIRS_NormDer_Wood_Poplar_SAV <- scale(NIRS_NormDer_Wood_Poplar_SAV)

str(Phenotyping_Poplar)
dim(Genotyping_Poplar)
Genotyping_Poplar[1:5, 1:5]
dim(NIRS_NormDer_Wood_Poplar_ORL)
NIRS_NormDer_Wood_Poplar_ORL[1:5, 1:5]
NIRS_NormDer_Wood_Poplar_SAV[1:5, 1:5]

kinship <- Genotyping_Poplar %*% t(Genotyping_Poplar) / nrow(Genotyping_Poplar) #####

set.seed(1) 
ORL_names <- colnames(Phenotyping_Poplar)[grepl("ORL", colnames(Phenotyping_Poplar))]
SAV_names <- colnames(Phenotyping_Poplar)[!grepl("ORL", colnames(Phenotyping_Poplar))]
five_folds <- sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1]) ### 3, 4, 5 folds
five_folds <- split(five_folds, 1:5)
four_folds <- sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1])
four_folds <- split(four_folds, 1:4)
three_folds <- sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1])
three_folds <- split(three_folds, 1:3)

accuracy_table_five <- data.frame(fold="five", 
                             trait=rep(colnames(Phenotyping_Poplar), each=5), 
                             jth_fold=rep(1:5, length(colnames(Phenotyping_Poplar))), 
                             g_acc_pear=NA, 
                             g_acc_gcor=NA, 
                             p_acc_pear=NA, 
                             p_acc_gcor=NA)
i
# [1] 2
j
# [1] 4
for (i in 1:ncol(Phenotyping_Poplar)){
  for (j in 1:5){
    tryCatch({
      Y <- data.frame(ID=rownames(Genotyping_Poplar[five_folds[[j]], ]), 
                      obs=Phenotyping_Poplar[five_folds[[j]], colnames(Phenotyping_Poplar)[i]])
      kin = kinship[Y$ID, Y$ID]
      rownames(kin) = Y$ID
      colnames(kin) = Y$ID
      
      gmodel = mixed.solve(Phenotyping_Poplar[-five_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                           Genotyping_Poplar[-five_folds[[j]], ])
      Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
      
      accuracy_table_five$g_acc_pear[((i-1)*5 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_five$g_acc_gcor[((i-1)*5 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
        pmodel = mixed.solve(Phenotyping_Poplar[-five_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_ORL[-five_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_ORL[five_folds[[j]], ] %*% pmodel$u
      } else {
        pmodel = mixed.solve(Phenotyping_Poplar[-five_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_SAV[-five_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_SAV[five_folds[[j]], ] %*% pmodel$u
      }
      
      accuracy_table_five$p_acc_pear[((i-1)*5 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_five$p_acc_gcor[((i-1)*5 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
    }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
    }
}

accuracy_table_four <- data.frame(fold="four", 
                                  trait=rep(colnames(Phenotyping_Poplar), each=4), 
                                  jth_fold=rep(1:4, length(colnames(Phenotyping_Poplar))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear=NA, 
                                  p_acc_gcor=NA)
for (i in 1:ncol(Phenotyping_Poplar)){
  for (j in 1:4){
    tryCatch({
      Y <- data.frame(ID=rownames(Genotyping_Poplar[four_folds[[j]], ]), 
                      obs=Phenotyping_Poplar[four_folds[[j]], colnames(Phenotyping_Poplar)[i]])
      kin = kinship[Y$ID, Y$ID]
      rownames(kin) = Y$ID
      colnames(kin) = Y$ID
      
      gmodel = mixed.solve(Phenotyping_Poplar[-four_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                           Genotyping_Poplar[-four_folds[[j]], ])
      Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
      
      accuracy_table_four$g_acc_pear[((i-1)*4 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_four$g_acc_gcor[((i-1)*4 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
        pmodel = mixed.solve(Phenotyping_Poplar[-four_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_ORL[-four_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_ORL[four_folds[[j]], ] %*% pmodel$u
      } else {
        pmodel = mixed.solve(Phenotyping_Poplar[-four_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_SAV[-four_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_SAV[four_folds[[j]], ] %*% pmodel$u
      }
      
      accuracy_table_four$p_acc_pear[((i-1)*4 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_four$p_acc_gcor[((i-1)*4 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
    }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
  }
}

accuracy_table_three <- data.frame(fold="three", 
                                  trait=rep(colnames(Phenotyping_Poplar), each=3), 
                                  jth_fold=rep(1:3, length(colnames(Phenotyping_Poplar))), 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear=NA, 
                                  p_acc_gcor=NA)
for (i in 1:ncol(Phenotyping_Poplar)){
  for (j in 1:3){
    tryCatch({
      Y <- data.frame(ID=rownames(Genotyping_Poplar[three_folds[[j]], ]), 
                      obs=Phenotyping_Poplar[three_folds[[j]], colnames(Phenotyping_Poplar)[i]])
      kin = kinship[Y$ID, Y$ID]
      rownames(kin) = Y$ID
      colnames(kin) = Y$ID
      
      gmodel = mixed.solve(Phenotyping_Poplar[-three_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                           Genotyping_Poplar[-three_folds[[j]], ])
      Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
      
      accuracy_table_three$g_acc_pear[((i-1)*3 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_three$g_acc_gcor[((i-1)*3 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
      if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
        pmodel = mixed.solve(Phenotyping_Poplar[-three_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_ORL[-three_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_ORL[three_folds[[j]], ] %*% pmodel$u
      } else {
        pmodel = mixed.solve(Phenotyping_Poplar[-three_folds[[j]], colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_SAV[-three_folds[[j]], ])
        Y$pred = NIRS_NormDer_Wood_Poplar_SAV[three_folds[[j]], ] %*% pmodel$u
      }
      
      accuracy_table_three$p_acc_pear[((i-1)*3 + j)] = cor(Y$obs, Y$pred)
      accuracy_table_three$p_acc_gcor[((i-1)*3 + j)] = 
        estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)[1]
      
    }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
  }
}



saveRDS(accuracy_table_five, "estimate_accuracy/accuracy_table_five.rds")
saveRDS(accuracy_table_four, "estimate_accuracy/accuracy_table_four.rds")
saveRDS(accuracy_table_three, "estimate_accuracy/accuracy_table_three.rds")



accuracy_table <- rbind(rbind(accuracy_table_five, accuracy_table_four), accuracy_table_three)
accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=4), 1:3]
accuracy_table_long$type <- rep(colnames(accuracy_table)[4:7], nrow(accuracy_table))
accuracy_table_long$accuracy <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*4+1) : ((i-1)*4+4)] = 
    accuracy_table[i, 4:7]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
accuracy_table_long <- na.omit(accuracy_table_long)



p1 <- ggplot(accuracy_table_long, aes(x=type, y=accuracy)) + 
  geom_boxplot() + 
  theme_minimal_grid(font_size=10) 

save_plot(paste("estimate_accuracy/plots/", "accuracy.pdf", sep=""), 
          p1, 
          base_width=6.5)

  
  