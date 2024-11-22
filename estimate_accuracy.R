# Renaud Rincent on poplar
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
ls()
# [1] "estimate_gcor"                "file.names"                  
# [3] "file.names.poplar"            "file.names.wheat"            
# [5] "Genotyping_Poplar"            "i"                           
# [7] "NIRS_NormDer_Wood_Poplar_ORL" "NIRS_NormDer_Wood_Poplar_SAV"
# [9] "Phenotyping_Poplar"   

# for (i in 1:length(file.names.wheat)){
#   eval(call("<-", as.name(file.names.wheat[i]), 
#             read.table(paste("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", 
#                              file.names.wheat[i], 
#                              ".txt", sep=""), 
#                        header=T)))
# }


rownames(Phenotyping_Poplar) <- Phenotyping_Poplar$Accession
Phenotyping_Poplar <- Phenotyping_Poplar[, -1]
rownames(Genotyping_Poplar) <- Genotyping_Poplar$Accession
Genotyping_Poplar <- Genotyping_Poplar[, -1]
rownames(NIRS_NormDer_Wood_Poplar_ORL) <- NIRS_NormDer_Wood_Poplar_ORL$Accession
NIRS_NormDer_Wood_Poplar_ORL <- NIRS_NormDer_Wood_Poplar_ORL[, -1]
rownames(NIRS_NormDer_Wood_Poplar_SAV) <- NIRS_NormDer_Wood_Poplar_SAV$Accession
NIRS_NormDer_Wood_Poplar_SAV <- NIRS_NormDer_Wood_Poplar_SAV[, -1]

# Genotyping_Poplar <- scale(Genotyping_Poplar)
Genotyping_Poplar <- Genotyping_Poplar*2 - 1
Genotyping_Poplar <- as.matrix(Genotyping_Poplar)
NIRS_NormDer_Wood_Poplar_ORL <- scale(NIRS_NormDer_Wood_Poplar_ORL)
NIRS_NormDer_Wood_Poplar_SAV <- scale(NIRS_NormDer_Wood_Poplar_SAV)

str(Phenotyping_Poplar)
dim(Genotyping_Poplar)
Genotyping_Poplar[1:5, 1:5]
dim(NIRS_NormDer_Wood_Poplar_ORL)
NIRS_NormDer_Wood_Poplar_ORL[1:5, 1:5]
NIRS_NormDer_Wood_Poplar_SAV[1:5, 1:5]

kinship <- A.mat(Genotyping_Poplar)

set.seed(1) 
ORL_names <- colnames(Phenotyping_Poplar)[grepl("ORL", colnames(Phenotyping_Poplar))]
SAV_names <- colnames(Phenotyping_Poplar)[!grepl("ORL", colnames(Phenotyping_Poplar))]
five_folds <- vector(mode="list", length=ncol(Phenotyping_Poplar))
five_folds <- lapply(five_folds, FUN=function(x){
  vector(mode="list", length=20)
})
for (i in 1:ncol(Phenotyping_Poplar)){
  five_folds[[i]] = lapply(five_folds[[i]], FUN=function(x){
    sample(1:dim(Phenotyping_Poplar)[1], round(dim(Phenotyping_Poplar)[1]/5))
  })
}

#   sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1]) ### 3, 4, 5 folds
# five_folds <- split(five_folds, 1:5)
# four_folds <- sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1])
# four_folds <- split(four_folds, 1:4)
# three_folds <- sample(1:dim(Phenotyping_Poplar)[1], dim(Phenotyping_Poplar)[1])
# three_folds <- split(three_folds, 1:3)

accuracy_table_five <- data.frame(fold="five folds", 
                             trait=rep(colnames(Phenotyping_Poplar), each=20), 
                             jth_fold=rep(1:20, length(colnames(Phenotyping_Poplar))), 
                             g_acc_pear=NA, 
                             g_acc_gcor=NA, 
                             p_acc_pear=NA, 
                             p_acc_gcor=NA, 
                             h2 = NA)

for (i in 1:ncol(Phenotyping_Poplar)){
  for (j in 1:20){
    tryCatch({
      selected = five_folds[[i]][[j]][!is.na(Phenotyping_Poplar[five_folds[[i]][[j]], 
                                                                colnames(Phenotyping_Poplar)[i]])]
      
      Y = data.frame(ID=rownames(Genotyping_Poplar[selected, ]), 
                      obs=Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]])
      
      kin = kinship[Y$ID, Y$ID]
      scale_by = mean(diag(kin))
      K = kinship/scale_by
      
      K_train = K[rownames(Genotyping_Poplar[-selected, ]), rownames(Genotyping_Poplar[-selected, ])]
      K_test = K[Y$ID, Y$ID]
      
      gmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
                           Genotyping_Poplar[-selected, ])
      Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
      
      # pdf("estimate_accuracy/plots/marker_effects.pdf")
      # hist(gmodel$u)
      # dev.off()
      # 
      accuracy_table_five$g_acc_pear[((i-1)*20 + j)] = cor(Y$obs, Y$pred)
      
      gmodel2 = mixed.solve(Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]], K=K_test)
      h2 = gmodel2$Vu / var(Y$obs)
      accuracy_table_five$g_acc_gcor[((i-1)*20 + j)] = cor(Y$obs, Y$pred) / sqrt(h2)
      accuracy_table_five$h2[((i-1)*20 + j)] = h2
      
      Y2 = Y
      if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
        pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_ORL[-selected, ])
        Y2$pred = NIRS_NormDer_Wood_Poplar_ORL[selected, ] %*% pmodel$u
      } else {
        pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
                             NIRS_NormDer_Wood_Poplar_SAV[-selected, ])
        Y2$pred = NIRS_NormDer_Wood_Poplar_SAV[selected, ] %*% pmodel$u
      }
      
      accuracy_table_five$p_acc_pear[((i-1)*20 + j)] = cor(Y2$obs, Y2$pred)
      accuracy_table_five$p_acc_gcor[((i-1)*20 + j)] = 
        estimate_gcor(data=Y2, Knn=K_test, method="MCMCglmm", normalize=F)[1]
      
    }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
    }
}
# saveRDS(Y, "Y.rds")
# saveRDS(Y2, "Y2.rds")
# saveRDS(kin, "kin.rds")
# 
# set.seed(11)
# X = sample(1:ncol(Genotyping_Poplar), 1000)
# test_table <- data.frame(geno=X, 
#                          Vu=NA,
#                          Ve=NA)
# for (k in 1:1000){
#   pheno = Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]]
#   geno = as.factor(Genotyping_Poplar[-selected, X[k]])
#   model = lmer(pheno ~ (1 | geno))
#   variance = data.frame(VarCorr(model))
#   test_table$Vu[k] = variance[1, "vcov"]
#   test_table$Ve[k] = variance[2, "vcov"] / mean(table(geno))
# }
# test_table$h2 <- test_table$Vu / (test_table$Vu + test_table$Ve)
# saveRDS(test_table, "estimate_accuracy/test_table.rds")
# pdf("estimate_accuracy/plots/h2.pdf")
# hist(test_table$h2)
# dev.off()


# 
# accuracy_table_four <- data.frame(fold="four folds", 
#                                   trait=rep(colnames(Phenotyping_Poplar), each=4), 
#                                   jth_fold=rep(1:4, length(colnames(Phenotyping_Poplar))), 
#                                   g_acc_pear=NA, 
#                                   g_acc_gcor=NA, 
#                                   p_acc_pear=NA, 
#                                   p_acc_gcor=NA, 
#                                   h2=NA)
# for (i in 1:ncol(Phenotyping_Poplar)){
#   for (j in 1:4){
#     tryCatch({
#       selected = four_folds[[j]][!is.na(Phenotyping_Poplar[four_folds[[j]], colnames(Phenotyping_Poplar)[i]])]
#       
#       Y = data.frame(ID=rownames(Genotyping_Poplar[selected, ]), 
#                       obs=Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]])
#       
#       kin = kinship[Y$ID, Y$ID]
#       scale_by = mean(diag(kin))
#       K = kinship/scale_by
#       
#       K_train = K[rownames(Genotyping_Poplar[-selected, ]), rownames(Genotyping_Poplar[-selected, ])]
#       K_test = K[Y$ID, Y$ID]
#       
#       gmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                            Genotyping_Poplar[-selected, ])
#       Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
#       
#       accuracy_table_four$g_acc_pear[((i-1)*4 + j)] = cor(Y$obs, Y$pred)
#       
#       # K = kinship[selected, selected]
#       gmodel2 = mixed.solve(Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]], K=K_test)
#       h2 = gmodel2$Vu / var(Y$obs)
#       accuracy_table_four$g_acc_gcor[((i-1)*4 + j)] = accuracy_table_four$g_acc_pear[((i-1)*4 + j)] / sqrt(h2)
#       accuracy_table_four$h2[((i-1)*4 + j)] = h2
#       
#       if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
#         pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                              NIRS_NormDer_Wood_Poplar_ORL[-selected, ])
#         Y$pred = NIRS_NormDer_Wood_Poplar_ORL[selected, ] %*% pmodel$u
#       } else {
#         pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                              NIRS_NormDer_Wood_Poplar_SAV[-selected, ])
#         Y$pred = NIRS_NormDer_Wood_Poplar_SAV[selected, ] %*% pmodel$u
#       }
#       
#       accuracy_table_four$p_acc_pear[((i-1)*4 + j)] = cor(Y$obs, Y$pred)
#       accuracy_table_four$p_acc_gcor[((i-1)*4 + j)] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#     }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
#   }
# }
# 
# accuracy_table_three <- data.frame(fold="three folds", 
#                                   trait=rep(colnames(Phenotyping_Poplar), each=3), 
#                                   jth_fold=rep(1:3, length(colnames(Phenotyping_Poplar))), 
#                                   g_acc_pear=NA, 
#                                   g_acc_gcor=NA, 
#                                   p_acc_pear=NA, 
#                                   p_acc_gcor=NA, 
#                                   h2=NA)
# for (i in 1:ncol(Phenotyping_Poplar)){
#   for (j in 1:3){
#     tryCatch({
#       selected = three_folds[[j]][!is.na(Phenotyping_Poplar[three_folds[[j]], 
#                                                             colnames(Phenotyping_Poplar)[i]])]
#       
#       Y = data.frame(ID=rownames(Genotyping_Poplar[selected, ]), 
#                       obs=Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]])
#       
#       kin = kinship[Y$ID, Y$ID]
#       scale_by = mean(diag(kin))
#       K = kinship/scale_by
#       
#       K_train = K[rownames(Genotyping_Poplar[-selected, ]), rownames(Genotyping_Poplar[-selected, ])]
#       K_test = K[Y$ID, Y$ID]
#       
#       gmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                            Genotyping_Poplar[-selected, ])
#       Y$pred = Genotyping_Poplar[Y$ID, ] %*% gmodel$u
#       
#       accuracy_table_three$g_acc_pear[((i-1)*3 + j)] = cor(Y$obs, Y$pred)
#       
#       # K = kinship[selected, selected]
#       gmodel2 = mixed.solve(Phenotyping_Poplar[selected, colnames(Phenotyping_Poplar)[i]], K=K_test)
#       h2 = gmodel2$Vu / var(Y$obs)
#       accuracy_table_three$g_acc_gcor[((i-1)*3 + j)] = accuracy_table_three$g_acc_pear[((i-1)*3 + j)] / sqrt(h2)
#       accuracy_table_three$h2[((i-1)*3 + j)] = h2
#       
#       if (colnames(Phenotyping_Poplar)[i] %in% ORL_names){
#         pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                              NIRS_NormDer_Wood_Poplar_ORL[-selected, ])
#         Y$pred = NIRS_NormDer_Wood_Poplar_ORL[selected, ] %*% pmodel$u
#       } else {
#         pmodel = mixed.solve(Phenotyping_Poplar[-selected, colnames(Phenotyping_Poplar)[i]], 
#                              NIRS_NormDer_Wood_Poplar_SAV[-selected, ])
#         Y$pred = NIRS_NormDer_Wood_Poplar_SAV[selected, ] %*% pmodel$u
#       }
#       
#       accuracy_table_three$p_acc_pear[((i-1)*3 + j)] = cor(Y$obs, Y$pred)
#       accuracy_table_three$p_acc_gcor[((i-1)*3 + j)] = 
#         estimate_gcor(data=Y, Knn=K_test, method="MCMCglmm", normalize=F)[1]
#       
#     }, error=function(e){cat(i, j, "\n", "Error:", conditionMessage(e), "\n")})
#   }
# }



saveRDS(accuracy_table_five, "estimate_accuracy/accuracy_table_five.rds")
# saveRDS(accuracy_table_four, "estimate_accuracy/accuracy_table_four.rds")
# saveRDS(accuracy_table_three, "estimate_accuracy/accuracy_table_three.rds")
accuracy_table_five <- readRDS("estimate_accuracy/accuracy_table_five.rds")
# accuracy_table_four <- readRDS("estimate_accuracy/accuracy_table_four.rds")
# accuracy_table_three <- readRDS("estimate_accuracy/accuracy_table_three.rds")



accuracy_table <- accuracy_table_five

accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=4), 2:3]
accuracy_table_long$type <- rep(colnames(accuracy_table)[4:7], nrow(accuracy_table))
accuracy_table_long$accuracy <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*4+1) : ((i-1)*4+4)] = 
    accuracy_table[i, 4:7]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
h2_table <- aggregate(accuracy_table$h2, list(accuracy_table$trait), mean)
rownames(h2_table) <- h2_table$Group.1
accuracy_table_long$h2 <- h2_table[accuracy_table_long$trait, "x"]
accuracy_table_long <- na.omit(accuracy_table_long)
accuracy_table_long$h2 <- as.character(round(accuracy_table_long$h2, 5))
accuracy_table_long$type <- factor(accuracy_table_long$type, levels=c("g_acc_pear", "p_acc_pear", 
                                                                      "g_acc_gcor", "p_acc_gcor"))

accuracy_table_long$adjustment <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "p_acc_pear"), 
                                         "before correction", "after correction")
accuracy_table_long$adjustment <- factor(accuracy_table_long$adjustment, 
                                         levels=c("before correction", "after correction"))
accuracy_table_long$`prediction method` <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "g_acc_gcor"), 
                                         "genomic", "phenomic")

saveRDS(accuracy_table_long, "estimate_accuracy/accuracy_table_long.rds")
accuracy_table_long <- readRDS("estimate_accuracy/accuracy_table_long.rds")



traits <- unique(accuracy_table_long$trait)
h2 <- unique(accuracy_table_long$h2)



# pdf(paste("estimate_accuracy/plots/", "accuracy.pdf", sep=""), width=10)
# for (i in 1:length(traits)){
#   print(ggplot(accuracy_table_long[accuracy_table_long$trait == traits[i], ], 
#                aes(x=type, y=accuracy, fill=type)) + 
#     geom_boxplot() + 
#     # facet_wrap(~fold) + 
#     scale_fill_hue(labels=c("g_acc_pear"="genomic accuracy Pearson", 
#                             "g_acc_gcor"="genomic accuracy adjusted", 
#                             "p_acc_pear"="phenomic accuracy Pearson", 
#                             "p_acc_gcor"="phenomic accuracy adjusted")) + 
#     ylim(0, 1) + 
#     ylab("predictive ability/prediction accuracy") + 
#     theme_minimal_grid(font_size=12) + 
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank(), 
#           legend.position="bottom") + 
#     ggtitle(paste("trait: ", traits[i], ", h2: ", h2[i], sep=""))
# )
#   
#   # save_plot(paste("estimate_accuracy/plots/", traits[i], "_accuracy.pdf", sep=""), 
#   #           p1, 
#   #           base_width=6.5)
#   
# }
# dev.off()





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
  
  save_plot(paste("estimate_accuracy/plots/", traits[i], "_accuracy.pdf", sep=""), 
            plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))
}








 