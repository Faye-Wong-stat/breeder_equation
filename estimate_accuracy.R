setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)



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
for (i in 1:length(file.names.wheat)){
  eval(call("<-", as.name(file.names.wheat[i]), 
            read.table(paste("data/phenomic_selection_is_a_low_cost_and_high_throughput_method/", 
                             file.names.wheat[i], 
                             ".txt", sep=""), 
                       header=T)))
}

ls()
# [1] "estimate_gcor"                "file.names"                  
# [3] "file.names.poplar"            "Genotyping_Poplar"           
# [5] "i"                            "NIRS_NormDer_Wood_Poplar_ORL"
# [7] "NIRS_NormDer_Wood_Poplar_SAV" "Phenotyping_Poplar"  

rownames(Phenotyping_Poplar) <- Phenotyping_Poplar$Accession
Phenotyping_Poplar <- Phenotyping_Poplar[, -1]
rownames(Genotyping_Poplar) <- Genotyping_Poplar$Accession
Genotyping_Poplar <- Genotyping_Poplar[, -1]
rownames(NIRS_NormDer_Wood_Poplar_ORL) <- NIRS_NormDer_Wood_Poplar_ORL$Accession
NIRS_NormDer_Wood_Poplar_ORL <- NIRS_NormDer_Wood_Poplar_ORL[, -1]

Genotyping_Poplar <- scale(Genotyping_Poplar)
NIRS_NormDer_Wood_Poplar_ORL <- scale(NIRS_NormDer_Wood_Poplar_ORL)

str(Phenotyping_Poplar)
dim(Genotyping_Poplar)
Genotyping_Poplar[1:5, 1:5]
dim(NIRS_NormDer_Wood_Poplar_ORL)
NIRS_NormDer_Wood_Poplar_ORL[1:5, 1:5]

kinship <- Genotyping_Poplar %*% t(Genotyping_Poplar) / nrow(Genotyping_Poplar) #####

set.seed(1) #####
eight_fold <- sample(1:dim(Phenotyping_Poplar)[1], round(dim(Phenotyping_Poplar)[1] / 8, 0))

gmodel <- mixed.solve(Phenotyping_Poplar$HT.ORL[-eight_fold], 
                      Genotyping_Poplar[-eight_fold, ])
pmodel <- mixed.solve(Phenotyping_Poplar$HT.ORL[-eight_fold], 
                      NIRS_NormDer_Wood_Poplar_ORL[-eight_fold, ])
  
pred.gmodel <- Genotyping_Poplar[eight_fold, ] %*% gmodel$u
pred.pmodel <- NIRS_NormDer_Wood_Poplar_ORL[eight_fold, ] %*% pmodel$u
  
Y <- data.frame(ID=rownames(Genotyping_Poplar[eight_fold, ]), 
                obs=Phenotyping_Poplar$HT.ORL[eight_fold], 
                pred=pred.gmodel)
kin <- matrix(NA, nrow=nrow(Y), ncol=nrow(Y))
rownames(kin) <- Y$ID
colnames(kin) <- Y$ID
for (i in 1:nrow(kin)){
  for (j in 1:ncol(kin)){
    kin[i, j] = kinship[Y$ID[i], Y$ID[j]]
      Y$ID[i]
  }
}

A <- estimate_gcor(data=Y, Knn=kin, method="MCMCglmm", normalize=F)
A
# g_cor      Rhat 
# 0.3582751 1.0133457 





  
  