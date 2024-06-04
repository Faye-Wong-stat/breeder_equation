# Zachary Winn on wheat
setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(lme4)
library(emmeans)
library(rrBLUP)
library(ggplot2)
library(cowplot)



geno <- read.table("data/phenomic_versus_genomic_prediction/geno_coded.txt")
pheno <- read.table("data/phenomic_versus_genomic_prediction/pheno_coded.txt", header=T)



# check allele frequency and heterozygous rate
dim(geno)
# [1]  3442 23130

alle_freq <- apply(geno, 2, FUN=function(x){
  sum(x)/(length(x)*2)
})
minor_alle_freq <- ifelse(alle_freq < 0.5, alle_freq, 1-alle_freq)
sum(minor_alle_freq<0.05)
# [1] 5117

heter_freq <- apply(geno, 2, function(x){
  sum(x==1)/23130
})
sum(heter_freq>0.1)
# [1] 0

geno <- geno[, ifelse(minor_alle_freq<0.05, F, T)]
geno <- as.matrix(geno)

saveRDS(geno, "estimate_accuracy5/geno.rds")
# geno <- readRDS("estimate_accuracy5/geno.rds")

kinship <- A.mat(geno-1)
saveRDS(kinship, "estimate_accuracy5/kinship.rds")



# organize phenotyic/phenomic data
model <- lmer(grain_yield ~ fullsamplename + (1|year) + (1|nursery), pheno)
emm_summ <- summary(emmeans(model, "fullsamplename"))

phenotype <- data.frame(grain_yield = emm_summ[, 2][match(rownames(geno), emm_summ[, 1])])
rownames(phenotype) <- emm_summ[, 1][match(rownames(geno), emm_summ[, 1])]
phenotype <- as.matrix(phenotype)

saveRDS(phenotype, "estimate_accuracy5/phenotype.rds")
# phenotype <- readRDS("estimate_accuracy5/phenotype.rds")

models <- list()
for (i in colnames(pheno[, 21:44])){
  models[[i]] = lmer(get(i) ~ fullsamplename + (1|year), pheno)
}
emm_summs <- list()
for (i in 1:length(models)){
  emm_summs[[i]] = summary(emmeans(models[[i]], "fullsamplename"))
}

saveRDS(emm_summs, "estimate_accuracy5/emm_summs.rds")
# emm_summs <- readRDS("estimate_accuracy5/emm_summs.rds")

phenomic <- data.frame(matrix(NA, nrow=nrow(emm_summs[[1]]), ncol=length(emm_summs)+1))
phenomic[, 1:2]  <- emm_summs[[1]][, 1:2]
for (i in 2:length(emm_summs)){
  phenomic[, i+1] = emm_summs[[i]][, 2]
}

colnames(phenomic) <- c("fullsamplename", colnames(pheno)[21:44])
rownames(phenomic) <- phenomic$fullsamplename
phenomic <- phenomic[match(rownames(geno), rownames(phenomic)), ]
phenomic <- phenomic[, -1]
phenomic <- as.matrix(phenomic)

saveRDS(phenomic, "estimate_accuracy5/phenomic.rds")



# genomic prediction
set.seed(1) 
five_folds <- sample(1:nrow(geno), nrow(geno))
five_folds <- split(five_folds, 1:5)



accuracy_table_five <- data.frame(fold="five folds", 
                                  jth_fold=1:5, 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  p_acc_pear_fix=NA, 
                                  p_acc_pear_mix=NA, 
                                  p_acc_gcor_fix=NA, 
                                  p_acc_gcor_mix=NA, 
                                  h2=NA)
for (i in 1:nrow(accuracy_table_five)){
  print(i)
  tryCatch({
    jth_fold = accuracy_table_five$jth_fold[i]
    selected = five_folds[[jth_fold]][!is.na(phenotype[five_folds[[jth_fold]], "grain_yield"])]
    
    Y <- data.frame(ID=rownames(phenotype)[selected], 
                    obs=phenotype[selected, "grain_yield"])
    kin = kinship[Y$ID, Y$ID]
    rownames(kin) = Y$ID
    colnames(kin) = Y$ID
    
    gmodel = mixed.solve(phenotype[-selected, "grain_yield"], 
                         geno[-selected, ])
    Y$pred = geno[Y$ID, ] %*% gmodel$u
    accuracy_table_five$g_acc_pear[i] = cor(Y$obs, Y$pred)
    
    K = kinship[-selected, -selected]
    gmodel2 = mixed.solve(phenotype[-selected, "grain_yield"], K=K)
    h2 = gmodel2$Vu / (gmodel2$Vu+gmodel2$Ve)
    accuracy_table_five$g_acc_gcor[i] = accuracy_table_five$g_acc_pear[i] / sqrt(h2)
    accuracy_table_five$h2[i] = h2
    
    
    Y2 = Y
    pmodel_mix = mixed.solve(phenotype[-selected, "grain_yield"], 
                               phenomic[-selected, ])
    Y2$pred = phenomic[selected, ] %*% pmodel_mix$u
    accuracy_table_five$p_acc_pear_mix[i] = cor(Y2$obs, Y2$pred)
    accuracy_table_five$p_acc_gcor_mix[i] = 
      estimate_gcor(data=Y2, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
    Y3 = Y
    data = data.frame(grain_yield = phenotype[-selected, ])
    data = cbind(data, phenomic[-selected, ])
    pmodel_fix = lm(grain_yield ~ ., data)
    Y3$pred = phenomic[selected, ] %*% pmodel_fix$coefficients[2:length(pmodel_fix$coefficients)] + 
      pmodel_fix$coefficients[1]
    accuracy_table_five$p_acc_pear_fix[i] = cor(Y3$obs, Y3$pred)
    accuracy_table_five$p_acc_gcor_fix[i] = 
      estimate_gcor(data=Y3, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
  }, error=function(e){cat(i, "\n", "Error:", conditionMessage(e), "\n")})
}

saveRDS(accuracy_table_five, "estimate_accuracy5/accuracy_table_five.rds")



accuracy_table <- accuracy_table_five
accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=6), 1:2]
accuracy_table_long$type <- rep(colnames(accuracy_table)[3:8], nrow(accuracy_table))
accuracy_table_long$accuracy <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*6+1) : ((i-1)*6+6)] = 
    accuracy_table[i, 3:8]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
accuracy_table_long <- na.omit(accuracy_table_long)
accuracy_table_long$type <- factor(accuracy_table_long$type, levels=c("g_acc_pear", "p_acc_pear_mix", 
                                                                      "p_acc_pear_fix", "g_acc_gcor", 
                                                                      "p_acc_gcor_mix", "p_acc_gcor_fix"))
accuracy_table_long$adjustment <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "p_acc_pear"), 
                                         "before adjustment", "after adjustment")
accuracy_table_long$adjustment <- factor(accuracy_table_long$adjustment, 
                                         levels=c("before adjustment", "after adjustment"))
accuracy_table_long$prediction <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "g_acc_gcor"), 
                                         "genomic", "phenomic")

saveRDS(accuracy_table_long, "estimate_accuracy5/accuracy_table_long.rds") 



pdf(paste("estimate_accuracy5/plots/", "accuracy.pdf", sep=""), width=14)
ggplot(accuracy_table_long, 
       aes(x=type, y=accuracy, fill=type)) + 
  geom_boxplot() + 
  facet_wrap(~fold) + 
  scale_fill_hue(labels=c("g_acc_pear"="genomic accuracy Pearson", 
                          "g_acc_gcor"="genomic accuracy adjusted", 
                          "p_acc_pear_mix"="phenomic accuracy Pearson, mixed model", 
                          "p_acc_pear_fix"="phenomic accuracy Pearson, fixed model", 
                          "p_acc_gcor_mix"="phenomic accuracy adjusted, mixed model", 
                          "p_acc_gcor_fix"="phenomic accuracy adjusted, fixed model")) + 
  theme_minimal_grid(font_size=10) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        legend.position="bottom")

dev.off()











