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
dim(geno)
# [1]  3442 18013

saveRDS(geno, "estimate_accuracy5/geno.rds")
geno <- readRDS("estimate_accuracy5/geno.rds")

kinship <- A.mat(geno-1)
saveRDS(kinship, "estimate_accuracy5/kinship.rds")
kinship <- readRDS("estimate_accuracy5/kinship.rds")



# organize phenotyic/phenomic data
dim(pheno)
# [1] 4347   44
length(unique(pheno$fullsamplename))
# [1] 3442
model <- lmer(grain_yield ~ fullsamplename + (1|environment), pheno)
emm_summ <- summary(emmeans(model, "fullsamplename"))

phenotype <- data.frame(grain_yield = emm_summ[, 2][match(rownames(geno), emm_summ[, 1])])
rownames(phenotype) <- emm_summ[, 1][match(rownames(geno), emm_summ[, 1])]
phenotype <- as.matrix(phenotype)

saveRDS(phenotype, "estimate_accuracy5/phenotype.rds")
phenotype <- readRDS("estimate_accuracy5/phenotype.rds")

models <- list()
for (i in colnames(pheno[, 21:44])){
  models[[i]] = lmer(get(i) ~ fullsamplename + (1|environment), pheno)
}
emm_summs <- list()
for (i in 1:length(models)){
  emm_summs[[i]] = summary(emmeans(models[[i]], "fullsamplename"))
}

saveRDS(emm_summs, "estimate_accuracy5/emm_summs.rds")
emm_summs <- readRDS("estimate_accuracy5/emm_summs.rds")

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
phenomic <- readRDS("estimate_accuracy5/phenomic.rds")



# genomic prediction
set.seed(1) 
# five_folds <- vector(mode="list", length=30)
# five_folds <- lapply(five_folds, FUN=function(x){
#   sample(1:nrow(geno), round(nrow(geno)/5))
# })
five_folds <- sample(1:nrow(geno), nrow(geno)) ### 3, 4, 5 folds
five_folds <- split(five_folds, 1:5)



accuracy_table_five <- data.frame(fold="five folds", 
                                  jth_fold=1:5, 
                                  g_acc_pear=NA, 
                                  g_acc_gcor=NA, 
                                  # p_acc_pear=NA, 
                                  # p_acc_gcor=NA, 
                                  p_acc_pear_fix=NA,
                                  p_acc_pear_mix=NA,
                                  p_acc_gcor_fix=NA,
                                  p_acc_gcor_mix=NA,
                                  h2=NA)

for (j in 1:5){
  tryCatch({
    selected = five_folds[[j]][!is.na(phenotype[five_folds[[j]], "grain_yield"])]
    
    Y = data.frame(ID=rownames(phenotype)[selected], 
                    obs=phenotype[selected, "grain_yield"])
    
    kin = kinship[Y$ID, Y$ID]
    scale_by = mean(diag(kin))
    K = kinship/scale_by
    
    K_train = K[rownames(geno[-selected, ]), rownames(geno[-selected, ])]
    K_test = K[Y$ID, Y$ID]
    
    gmodel = mixed.solve(phenotype[-selected, "grain_yield"], 
                         geno[-selected, ])
    Y$pred = geno[Y$ID, ] %*% gmodel$u
    
    accuracy_table_five$g_acc_pear[j] = cor(Y$obs, Y$pred)
    
    gmodel2 = mixed.solve(phenotype[selected, "grain_yield"], K=K_test)
    h2 = gmodel2$Vu / var(Y$obs)
    accuracy_table_five$g_acc_gcor[j] = cor(Y$obs, Y$pred) / sqrt(h2)
    accuracy_table_five$h2[j] = h2
    
    Y2 = Y
    pmodel = mixed.solve(phenotype[-selected, "grain_yield"], 
                         phenomic[-selected, ])
    Y2$pred = phenomic[selected, ] %*% pmodel$u
    
    accuracy_table_five$p_acc_pear_mix[j] = cor(Y2$obs, Y2$pred)
    accuracy_table_five$p_acc_gcor_mix[j] = 
      estimate_gcor(data=Y2, Knn=K_test, method="MCMCglmm", normalize=F)[1]
    
    Y3 = Y
    data = data.frame(grain_yield = phenotype[-selected, ])
    data = cbind(data, phenomic[-selected, ])
    pmodel_fix = lm(grain_yield ~ ., data)
    Y3$pred = phenomic[selected, ] %*% pmodel_fix$coefficients[2:length(pmodel_fix$coefficients)] +
      pmodel_fix$coefficients[1]
    accuracy_table_five$p_acc_pear_fix[j] = cor(Y3$obs, Y3$pred)
    accuracy_table_five$p_acc_gcor_fix[j] =
      estimate_gcor(data=Y3, Knn=kin, method="MCMCglmm", normalize=F)[1]
    
  }, error=function(e){cat(j, "\n", "Error:", conditionMessage(e), "\n")})
}

# for (i in 1:nrow(accuracy_table_five)){
#   print(i)
#   tryCatch({
#     jth_fold = accuracy_table_five$jth_fold[i]
#     selected = five_folds[[jth_fold]][!is.na(phenotype[five_folds[[jth_fold]], "grain_yield"])]
#     
#     Y <- data.frame(ID=rownames(phenotype)[selected], 
#                     obs=phenotype[selected, "grain_yield"])
#     kin = kinship[Y$ID, Y$ID]
#     rownames(kin) = Y$ID
#     colnames(kin) = Y$ID
#     
#     gmodel = mixed.solve(phenotype[-selected, "grain_yield"], 
#                          geno[-selected, ])
#     Y$pred = geno[Y$ID, ] %*% gmodel$u
#     accuracy_table_five$g_acc_pear[i] = cor(Y$obs, Y$pred)
#     
#     K = kinship[-selected, -selected]
#     gmodel2 = mixed.solve(phenotype[-selected, "grain_yield"], K=K)
#     h2 = gmodel2$Vu / (gmodel2$Vu+gmodel2$Ve)
#     accuracy_table_five$g_acc_gcor[i] = accuracy_table_five$g_acc_pear[i] / sqrt(h2)
#     accuracy_table_five$h2[i] = h2
#     
#     
#     Y2 = Y
#     pmodel_mix = mixed.solve(phenotype[-selected, "grain_yield"], 
#                                phenomic[-selected, ])
#     Y2$pred = phenomic[selected, ] %*% pmodel_mix$u
#     accuracy_table_five$p_acc_pear_mix[i] = cor(Y2$obs, Y2$pred)
#     accuracy_table_five$p_acc_gcor_mix[i] = 
#       estimate_gcor(data=Y2, Knn=kin, method="MCMCglmm", normalize=F)[1]
#     
#     Y3 = Y
#     data = data.frame(grain_yield = phenotype[-selected, ])
#     data = cbind(data, phenomic[-selected, ])
#     pmodel_fix = lm(grain_yield ~ ., data)
#     Y3$pred = phenomic[selected, ] %*% pmodel_fix$coefficients[2:length(pmodel_fix$coefficients)] + 
#       pmodel_fix$coefficients[1]
#     accuracy_table_five$p_acc_pear_fix[i] = cor(Y3$obs, Y3$pred)
#     accuracy_table_five$p_acc_gcor_fix[i] = 
#       estimate_gcor(data=Y3, Knn=kin, method="MCMCglmm", normalize=F)[1]
#     
#   }, error=function(e){cat(i, "\n", "Error:", conditionMessage(e), "\n")})
# }

saveRDS(accuracy_table_five, "estimate_accuracy5/accuracy_table_five.rds")
accuracy_table_five <- readRDS("estimate_accuracy5/accuracy_table_five.rds")



accuracy_table <- accuracy_table_five
# accuracy_table_long <- accuracy_table[rep(1:nrow(accuracy_table), each=6), 2]
# accuracy_table_long$type <- rep(colnames(accuracy_table)[3:8], nrow(accuracy_table))
accuracy_table_long <- data.frame(jth_fold=accuracy_table[rep(1:nrow(accuracy_table), each=6), 2], 
                                  type=rep(colnames(accuracy_table)[3:8], nrow(accuracy_table)))
accuracy_table_long$accuracy <- NA
accuracy_table_long$h2 <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*6+1) : ((i-1)*6+6)] = 
    accuracy_table[i, 3:8]
  accuracy_table_long$h2[((i-1)*6+1) : ((i-1)*6+6)] = 
    accuracy_table[i, 9]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)
# accuracy_table_long <- na.omit(accuracy_table_long)
# accuracy_table_long$h2 <- mean(accuracy_table$h2)
accuracy_table_long$type <- factor(accuracy_table_long$type, levels=c("g_acc_pear", "p_acc_pear_mix", 
                                                                      "p_acc_pear_fix", "g_acc_gcor", 
                                                                      "p_acc_gcor_mix", "p_acc_gcor_fix"))
accuracy_table_long$adjustment <- ifelse(accuracy_table_long$type %in% 
                                           c("g_acc_pear", "p_acc_pear_mix", "p_acc_pear_fix"), 
                                         "predictive ability", "prediction accuracy")
accuracy_table_long$adjustment <- factor(accuracy_table_long$adjustment, 
                                         levels=c("predictive ability", "prediction accuracy"))
accuracy_table_long$`prediction method` <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "g_acc_gcor"), 
                                         "genomic", ifelse(accuracy_table_long$type %in% 
                                                             c("p_acc_pear_mix", "p_acc_gcor_mix"), 
                                                           "phenomic-mixed", "phenomic-fixed"))
accuracy_table_long$`prediction method` <- factor(accuracy_table_long$`prediction method`, 
                                                  levels=c("genomic", "phenomic-mixed", "phenomic-fixed"))

saveRDS(accuracy_table_long, "estimate_accuracy5/accuracy_table_long.rds") 

p1 = ggplot(accuracy_table_long[accuracy_table_long$adjustment == "before correction", ],
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
  ggtitle(paste("trait: grain yield", ", h2: ", accuracy_table_long$h2[1], sep=""))

p2 = ggplot(accuracy_table_long[accuracy_table_long$adjustment == "after correction", ],
            aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) +
  geom_boxplot() +
  facet_wrap(~adjustment) +
  ylim(0, 1) + 
  ylab("prediction accuracy") + 
  theme_minimal_grid(font_size=12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none") 
# +
  # ggtitle(paste("trait: grain yield", ", h2: ", h2[i], sep=""))

save_plot(paste("estimate_accuracy5/plots/", "accuracy_fixed.pdf", sep=""), 
          plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))

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











