# Margaret Krause on wheat
setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(BGLR)
library(ggplot2)
library(cowplot)



data_path <- "data/hyperspectral_reflectance_derived_relationship_matrices/"


data_DTHD <- read.csv(paste(data_path, "Krause_et_al_2018_Yield_iid_BLUPs_DTHD.csv", sep=""), 
                      check.names = F)


# load pedigree relationship matrix
pedi_rel <- read.csv(paste(data_path, "Krause_et_al_2018_Pedigree.csv", sep=""), row.names = 1, 
                     check.names = F)
dim(pedi_rel)
# [1] 3771 3771
pedi_rel <- as.matrix(pedi_rel)



# load the phenomic data
pheno_stages <- read.csv(paste(data_path, "Krause_et_al_2018_Hyper_BLUEs_Growth_Stages_VEG_HEAD_GF_ALL.csv", 
                               sep=""), check.names = F)
dim(pheno_stages)
# [1] 65271    66
pheno_stages[1:5, 1:5]
#       GID Breeding_Cycle Managed_Treatment Growth_Stage Hyper_BLUE_398nm
# 1 6569128        2013-14  Moderate Drought          VEG         91.93150
# 2 6569128        2013-14  Moderate Drought         HEAD         70.02356
# 3 6569128        2013-14  Moderate Drought           GF         61.04231
# 4 6569128        2013-14  Moderate Drought          ALL         76.76667
# 5 6569128        2013-14       Optimal Bed          VEG         84.88888

phenomic <- pheno_stages[pheno_stages$Growth_Stage=="ALL", ]
dim(phenomic)
# [1] 17822    66
# str(phenomic)
phenomic$GID <- as.character(phenomic$GID)
length(unique(phenomic$GID))
# [1] 3771

phenomic$site_year <- paste(phenomic$Breeding_Cycle, phenomic$Managed_Treatment, sep="_")
phenomic <- phenomic[, c(1, 67, 5:66)]
dim(phenomic)
# [1] 17822    64
phenomic[1:5, 1:5]
#        GID                site_year Hyper_BLUE_398nm Hyper_BLUE_405nm
# 4  6569128 2013-14_Moderate Drought         76.76667        104.09324
# 8  6569128      2013-14_Optimal Bed         63.84228         85.44301
# 11 6569128             2013-14_Heat         93.74714        120.91550
# 15 6569128   2013-14_Severe Drought         68.40447         98.56469
# 19 6569128     2013-14_Optimal Flat         65.08903         85.74871
#    Hyper_BLUE_412nm
# 4          121.3816
# 8          106.0869
# 11         150.8088
# 15         110.9700
# 19          93.1260
length(unique(interaction(phenomic$GID, phenomic$site_year)))
# [1] 17822
phenomic[, 3:64] <- scale(phenomic[, 3:64])

phenomic_rel <- as.matrix(phenomic[, 3:64]) %*% t(as.matrix(phenomic[, 3:64]))
phenomic_rel <- phenomic_rel/62

saveRDS(phenomic_rel, "phenomic_relationship.rds")
phenomic_rel <- readRDS("phenomic_relationship.rds")



# load the phenotypic data
phenotypes <- read.csv(paste(data_path, "Krause_et_al_2018_Yield_BLUEs.csv", sep=""), check.names = F)
dim(phenotypes)
# [1] 18855     4
str(phenotypes)
# 'data.frame':   18855 obs. of  4 variables:
# $ GID              : int  6569128 6569128 6569128 6569128 6569128 6688880 6688880 6688880 6688880 6688880 ...
# $ Breeding Cycle   : chr  "2013-14" "2013-14" "2013-14" "2013-14" ...
# $ Managed_Treatment: chr  "Moderate Drought" "Optimal Bed" "Heat" "Severe Drought" ...
# $ Grain_Yield_BLUE : num  3.43 6.07 1.54 2.02 6.12 ...
phenotypes$GID <- as.character(phenotypes$GID)
length(unique(phenotypes$GID))
# [1] 3771

phenotypes$site_year <- paste(phenotypes$`Breeding Cycle`, phenotypes$Managed_Treatment, sep="_")
phenotypes <- phenotypes[, c(1, 5, 4)]
dim(phenotypes)
# [1] 18855     3
length(unique(interaction(phenotypes$GID, phenotypes$site_year)))
# [1] 18855

sum(! unique(interaction(phenomic$GID, phenomic$site_year)) %in% 
      unique(interaction(phenotypes$GID, phenotypes$site_year)))
# [1] 0



# there are more phenotype observations than phenomic observations
# which phenomic observation is not present in phenotype data?

# select only the phenotypes with combinations shown in phenomic data
phenomic$ID <- paste(phenomic$GID, phenomic$site_year, sep=".")
phenotypes$ID <- paste(phenotypes$GID, phenotypes$site_year, sep=".")
# which phenomic data missing? 
missing_phenomic <- phenotypes[! phenotypes$ID %in% phenomic$ID, ]
dim(missing_phenomic)
# [1] 1033    4
phenotypes <- phenotypes[phenotypes$ID %in% phenomic$ID, ]
phenotypes <- phenotypes[match(phenotypes$ID, phenomic$ID), ]
phenotypes <- phenotypes[, ! colnames(phenotypes) %in% "ID"]
phenomic <- phenomic[, ! colnames(phenomic) %in% "ID"]



# # make the fixed effect design matrix
# site_year <- unique(phenotypes$site_year)
# fixed_effect <- matrix(NA, nrow=nrow(phenotypes), ncol=18)
# # fixed_effect <- as.data.frame(fixed_effect)
# for (i in 1:nrow(phenotypes)) {
#   fixed_effect[i, ] = as.numeric(phenotypes$site_year[i] == site_year[2:19])
# }
# dim(fixed_effect)
# # [1] 17822    18



# set.seed(1)
# five_folds <- sample(1:dim(phenotypes)[1], dim(phenotypes)[1]) ### 3, 4, 5 folds
# five_folds <- split(five_folds, 1:5)
# four_folds <- sample(1:dim(phenotypes)[1], dim(phenotypes)[1])
# four_folds <- split(four_folds, 1:4)
# three_folds <- sample(1:dim(phenotypes)[1], dim(phenotypes)[1])
# three_folds <- split(three_folds, 1:3)

site_year <- unique(phenotypes$site_year)
# No_geno <- rep(NA, 19)
# for (i in 1:19){
#   No_geno[i] = length(unique(phenotypes[phenotypes$site_year == site_year[i], "GID"]))
# }
accuracy_table <- data.frame(site_year = site_year, 
                             g_acc_pear=NA, 
                             g_acc_gcor=NA, 
                             p_acc_pear=NA, 
                             p_acc_gcor=NA, 
                             h2 = NA)

set.seed(1)
for (i in 1:19){
  tryCatch({
    testing_data = phenotypes[phenotypes$site_year == site_year[i], ]
  training_data = phenotypes[phenotypes$site_year != site_year[i], ]
  
  chol_K = t(chol(pedi_rel))
  ZK = chol_K[training_data$GID, ]
  
  gmodel = BGLR(y = training_data$Grain_Yield_BLUE, 
                ETA = list(list(~factor(site_year), data=training_data, model="FIXED"), 
                           list(X=ZK, model="BRR")), 
                verbose=F, 
                saveAt="estimate_accuracy4/")
  
  Y = data.frame(ID = testing_data[, "GID"], 
                 obs = testing_data[, "Grain_Yield_BLUE"])
  Y$pred = gmodel$ETA[[2]]$b[Y$ID]
  accuracy_table[i, "g_acc_pear"] = cor(Y$pred, Y$obs)
  
  h2 = gmodel$ETA[[2]]$varB / var(Y$obs)
  accuracy_table[i, "g_acc_gcor"] = accuracy_table[i, "g_acc_pear"] / sqrt(h2)
  accuracy_table[i, "h2"] = h2
  
  
  
  pmodel = BGLR(y = training_data$Grain_Yield_BLUE, 
                ETA = list(list(~factor(site_year), data=training_data, model="FIXED"), 
                           list(X=phenomic[phenotypes$site_year != site_year[i], 3:64], model="BRR")), 
                verbose=F, 
                saveAt="estimate_accuracy4/")
  
  Y2 = Y
  Y2$pred = as.matrix(phenomic[phenotypes$site_year == site_year[i], 3:64]) %*% pmodel$ETA[[2]]$b
  
  accuracy_table[i, "p_acc_pear"] = cor(Y2$pred, Y2$obs)
  
  kin = pedi_rel[testing_data$GID, testing_data$GID]
  rownames(kin) = Y2$ID
  colnames(kin) = Y2$ID
  
  accuracy_table[i, "p_acc_gcor"] = estimate_gcor(data=Y2, Knn=kin, method="MCMCglmm", normalize=F)[1]
  }, error=function(e){cat(i, "\n", "Error:", conditionMessage(e), "\n")})
}

saveRDS(accuracy_table, "estimate_accuracy4/accuracy_table.rds")
accuracy_table <- readRDS("estimate_accuracy4/accuracy_table.rds")



accuracy_table_long <- data.frame(site_year = accuracy_table[rep(1:nrow(accuracy_table), each=4), 1], 
                                  type = rep(colnames(accuracy_table)[2:5], nrow(accuracy_table)))
accuracy_table_long$accuracy <- NA
accuracy_table_long$h2 <- NA
for (i in 1:nrow(accuracy_table)){
  accuracy_table_long$accuracy[((i-1)*4+1) : ((i-1)*4+4)] = 
    accuracy_table[i, 2:5]
  accuracy_table_long$h2[((i-1)*4+1) : ((i-1)*4+4)] = 
    accuracy_table[i, 6]
}
accuracy_table_long$accuracy <- unlist(accuracy_table_long$accuracy)

# accuracy_table_long$h2 <- rep(accuracy_table$h2, each=4)
# accuracy_table_long$h2 <- as.character(round(accuracy_table_long$h2, 5))
accuracy_table_long$type <- factor(accuracy_table_long$type, levels=c("g_acc_pear", "p_acc_pear", 
                                                                      "g_acc_gcor", "p_acc_gcor"))
accuracy_table_long$adjustment <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "p_acc_pear"), 
                                         "predictive ability", "prediction accuracy")
accuracy_table_long$adjustment <- factor(accuracy_table_long$adjustment, 
                                         levels=c("predictive ability", "prediction accuracy"))
accuracy_table_long$`prediction method` <- ifelse(accuracy_table_long$type %in% c("g_acc_pear", "g_acc_gcor"), 
                                                  "genomic", "phenomic")

saveRDS(accuracy_table_long, "estimate_accuracy4/accuracy_table_long.rds")
accuracy_table_long <- readRDS("estimate_accurac4/accuracy_table_long.rds")



h2 <- mean(accuracy_table$h2)

p1 <- ggplot(accuracy_table_long[accuracy_table_long$adjustment == "before correction", ], 
             aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) + 
  geom_boxplot() +
  facet_wrap(~adjustment) +
  ylim(-0.6, 1) +
  ylab("predictive ability") + 
  theme_minimal_grid(font_size=12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="bottom") + 
  ggtitle(paste("h2: ", h2, sep=""))
p2 <- ggplot(accuracy_table_long[accuracy_table_long$adjustment == "after correction", ],
                  aes(x=`prediction method`, y=accuracy, fill=`prediction method`)) +
  geom_boxplot() +
  facet_wrap(~adjustment) +
  ylim(-0.6, 1) +
  ylab("prediction accuracy") + 
  theme_minimal_grid(font_size=12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none") 

save_plot(paste("estimate_accuracy4/plots/", "accuracy.pdf", sep=""), 
          plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))

for (i in 1:19){
  p1 = ggplot(accuracy_table_long[accuracy_table_long$site_year == site_year[i] & 
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
    ggtitle(paste("treatment: ", site_year[i], ", h2: ", h2[i], sep=""))
  
  p2 = ggplot(accuracy_table_long[accuracy_table_long$trait == site_year[i] & 
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
          legend.position="none") 
  
  save_plot(paste("estimate_accuracy4/plots/", site_year[i], "_accuracy.pdf", sep=""), 
            plot_grid(p1, p2, nrow=1, align = "vh", axis="btlr"))
}


