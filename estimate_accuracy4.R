# Margaret Krause on wheat
setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)



data_path <- "data/hyperspectral_reflectance_derived_relationship_matrices/"

# load pedigree relationship matrix
pedi_rel <- read.csv(paste(data_path, "Krause_et_al_2018_Pedigree.csv", sep=""), row.names = 1, check.names = F)
dim(pedi_rel)
# [1] 3771 3771

# load the phenomic data
# pheno_indi <- read.csv(paste(data_path, "Krause_et_al_2018_Hyper_BLUEs_Individual_Time_Points.csv", 
#                              sep=""), check.names = F)
# dim(pheno_indi)
# # [1] 133639     66
pheno_stages <- read.csv(paste(data_path, "Krause_et_al_2018_Hyper_BLUEs_Growth_Stages_VEG_HEAD_GF_ALL.csv", 
                               sep=""), check.names = F)
dim(pheno_stages)
# [1] 65271    66
phenomic <- pheno_stages[pheno_stages$Growth_Stage=="ALL", ]
dim(phenomic)
# [1] 17822    66
str(phenomic)
phenomic$GID <- as.character(phenomic$GID)
length(unique(phenomic$GID))
# [1] 3771

# load the phenotypic data
phenotypes <- read.csv(paste(data_path, "Krause_et_al_2018_Yield_BLUEs.csv", sep=""), check.names = F)
dim(phenotypes)
# [1] 18855     4
str(phenotypes)
phenotypes$GID <- as.character(phenotypes$GID)
length(unique(phenotypes$GID))
# [1] 3771



# there are more phenotype observations than phenomic observations
# which phenomic observation is not present in phenotype data?
# this is taking too long, don't do it
# for (i in 1:nrow(phenomic)){
#   if (!any( apply(phenotypes[, c("GID", "Managed_Treatment")], 1, 
#                   setequal, phenomic[, c("GID", "Managed_Treatment")]) )){
#     print(i)
#   }
# }



# treatments 
treatments <- unique(phenotypes$Managed_Treatment)



geno_model




