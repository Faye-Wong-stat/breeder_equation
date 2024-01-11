setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)
library(vcfR)



path <- "data/interest_of_phenomic_prediction_as_an_alternative_to_genomic_predition/data/phenomic_pred/"
path2 <- "data/interest_of_phenomic_prediction_as_an_alternative_to_genomic_predition/script/phenomic_pred/"
geno_diversity <- read.table(paste(path, "data/genotypic/CC279_GBS/279_SNPs_MAF5.csv", sep=""))
# info_diversity <- read.delim(paste(path, "genotypic/CC279_GBS/origin_CC279_UTF8.tsv", sep=""), sep="\t")

# file_names_geno_diallel <- list.files(paste(path, "genotypic/diallel_GBS", sep=""), ".vcf")
# geno_diallel_vcf <- list()
# for (i in 1:length(file_names_geno_diallel)){
#   geno_diallel_vcf[[i]] = read.vcfR(paste(path, "genotypic/diallel_GBS/", file_names_geno_diallel[i], 
#                                           sep=""))
# }
# 
# A <- as.data.frame(geno_diallel_vcf[[1]]@fix)
# B <- as.data.frame(geno_diallel_vcf[[1]]@gt)



