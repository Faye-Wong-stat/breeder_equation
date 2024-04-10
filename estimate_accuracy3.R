# Charlotte Brault on grapevine
setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)
# library(vcfR)



# diversity panel 
geno_diver <- read.csv("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/genotypic/CC279_GBS/279_SNPs_MAF5.csv", sep="\t")
names_diver <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/genotypic/CC279_GBS/origin_CC279_UTF8.tsv", sep="\t")

malnorm_diver <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/phenotypic/CC279/asso-279_DLVitis_malnorm_lme4_geno-values.tsv", sep="\t")



# half_diallelic population

phenotypic_half <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/phenotypic/diallel/pheno-diallel_with-outlier-columns.tsv", sep="\t")
