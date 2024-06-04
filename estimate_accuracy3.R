# Charlotte Brault on grapevine
setwd("~/breeder_equation_project/")
source("codes/Estimate_gcor_prediction.R")
library(rrBLUP)
library(ggplot2)
library(cowplot)
# library(vcfR)



# diversity panel 
# genotype data
geno_diver <- read.csv("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/results/geno-p279_formated.tsv", sep="\t")
# names_diver <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/genotypic/CC279_GBS/origin_CC279_UTF8.tsv", sep="\t")

# malnorm_diver <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/phenotypic/CC279/asso-279_DLVitis_malnorm_lme4_geno-values.tsv", sep="\t")

# phenotypic data
phenotypic_diver <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/results/BLUPs-p279_formated-scaled.tsv", sep="\t")

# phenomic data
load("data/interest_of_phenomic_prediction/phenomic_pred/results/spectrum_correction/wood_p279_2020_nirs-list-2-blocks.Rdata")

wood_2020_diver <- nirs.list[["der1"]][["NIRS"]]



# half_diallelic population

# phenotypic_half <- read.delim("data/across_population_genomic_prediction_in_grapevine/phenomic_pred/data/phenotypic/diallel/pheno-diallel_with-outlier-columns.tsv", sep="\t")








































