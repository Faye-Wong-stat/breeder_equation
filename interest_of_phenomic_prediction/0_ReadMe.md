---
title: ReadMe
author: Brault Charlotte
bibliography: /home/brault/Documents/Articles/Article2/dataverse/biblio_zotero.bib
link-citations: true
date: "" 
---


This document describes the input files, scripts and output files generated for the publication "Across-population genomic prediction in grapevine opens up promising prospects for breeding". This document is in Markdown format.

A global overview of the study and associated scripts (mind-map) is available [here](https://framindmap.org/c/maps/1099842/public)

Table in **.tsv** format were automatically converted into **.tab** by data.inrae.fr.


# Input data 

## Marker data

### Diversity panel population

* `279_SNPs_MAF5.csv` This table is the SNP genotypic matrix for the diversity panel population composed of 277 *Vitis vinifera* cultivars and 83,264 SNP markers. These data were derived from [@flutreGenomewideAssociationPrediction2020] analysis, with a MAF filter at 5% and without LD pruning.

* `origin_CC279_UTF8.tsv` This table is the correspondence between cultivar number in `279_SNPs_MAF5.csv` and cultivar name. It also contains information about subpopulations (WW, WE, TE) composition, based on both the original SSR classification by [@nicolasGeneticDiversityLinkage2016] and the SNP classification by [@flutreGenomewideAssociationPrediction2020].

### Half-diallel population

SNP genotypic data for the half-diallel population were previously derived by [@telloNovelHighdensityGrapevine2019]; we used the vcf files (one per cross) from the associated [dataverse](https://doi.org/10.15454/99SLGP).

* `Variant_callingPOP31.filtered.gq40.env.sansdoublon3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Syrah x Pinot Noir cross, with 60 genotypes and 173,609 markers (after filtering).

* `Variant_callingPOP32.filtered.gq40.env.noreplicate2.chr0.vcf` This is the file with genotypic data in vcf format for Cabernet-Sauvignon x Pinot Noir cross, with 68 genotypes and 157,131 markers (after filtering).

* `Variant_callingPOP33.filtered.gq40.env.noreplicate2.chr0.vcf` This is the file with genotypic data in vcf format for Grenache x Pinot Noir cross, with 65 genotypes and 157,875 markers (after filtering).

* `Variant_callingPOP34.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Terret Noir x Pinot Noir cross, with 62 genotypes and 207,203 markers (after filtering).

* `Variant_callingPOP35.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Cabernet-Sauvignon x Syrah cross, with 66 genotypes and 250,749 markers (after filtering).

* `Variant_callingPOP36.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Grenache x Syrah cross, with 61 genotypes and 174,765 markers (after filtering).

* `Variant_callingPOP37.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Terret Noir x Syrah cross, with 65 genotypes and 184,129 markers (after filtering).

* `Variant_callingPOP38.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Grenache x Cabernet-Sauvignon cross, with 65 genotypes and 112,751 markers (after filtering).

* `Variant_callingPOP39.filtered.gq40.env.noreplicate3.chr0.recode.vcf` This is the file with genotypic data in vcf format for Terret Noir x Cabernet-Sauvignon cross, with 63 genotypes and 180,423 markers (after filtering).

* `Variant_callingPOP40.filtered.gq40.env.noreplicate2.chr0.vcf` This is the file with genotypic data in vcf format for Terret Noir x Grenache cross, with 69 genotypes and 153,931 markers (after filtering).

* `corresp-genotype-names_labo-code-intro.tsv` This table is the correspondence between half-diallel genotype name in vcf files and genotype names used in further analysis.


## Phenotypic data and BLUPs of genotypic values 

* For the diversity panel, BLUPs of genotypic values for the 15 studied traits were already available [@flutreGenomewideAssociationPrediction2020].


	* `asso-279_DLVitis_malnorm_lme4_geno-values.tsv` for malic acid at ripe stage trait.
	* `asso-279_DLVitis_tarnorm_lme4_geno-values.tsv` for tartaric acid at ripe stage trait.
 	* `asso-279_DLVitis_shinorm_lme4_geno-values.tsv` for shikimic acid at ripe stage trait.
 	* `asso-279_DLVitis_shiktar_lme4_geno-values.tsv` for shikimic/tartaric acid at ripe stage trait.
 	* `asso-279_DLVitis_maltar_lme4_geno-values.tsv` for malic/tartaric acid at ripe stage trait.
 	* `asso-279_DLVitis_verday_lme4_geno-values.tsv` for veraison date (onset of ripening) trait.
 	* `asso-279_DLVitis_samplday_lme4_geno-values.tsv` for maturity date trait.
 	* `asso-279_DLVitis_vermatu_lme4_geno-values.tsv` for interval between veraison and maturity trait.
 	* `asso-279_DLVitis_clucomp_lme4_geno-values.tsv` for cluster compactness trait.
 	* `asso-279_DLVitis_nbclu_lme4_geno-values.tsv` for number of clusters trait.
 	* `asso-279_DLVitis_mcl_lme4_geno-values.tsv` for mean cluster length trait.
 	* `asso-279_DLVitis_mcwi_lme4_geno-values.tsv` for mean cluster width trait.
 	* `asso-279_DLVitis_mcw_lme4_geno-values.tsv` for mean cluster weight trait.
 	* `asso-279_DLVitis_mbw_lme4_geno-values.tsv` for mean berry weight trait.
 	* `asso-279_DLVitis_vigour_lme4_geno-values.tsv` for vigour trait.
 	
 	
 	
 	
* For the half-diallel, the raw phenotypic data for the 15 studied traits is in the table `pheno-diallel_with-outlier-columns.tsv`. We visually checked for outlier values and we added a logical column *trait.out*, this indicating lines that should be removed from the analysis. Spatial (field layout) and temporal (phenotyping year) indications are also in this file.


# Scripts

## Source scripts

* `get_colours_diallel.R` is a small R-script with a function defining a specific color for each half-diallel cross, panel subpopulation and half-diallel parent.

* `half-diallel_analysis_functions.R` is a R-script gathering useful functions for half-diallel analysis.


## Imputation of half-diallel marker data

### Impute marker data in half-diallel population

* `format_diallel_markers.Rmd` This script takes as input genotypic data in vcf format (see above) and generates `genomic_data_gene-dose_diallel-all-pops_for_runFimpute_missdatmax-marker-0.8.Rdata`, the genomic data for the 10 crosses with all available markers. Markers were checked for segregation distortion and markers with more than 80% missing data were removed.

* `pedigree_diallel_parents_FImpute.tsv` is the pedigree file for the half-diallel population used in FImpute3.

* `diallel_markers_fimpute.Rmd` This script takes as inputs `pedigree_diallel_parents_FImpute.tsv` and `genomic_data_gene-dose_diallel-all-pops_for_runFimpute_missdatmax-marker-0.8.Rdata` and impute missing marker data with FImpute3 software version 2.2, [@sargolzaeiNewApproachEfficient2014]. It generates the file `genotypes_diallel_ped_86k_fimpute_missdatmax-marker-0.8.tsv`.


## Phenotypic data analysis


### Half-diallel

* First, raw phenotypic data were combined into a single file in a consistent format and raw phenotypic data were plotted to detect outlier values. All raw phenotypic data are combined in `pheno-diallel_with-outlier-columns.tsv`.  


* Second, `diallel_with-outlier-columns.tsv` was supplied as input to a separate script for each of the 15 studied traits to perform mixed model selection (after raw data transformation when necessary), to check underluying assumptions, and to estimate variance components, heritability values and BLUPs of genotypic values. For example for **mcl** trait, the script is named `pheno_lme4_diallel_mcl.Rmd`. Each script generated two tables: one with mixed model fitting information (`mcl_lme4_fit-info.tsv`) and one with the BLUPs of genotypic values (`mcl_geno-cross_lme4-values.tsv`).

* Then, tables of mixed model fitting information and BLUPs of genotypic values for each trait and for both populations were compiled in `compile_BLUPs-info_diallel-p279.Rmd`. Fitting information from the diversity panel were retrieved from Table S2 from [@flutreGenomewideAssociationPrediction2020]. This script generates one table with all BLUPs of genotypic values (`diallel_blup_geno-wide_chapitre.tsv`), as well as fit information `diallel_fit-info_blups_chapitre.tsv` (Table S1 and S5).

* `estim_H2_per-cross.Rmd` This script takes as input `pheno-diallel_with-outlier-columns.tsv` and model fitting information `diallel_fit-info_blups_chapitre.tsv`, estimates heritability per cross (`diallel_H2-per-cross_all-traits.tsv`) and generates Figure S2.


### Diversity panel

* `compile_explo-BLUPs_p279.Rmd` compiles the BLUPs of genotypic values derived by [@flutreGenomewideAssociationPrediction2020] for the 15 traits studied (listed in « Phenotypic data »). It generates `blups_p279_chapitre_all-traits_Flutre2020.tsv`.

## Combining half-diallel and p279 population information

* `prepare_files_inter-pop_pred_diall_p279.Rmd`. This script combines, within each population, marker data and estimated genotypic values (BLUPs) (e.g., `diallel_blup_geno-wide_chapitre.tsv` and `blups_p279_chapitre_all-traits_Flutre2020.tsv` tables), keeping only common individuals and scaling BLUPs. It also keeps only markers that are common between the two populations. It generates two files per population:

	* `geno-diallel_formated.tsv` for marker data (32,894 SNPs) in the half-diallel, with 627 genotypes, comprising the five half-diallel parents and offspring of the ten crosses.
	
	* `geno-p279_formated.tsv` for marker data (32,894 SNPS) in the diversity panel, with 277 genotypes. 

	* `BLUPs-diallel_formated-scaled.tsv` for BLUPs of genotypic values (15 traits) in the half-diallel, with 627 genotypes, comprising the five half-diallel parents and offspring of the ten crosses.

	* `BLUPs-p279_formated-scaled.tsv` for BLUPs of genotypic values (15 traits) in the diversity panel, with 277 genotypes. 



## Genomic data analysis

* `study_p279_diallel_relationship.Rmd`. This script analyzes genomic data from half-diallel and diversity panel populations. It takes as inputs `geno-diallel_formated.tsv`, `geno-p279_formated.tsv` and `origin_CC279_UTF8.tsv`. It calculates additive relationships between the two populations and genomic PCA of the diversity panel. It generates tables with additive relationship between half-diallel parents, between diversity panel and half-diallel crosses; and distances on the PCA along the first or first two axes. It plots Figures 1a and S1.

* `p279-diallel_explo-geno_SNPs_karyogram.Rmd`. This script analyzes the set of SNP markers of this study. It takes as input a vcf file for diversity panel population, generated in Flutre et al. 2020, and `geno-p279_formated.tsv` for the list of SNP markers. It generates a plot of marker density along the genome, which is the Figure S15.


## Analysis of phenotypic data for both populations

* `distribution_raw-pheno_diallel-p279_pops.Rmd`. This scripts aims at visualizing distribution of raw phenotypic data for both population (Figure S14).

* `compare_BLUPs_diallel-p279_pops.Rmd`. This script aims at study and highlight differences between both populations, by comparing broad-sense heritability (Figure 1b), distribution of genotypic values (Figure 1c-e, Figure S3).It also plots the mean offspring genotypic value vs parental average predicted genotypic value (Figure S12) and a PCA with half-diallel genotypic vallues for all traits for both half-diallel and diversity panel populations (Figure S4).



## Fit genomic prediction

### Scenario 1a: within half-diallel

* `genom-pred_intra-pop_S1a_RR.Rmd`. This script takes as input `geno-diallel_formated.tsv` and `BLUPs_diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 1a and RR method for all traits within the half-diallel population. For each cross-validation replicate and each trait, it saves predicted individual genotypic values (`GP_S1a_RR_pred-obs.Rdata`), predictive ability within each half-diallel cross (`GP_S1a_RR_per-cross-summary.tsv`) and estimated marker effects (`GP_S1a_RR_mrk-effects.tsv`).


* `genom-pred_intra-pop_S1a_LASSO.Rmd`. This script takes as input `geno-diallel_formated.tsv` and `BLUPs_diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 1a and LASSO method for all traits within the half-diallel population. For each cross-validation replicate and each trait, it saves predicted individual genotypic values (`GP_S1a_LASSO_pred-obs.Rdata`), predictive ability within each half-diallel cross (`GP_S1a_LASSO_per-cross-summary.tsv`) and estimated marker effects (`GP_S1a_LASSO_mrk-effects.tsv`).


### Scenario 1b: between half-sib families

* `genom-pred_intra-pop_S1b_RR.Rmd`. This script takes as input `geno-diallel_formated.tsv` and `BLUPs_diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 1a and RR method for all traits acros half-sib families of the half-diallel. For each cross-validation replicate and each trait, it saves predicted individual genotypic values (`GP_S1b_RR_pred-obs.Rdata`), predictive ability within each half-diallel cross (`GP_S1b_RR_per-cross-summary.tsv`) and estimated marker effects (`GP_S1b_RR_mrk-effects.tsv`).

* `genom-pred_intra-pop_S1b_LASSO.Rmd`. This script takes as input `geno-diallel_formated.tsv` and `BLUPs_diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 1a and LASSO method for all traits acros half-sib families of the half-diallel. For each cross-validation replicate and each trait, it saves predicted individual genotypic values (`GP_S1b_LASSO_pred-obs.Rdata`), predictive ability within each half-diallel cross (`GP_S1b_LASSO_per-cross-summary.tsv`) and estimated marker effects (`GP_S1b_LASSO_mrk-effects.tsv`).


### Scenario 2: across-population 

*  `genom-pred_inter-pop_S2_RR.Rmd`. This script takes as inputs `geno-p279_formated.tsv`, `BLUPs-p279_formated-scaled.tsv`, `geno-diallel_formated.tsv` and `BLUPs-diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 2 with RR method. For each trait, it saves predicted individual genotypic values (`GP_S2_RR_pred.tsv`), predictive ability within each half-diallel cross (`GP_S2_RR_per-cross-summary.tsv`) and estimated marker effects (`GP_S2_RR_mrk-effects.tsv`).

* `genom-pred_inter-pop_S2_LASSO.Rmd`. This script takes as inputs `geno-p279_formated.tsv`, `BLUPs-p279_formated-scaled.tsv`, `geno-diallel_formated.tsv` and `BLUPs-diallel_formated-scaled.tsv`. It applies genomic prediction for scenario 2 with RR method. For each trait, it saves predicted individual genotypic values (`GP_S2_LASSO_pred.tsv`), predictive ability within each half-diallel cross (`GP_S2_LASSO_per-cross-summary.tsv`) and estimated marker effects (`GP_S2_LASSO_mrk-effects.tsv`).

*  `genom-pred_inter-pop_S2WW_RR.Rmd`. This script takes as inputs `geno-p279_formated.tsv`, `BLUPs-p279_formated-scaled.tsv`, `geno-diallel_formated.tsv`, `BLUPs-diallel_formated-scaled.tsv` and `origin_CC279_UTF8.tsv`. It applies genomic prediction for scenario 2 with RR method, with only individuals from WW subpopulation as training set. For each trait, it saves predicted individual genotypic values (`GP_S2WW_RR_pred.tsv`), predictive ability within each half-diallel cross (`GP_S2WW_RR_per-cross-summary.tsv`) and estimated marker effects (`GP_S2WW_RR_mrk-effects.tsv`).


### Analysis of Mendelian sampling predictive ability

* `analyze_results_genom-pred_diallel.Rmd`. This script loads predictive abilities and predicted individual genotypic values from all scenarios and methods (e.g. `GP_S1a/1b/2_RR/LASSO_per-cross-summary.tsv` and `GP_S1a/1b/2_RR/LASSO_pred-obs.Rdata`). It combines results for all scenarios, traits and crosses for the best methods are in `table_results_best-pred-ability_all-scenarios_per-cross.tsv`. It generates Figure 4, S6, S7, S8 and S9.

* `compute_ensemble-genom-pred_diallel.Rmd`. This script loads predicted individual genotypic values from all scenarios and methods (`GP_S1a/1b/2_RR/LASSO_pred-obs.Rdata`) and it computes the best method and ensemble prediction (i.e., the averaged predicted value between RR and the LASSO) at each cross-validation repetition. It generates the following tables: `PA_S1a_per-cross-trait_RR-LASSO-best-ens.tsv`, `PA_S1b_per-cross-trait_RR-LASSO-best-ens.tsv` and `PA_S2_per-cross-trait_RR-LASSO-best-ens.tsv`.


* `analyze_marker-effects_p279-diallel.Rmd`. This script i) analyzes estimated marker effects in both populations and for both methods, ii) studies the proportion of homozygous loci within each cross of the half-diallel and links this to estimated BLUPs of genotypic values and predictive ability. It takes as input estimated marker effects (`GP_S1a/1b/2_RR/LASSO_mrk-effects.tsv`), `GP_S2_RR_per-cross-summary.tsv`, `geno-diallel_formated`, `BLUPs_diallel_formated-scaled.tsv` and `geno-p279_formated.tsv`. It generates a table of the proportion of non-segregating loci `diallel-cross_prop-nonsegreg-homo-loci.tsv` (Table S7).


* `compute_multiple-lin-reg_res_pred-ability.Rmd`. This script applies a multiple linear regression to several variables potentially affecting Mendelian sampling predictive ability for each cross and trait. It takes as inputs predictive abilities from RR, LASSO and the best method resulting from `compute_ensemble-genom-pred_diallel.Rmd` script; heritability and fitting information from `diallel_fit-info_blups_chapitre.tsv` and `diallel_H2-per-cross_all-traits.tsv`; the proportion of homozygous loci `diallel-cross_prop-nonsegreg-homo-loci.tsv`; measures of genetic distances generated in `study_p279_diallel_relationship.Rmd`.It generates Figure 5b and `AllResults_genom-pred-3-scenarios_15-traits_10-crosses_correlated-variables.tab`, a table which recapitulates Mendelian sampling prediction results for all scenarios, crosses and traits.
 


## Training population optimization

Applied to scenario S2 only.

* `Training-pop_p279_opti_diallel_STPGA.Rmd`. This script takes as input files `geno-diallel_formated.tsv`, `geno-p279_formated.tsv` and `origin_CC279_UTF8.tsv`. It applies training set optimization with three methods (CDmean, PEVmean and mean additive relationship), for prediction both inthe whole half-diallel population and within each cross of the half-diallel. As a control, it also applies random sampling of cultivars from the diversity panel. The optimized training sets (i.e., subsets of the diversity panel) are saved in three files: `training-pop-p279_non-optimized_STPGA-varsize.Rdata` (control random sampling), `training-pop-p279_optimized-diall-all_STPGA-varsize.Rdata` (optimization for prediction in whole half-diallel), and `training-pop-p279_optimized-diall-cross_STPGA-varsize.Rdata` (optimization for prediction within each cross).

* `genom-pred_inter-pop_S2_RR_optTRS.Rmd`. This script takes as inputs files generated in `prepare_files_inter-pop_pred_diallel.Rmd` and in `Training-pop-p279_opti_diallel_STPGA.Rmd`.It applies genomic prediction with RR method on the optimized training sets previously created and plots results. Results from each selected training set (i.e., each cross, whole half-diallel and random sampling) are respectively in `GP_S2_RR_optTRS_per-cross-summary.tsv`, `GP_S2_RR_optTRS_all-cross-summary.tsv` and `GP_S2_RR_randTRS_per-cross-summary.Rdata`.

* `genom-pred_inter-pop_S2_LASSO_optTRS.Rmd`.  This script takes as inputs files generated in `prepare_files_inter-pop_pred_diallel.Rmd` and in `Training-pop-p279_opti_diallel_STPGA.Rmd`.It applies genomic prediction with LASSO method on the optimized training sets previously created and plots results. Results from each selected training set (i.e., each cross, whole half-diallel and random sampling) are respectively in `GP_S2_LASSO_optTRS_per-cross-summary.tsv`, `GP_S2_LASSO_optTRS_all-cross-summary.tsv` and `GP_S2_LASSO_randTRS_per-cross-summary.Rdata`.

* `genom-pred_inter-pop_S2_bestMeth_optTRS.Rmd`. This script loads results from `genom-pred_inter-pop_S2_RR_optTRS.Rmd` and `genom-pred_inter-pop_S2_LASSO_optTRS.Rmd` and results of genomic prediction without training set optimization in `table_results_best-pred-ability_all-scenarios_per-cross.tsv`. It selects the best method in terms of predictive ability, both for Mendelian sampling genomic prediction and cross mean prediction (see below). This script generates the Figure S11.


## Genomic prediction of cross means


* `pred_mean_cross_parental-average_diallel.Rmd`. This script takes as inputs `geno-diallel_formated.tsv`, `BLUPs_diallel_formated-scaled.tsv`, estimated marker effects in all scenarios and for both methods in `GP_S1a/1b/2_RR_mrk-effects.tsv`. It computes the mean predicted genotypic value per cross, based on realized or parental average genotypes and previously estimated marker effects (with RR and LASSO methods). It ouputs table for parental average genotypes (`diallel_parental-average-genotypes.tsv`), tables for mean observed and predicted genotypic value based on existing and parental average genotypes, for each cross for all scenarios (`S1a-1b-2_mean.obs.mean.pred-per-cross.tsv`, `pred_cross-mean_S1a-1b-2_parent-average-diallel-cross`), observed and predicted cross mean for the best method between RR and LASSO based on per-cross predictive ability (`pred_cross-mean_S1a-1b-2_best-method-cross_parent-average-diallel-cross.tsv`) or per-trait predictive ability `pred_cross-mean_S1a-1b-2_best-method-cross_parent-average-diallel-trait.tsv`.

* `study_cross-mean_prediction.Rmd`. This script takes as input predicted genotypic values for the prediction of cross mean (`S1a-1b-2_mean.obs.mean.pred-per-cross.tsv`, `pred_cross-mean_S1a-1b-2_parent-average-diallel-cross.tsv`), per-cross and per-trait best method for observed and predicted cross mean (`pred_cross-mean_S1a-1b-2_best-method-cross_parent-average-diallel-cross.tsv`, `pred_cross-mean_S1a-1b-2_best-method-cross_parent-average-diallel-trait.tsv`). It computes predictive ability for the prediction of cross mean, on cross and trait basis. Then, it analyzes variables affecting either per-trait or per-cross predictive ability of cross mean, with same parameters tested as in `compute_multiple-lin-reg_res_pred-ability.Rmd` script. Finally, this script analyzes cross mean prediction with optimized training sets with results generated in `genom-pred_inter-pop_S2_RR_optTRS.Rmd` and `genom-pred_inter-pop_S2_LASSO_optTRS.Rmd`. This script generates Figures 2, 5a, 6, S5, S10, S16, and Table S3.

* `compute_crosses_among-p279.Rmd`. This script loads marker effects estimated in the whole diversity panel and computes parental average genotypes for all possible crosses among cultivars of this population. It generates predicted genotypic values of cross means for all possible crosses. Then, it performs a PCA on all these values on which half-diallel parents predicted values are projected. This script generates Figure S13.


# Bibliography



*This Markdown file can be converted into html with this command: `pandoc -s --toc --from markdown --to html 0_ReadMe.md --citeproc -o 0_ReadMe.html`*

