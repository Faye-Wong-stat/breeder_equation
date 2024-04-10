---
title: ReadMe
author: Brault Charlotte
bibliography: /home/brault/Documents/Articles/Article2/dataverse/biblio_zotero.bib
link-citations: true
date: "" 
---


This document describes the input files, scripts and output files generated for the publication "Interest of phenomic prediction as an alternative to genomic prediction in grapevine". This document is in Markdown format.


Table in **.tsv** format were automatically converted into **.tab** by data.inrae.fr.

*This Markdown file can be converted into html with this command: `pandoc -s --toc --from markdown --to html 0_ReadMe.md --citeproc -o 0_ReadMe.html`*


# Input data 

## Genotypic values and marker data

Genotypic values for phenotypic data and molecular markers were already available in an other dataverse: (dataverse)[https://doi.org/10.15454/PNQQUQ] .

* `origin_CC279.xlsx`. This table classifies diversity panel genotypes into the three subpopulations: WW (wine west), WE (wine east) and TE (table east), based on both the original SSR classification by [@nicolasGeneticDiversityLinkage2016] and the SNP classification by [@flutreGenomewideAssociationPrediction2020].

* `ids_279.xlsx`. This table provides diverse cultivar characteristics for genotypes from the diversity panel. 

* `diallel_layout_coord.tsv`. This table provides field coordinates of half-diallel trial in Villeneuve-lès-Maguelone.

## Spectra

### Half-diallel population

* `diallel_wood_2020.zip`. This directory contains all spectra collected on half-diallel wood samples in 2020. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-bois_final.xlsx` table.

* `diallel_wood_2021.zip`. This directory contains all spectra collected on half-diallel wood samples in 2021. Spectra are splitted in two files, one for each block. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-bois_2021.xlsx` table.

* `diallel_leaves_2020.zip`. This directory contains all spectra collected on half-diallel leaf samples in 2020. Spectra are splitted in two files, one for each block. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-feuilles_final.xlsx` table.

* `diallel_leaves_2021.zip`. This directory contains all spectra collected on half-diallel leaf samples in 2021. Spectra are splitted in two files, one for each block, and one supplementary file for spectra that have been remade. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-feuille_2021.xlsx` table.

### Diversity panel population

* `p279_wood_2020.zip`. This directory contains all spectra collected on diversity panel wood samples in 2020. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-bois_final.xlsx` table.

* `p279_wood_2021.zip`. This directory contains all spectra collected on diversity panel wood samples in 2021. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-bois_2021.xlsx` table.

* `p279_leaves_2020.zip`. This directory contains all spectra collected on diversity panel leaf samples in 2020. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-feuilles_final.xlsx` table.

* `p279_leaves_2021.zip`. This directory contains all spectra collected on diversity panel leaf samples in 2021. The correspondence between spectrum number and genotype information is in `tableau manip NIRS-feuille_2021.xlsx` table.

### Correspondence files

* `tableau manip NIRS-bois_final.xlsx`. This table makes the correspondence between spectrum number for **wood** tissue and **2020** and genotype information (genotype name, block, coordinates in the field).

* `tableau manip NIRS-feuilles_final.xlsx`. This table makes the correspondence between spectrum number for **leaf** tissue and **2020** and genotype information (genotype name, block, coordinates in the field).

* `tableau manip NIRS-bois_2021.xlsx`. This table makes the correspondence between spectrum number for **wood** tissue and **2021** and genotype information (genotype name, block, coordinates in the field).

* `tableau manip NIRS-feuille_2021.xlsx`. This table makes the correspondence between spectrum number for **leaf** tissue and **2021** and genotype information (genotype name, block, coordinates in the field).

# Scripts

## Source scripts

* `get_colours_diallel.R`. This is a small R-script with a function defining a specific color for each half-diallel cross, panel subpopulation and half-diallel parent.

* `half-diallel_analysis_functions.R`. This is a R-script gathering useful functions for half-diallel analysis.

* `ortho_nirs.R`. This is a R-script containing *ortho.nirs* function for computing orthogonalization of NIRS data, using REP-ASCA method.

* `PS_functions.R`. This is a R-script containing useful functions to read ASD files from spectrometer.

* `computeRV.R`. This is a R-script containing the function *computeRV* which calculates the RV coefficient from two matrices.

## Spectra data analysis

### Half-diallel population

* `wood_spectrum_pre-process_diallel_2020.Rmd`. This script takes as input `diallel_wood_2020` and `tableau manip NIRS-bois_final.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum, remove outlier spectra, average the four spectra per plot), and pre-processes for the half-diallel population and 2020 wood spectra. It generates a file containing all pre-process spectra at plot level: `wood_diallel_2020_nirs-list-2-blocks.Rdata` and spectra at genotype level: `wood_diallel_2020_all_spectra.Rdata`.

* `wood_spectrum_pre-process_diallel_2021.Rmd`. This script takes as input `diallel_wood_2021` and `tableau manip NIRS-bois_2021.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum, remove outlier spectra, average the four spectra per plot), and pre-processes for the half-diallel population and 2021 wood spectra. It generates a file containing all pre-process spectra at plot level: `wood_diallel_2021_nirs-list-2-blocks.Rdata` and spectra at genotype level: `wood_diallel_2021_all_spectra.Rdata`.

* `leaves_spectrum_pre-process_diallel_2020.Rmd`. This script takes as input `diallel_leaves_2020` and `tableau manip NIRS-feuilles-final.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum, remove outlier spectra, average the four spectra per plot), and pre-processes for the half-diallel population and 2020 leaves spectra. It generates a file containing all pre-process spectra at plot level: `leaves_diallel_2020_nirs-list-2-blocks.Rdata` and spectra at genotype level: `leaves_diallel_2020_all_spectra.Rdata`.

* `leaves_spectrum_pre-process_diallel_2021.Rmd`. This script takes as input `diallel_leaves_2021` and `tableau manip NIRS-feuille_2021.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum, remove outlier spectra, average the four spectra per plot), and pre-processes for the half-diallel population and 2021 leaves spectra. It generates a file containing all pre-process spectra at plot level: `leaves_diallel_2021_nirs-list-2-blocks.Rdata` and spectra at genotype level: `leaves_diallel_2021_all_spectra.Rdata`.

### Diversity panel population

* `wood_spectrum_pre-process_p279_2020.Rmd`. This script takes as input `p279_wood_2020`, `tableau manip NIRS-bois_final.xlsx`, `origin_CC279.xlsx` and `ids_279.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum), and pre-processes for the diversity panel population and 2020 wood spectra. It generates a file containing all pre-process spectra at plot level: `wood_p279_2020_nirs-list-2-blocks.Rdata` and spectra at genotype level: `wood_p279_2020_all_spectra.Rdata`.

* `wood_spectrum_pre-process_p279_2021.Rmd`. This script takes as input `p279_wood_2021`, `tableau manip NIRS-bois_2021.xlsx`, `origin_CC279.xlsx` and `ids_279.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum), and pre-processes for the diversity panel population and 2021 wood spectra. It generates a file containing all pre-process spectra at plot level: `wood_p279_2021_nirs-list-2-blocks.Rdata` and spectra at genotype level: `wood_p279_2021_all_spectra.Rdata`.

* `leaves_spectrum_pre-process_p279_2020.Rmd`. This script takes as input `p279_leaves_2020`, `tableau manip NIRS-feuilles_final.xlsx`, `origin_CC279.xlsx` and `ids_279.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum), and pre-processes for the diversity panel population and 2020 leaves spectra. It generates a file containing all pre-process spectra at plot level: `leaves_p279_2020_nirs-list-2-blocks.Rdata` and spectra at genotype level: `leaves_p279_2020_all_spectra.Rdata`.

* `leaves_spectrum_pre-process_p279_2021.Rmd`. This script takes as input `p279_leaves_2021`, `tableau manip NIRS-feuilles_final.xlsx`, `origin_CC279.xlsx` and `ids_279.xlsx`. It applies corrections (remove of jump, cut at the end of the spectrum), and pre-processes for the diversity panel population and 2021 leaves spectra. It generates a file containing all pre-process spectra at plot level: `leaves_p279_2021_nirs-list-2-blocks.Rdata` and spectra at genotype level: `leaves_p279_2021_all_spectra.Rdata`.

### Combine spectra and extract spectra BLUP for both populations


* `compile_spectra_for-phenom-pred.Rmd`. This script combines spectra from both tissues and years and it deduces the common set of genotypes for each population. Then, spectra in each modality for these common sets are saved in `nirs_all-pre-process_diallel-p279_common-genos.Rdata`. Genotypic values from phenotypic data for these common sets are saved in `diallel_blups_geno-plus-cross-wide_chapitre_common-nirs.tsv` and `p279_blups_geno_chapitre_common-nirs`.


* `extract_BLUPs-geno_spectra_lme4.Rmd`. This script takes as input `nirs_all-pre-process_diallel-p279_common-genos.Rdata`, and all spectra averaged at plot level such as `wood_diallel_2020_nirs-list-2-blocks.Rdata`. It applies different mixed models, combining one or several years and one or several tissues, for each population. For each model tested, it generates three files: one with spectra BLUP, one with heritability along the wavelengths, and variance components along the wavelengths. Genetic effect was modelled with **a genotype effect**, i.e., no information about cross or subpopulation was added in the model.

* `extract_BLUPs-genocross_spectra_lme4.Rmd`. This script takes as input `nirs_all-pre-process_diallel-p279_common-genos.Rdata`, and all spectra averaged at plot level such as `wood_diallel_2020_nirs-list-2-blocks.Rdata`. It applies different mixed models, combining one or several years and one or several tissues, for each population. For each model tested, it generates three files: one with spectra BLUP, one with heritability along the wavelengths, and variance components along the wavelengths. Genetic effect was modelled with **genotype and cross effects, these effects were summed** to constitute spectra BLUPs.

* `extract_BLUPs-genocross-geno_spectra_lme4.Rmd`. This script takes as input `nirs_all-pre-process_diallel-p279_common-genos.Rdata`, and all spectra averaged at plot level such as `wood_diallel_2020_nirs-list-2-blocks.Rdata`. It applies different mixed models, combining one or several years and one or several tissues, for each population. For each model tested, it generates three files: one with spectra BLUP, one with heritability along the wavelengths, and variance components along the wavelengths. Genetic effect was modelled with **genotype and cross effects, but only genotype effect** was used to constitute spectra BLUPs.

### Analyze and combine spectra for both populations

* `analyze_BLUP_spectra.Rmd`. This script takes as input variance components and heritability values from `extract_BLUPs-genocross_spectra_lme4.Rmd` script. It generates the Figure 1 and Figure S1.

* `study_coinertia_geno_rel-matrices.Rmd`. This script takes as input marker data from both diversity panel and half-diallel populations: `geno-p279_formated.tsv` and `geno-diallel_formated.tsv`; genotypic values from both populations: `BLUPs-p279_formated-scaled.tsv` and `BLUPs-diallel_formated-scaled.tsv`; spectra BLUPs from the script `extract_BLUPs-geno_spectra_lme4.Rmd`. It computes PCA analysis for a given year and pre-process and save PCA components.

* `study_coinertia_genocross_rel-matrices.Rmd`. This script takes as input marker data from both diversity panel and half-diallel populations: `geno-p279_formated.tsv` and `geno-diallel_formated.tsv`; genotypic values from both populations: `BLUPs-p279_formated-scaled.tsv` and `BLUPs-diallel_formated-scaled.tsv`; spectra BLUPs from the script `extract_BLUPs-genocross_spectra_lme4.Rmd`. It computes PCA analysis for a given year and pre-process and save PCA components.

* `study_coinertia_genocross-geno_rel-matrices.Rmd`. This script takes as input marker data from both diversity panel and half-diallel populations: `geno-p279_formated.tsv` and `geno-diallel_formated.tsv`; genotypic values from both populations: `BLUPs-p279_formated-scaled.tsv` and `BLUPs-diallel_formated-scaled.tsv`; spectra BLUPs from the script `extract_BLUPs-genocross-geno_spectra_lme4.Rmd`. It computes PCA analysis for a given year and pre-process and save PCA components.

* `study_multi-year-coinertia_rel-matrices.Rmd`. This script takes as input PCA outputs generated in the three previous scripts coinertia analysis. It computes co-inertia analysis for different years, BLUP models and for both populations. It generates Figure 2, Figure S2 and S3.


## Fit phenomic and genomic predictions

* `GP-PP_intra-pop_pred_glmnet` and `GP-PP_intra-pop_pred_lme4GS` for rrBLUP and GBLUP/HBLUP methods, respectively. This script is associated with `gpp.sh` bash script for launching genomic or phenomic prediction for all modalities: six pre-processes, two tissues, three years and two populations. It takes as input genotypic values for the 15 traits studied (`BLUPs-p279_formated-scaled.tsv` or `BLUPs-diallel_formated-scaled.tsv`) and marker or phenomic data. Results for all modalities are saved separately in a common directory.

* `GP-PP_intra-pop_pred_2mat-lme4GS.Rmd`. This script is associated with `gpp2mat.sh` for launching multi-matrix predictions for with two variance-covariance matrices for different modalities: SNPs + wood NIRS, SNPs + leaf NIRS or wood NIRS + leaf NIRS. It takes as input genotypic values for the 15 traits studied (`BLUPs-p279_formated-scaled.tsv` or `BLUPs-diallel_formated-scaled.tsv`) and marker or phenomic data. Results for all modalities are saved separately in a common directory.

* `GP-PP_intra-pop_pred_3mat-lme4GS.Rmd`. This script is associated with `gpp3mat.sh` for launching multi-matrix predictions for with three variance-covariance matrices for different modalities: SNPs + wood NIRS + leaf NIRS. It takes as input genotypic values for the 15 traits studied (`BLUPs-p279_formated-scaled.tsv` or `BLUPs-diallel_formated-scaled.tsv`), marker and phenomic data. Results for all modalities are saved separately in a common directory.

* `gpp.sh`. This is a bash script associated with `GP-PP_intra-pop_pred_glmnet.Rmd` and `GP-PP_intra-pop_pred_lme4GS.Rmd` for launching several scripts in a SLURM cluster.

* `gpp2mat.sh`. This is a bash script associated with `GP-PP_intra-pop_pred_2mat-lme4GS` for launching several scripts in a SLURM cluster.

* `gpp3mat.sh`. This is a bash script associated with `GP-PP_intra-pop_pred_3mat-lme4GS` for launching several scripts in a SLURM cluster.

## Analyze results

* `prelim-analysis_results_phenom-pred.Rmd`. This script takes as input results of phenomic and genomic predictions. It compares rrBLUP and GBLUP/HBLUP methods, as well as using base or BLUP spectra. It generates the Figure 3 and S6.

* `analyze_genom-phenom-pred_results.Rmd`. This script takes as input results of phenomic and genomic predictions. It combines these results and save them into tables: `results_GP-PP.1mat-RR-glmnet_p279.subpop_15-traits.tsv`, `results_GP-PP.1mat-RR-glmnet_diallel.cross_15-traits.tsv`, `results_GPP.2mat-lme4GS_p279.subpop_15-traits.tsv`, `results_GPP.2mat-lme4GS_diallel.cross_15-traits.tsv`, `results_GPP.3mat-lme4GS_p279.subpop_15-traits.tsv` and `results_GPP.3mat-lme4GS_diallel.cross_15-traits.tsv`. It compares phenomic prediction modalities (pre-processes, tissues, years) and generates Figure S4.

* `compare_genom-phenom-pred_results.Rmd`. This script takes as input the tables `results_GP-PP.1mat-RR-glmnet_p279.subpop_15-traits.tsv`, `results_GP-PP.1mat-RR-glmnet_diallel.cross_15-traits.tsv`, `results_GPP.2mat-lme4GS_p279.subpop_15-traits.tsv`, `results_GPP.2mat-lme4GS_diallel.cross_15-traits.tsv`, `results_GPP.3mat-lme4GS_p279.subpop_15-traits.tsv` and `results_GPP.3mat-lme4GS_diallel.cross_15-traits.tsv`. It generates the Figure 5, 6, S7 and S8.




# Bibliography





