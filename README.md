RetiPhenoAge downstream analyses
Overview
This repository contains R scripts used for downstream statistical analyses in the study:
From Human Retina to Health Risk: AI-enabled Oculomic Aging Predicts Phenomic Morbidities and Reveals Proteomic Signatures
RetiPhenoAge is a previously developed and externally validated deep learning–based retinal biological age biomarker. No AI model development or training was performed in this study. The analyses in this repository evaluate associations between RetiPhenoAge and disease outcomes, quantitative traits, and plasma proteins.

Contents
1_PheWAS_UKB.R
Phenome-wide association analyses in UK Biobank.
Main procedures:
•	Incident disease analyses using Cox proportional hazards models.
•	Exclusion of participants with prevalent diagnoses before baseline.
•	Estimation of hazard ratios and confidence intervals.
•	C-index and time-dependent discrimination analyses.
•	Calibration analyses.
•	Proportional hazards assumption testing using Schoenfeld residuals.
•	Restricted cubic spline analyses.
Disease outcomes were defined using phecode-based classifications derived from linked hospital inpatient and mortality records.
 
2_Quantiative_traits_UKB.R
Cross-sectional analyses of quantitative traits in UK Biobank.
Main procedures:
•	Multivariable linear regression.
•	Associations between RetiPhenoAge (per SD increment) and quantitative traits.
•	Adjustment for demographic, lifestyle, socioeconomic, medication, and clinical covariates.
•	False discovery rate correction within trait categories.
 
3_Proteomic_UKB.R
Olink proteomic data preprocessing.
Main procedures:
•	Quality control filtering.
•	Missing-value filtering.
•	Panel-specific processing.
•	Yeo–Johnson normalization.
•	Multiple imputation using Amelia.
•	Generation of processed proteomic datasets for downstream analyses.
•	RetiPhenoAge and proteomic association analyses
 
4_SEED.R
External replication analyses in the Singapore Epidemiology of Eye Diseases (SEED) cohort.
Main procedures:
•	Exclusion of participants with baseline disease.
•	Logistic regression analyses for selected incident outcomes.
•	Estimation of odds ratios and confidence intervals.
•	Ethnicity interaction analyses.
•	Ethnicity-stratified subgroup analyses.
Selected outcomes include diabetes, hypertension, chronic kidney disease, and hyperlipidaemia.
 
Data availability
Participant-level UK Biobank and SEED data are not included in this repository because redistribution is restricted by data access agreements and ethical approvals.
UK Biobank data are available through application to UK Biobank:
https://www.ukbiobank.ac.uk/
SEED data may be requested through the Singapore Eye Research Institute subject to institutional and ethical approval requirements.
 
Reproducibility
The scripts contain placeholder file names and paths that should be modified by approved data users.
Input datasets, phenotype definitions, and covariate variables must be generated from approved UK Biobank and SEED data resources before running the scripts.
 
Software
Analyses were conducted in R.
Key packages include:
•	survival
•	rms
•	pROC
•	epitools
•	dplyr
•	purrr
•	ggplot2
•	Amelia
•	bestNormalize
 
Citation
Please cite the associated manuscript when using or referencing this repository.
