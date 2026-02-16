# Analysis Code: Resilience as a Mediator in Life Event–Depression Pathways

This repository contains the R scripts used for the statistical analyses in the manuscript:

> **"Associations between adaptive life events and depressive symptoms are differentially mediated by resilience across the lifespan: path analyses in a naturalistic inpatient sample with major depression"** > *Under review at the Journal of Affective Disorders.*

## 📌 Project Overview
This study investigates the mediating role of resilience (CD-RISC) in the relationship between lifetime adaptive/adverse life events (TAQ) and depressive symptoms (BDI-II). The repository provides the code to reproduce the data cleaning steps, descriptive statistics, correlation matrices, and the structural equation models (path analyses) used in the study.

The analyses are divided into two parts:
1. **Discovery Sample:** A naturalistic cohort of psychiatric inpatients with Major Depressive Disorder (MDD).
2. **External Validation Sample:** An independent clinical cohort of patients with Persistent Depressive Disorder (PDD) and Borderline Personality Disorder (BPD) to test the robustness of the structural models.

## 📂 Repository Structure

* `01_main_discovery_analysis.R` - Script for data processing, descriptive statistics, and main path models for the MDD discovery sample.
* `02_external_validation_analysis.R` - Script for data processing and external validation of the path model in the independent clinical sample (PDD/BPD).

## 🔒 Data Availability Statement
To protect patient confidentiality and comply with the ethical regulations of the local ethics committee (Faculty of Medicine, LMU Munich), the raw clinical dataset (`.sav` files) is **not publicly shared** in this repository. 

The R scripts are provided here to ensure full transparency of the data transformation, scoring algorithms (e.g., TAQ subscales), and statistical modeling (`lavaan`). Researchers interested in collaborating or accessing the anonymized data should contact the corresponding author to discuss data-sharing agreements.

## 💻 Prerequisites & Dependencies
The analyses were conducted in R. To run these scripts, you will need to install the following R packages:


install.packages(c("tidyverse", "lavaan", "foreign", "psych", "broom", "easystats", "tidyr", "dplyr", "moments", "flextable", "Hmisc", "BBmisc"))


🚀 How to Use the Code
Clone or download this repository.

Open the scripts in R or RStudio.

Because the raw data is not provided, running the read.spss() functions at the top of the scripts will result in an error.

However, you can inspect the code to view the exact specifications of the path models (defined using lavaan syntax) and the custom functions used to calculate the psychometric sum scores and extract bootstrapped standard errors.

📄 License
This code is licensed under the MIT License. You are free to use, modify, and distribute this code for your own research, provided that proper attribution is given to the original authors.

✉️ Contact
For any questions regarding the analysis code or the manuscript, please contact the corresponding author: Prof. Dr. Stephan Goerigk Department of Psychology, Charlotte Fresenius Hochschule, Munich

Email: nico.steffen@charlotte-fresenius-uni.de or stephan.goerigk@charlotte-fresenius-uni.de


