Project-1-Dissertation-2024

Overview
This repository contains the R code used for the analyses in my dissertation entitled "Multi-Omics Biomarker Discovery: Using AI to Integrate Gene Expression and DNA Methylation Data in Lung Cancer".
The code is organized into sections that correspond to the methods and materials section of the dissertation.

Repository Structure
The R code is divided into several scripts, each corresponding to a specific section of the analysis. Below is an overview of the sections and their corresponding scripts:

2.3	RNA-Seq Data Preprocessing	= Step 1: Automated mRNA Filtering.R

2.4	Optimal Clustering Analysis =	Step 2: Optimal Clustering Algorithm.R

2.5	K-means Clustering Analysis	= Step 3: K-means Clustering.R

(Optional) Hierarchical Clustering	= Step 3.1: Hierarchical Clustering.R

2.5.1	RNA-Seq Heatmap	= Step 4: Gene Expression Heatmap.R

2.6	Survival Analysis = Step 5: Survival Analysis (OS and DFS).R

2.7	Clinical Association Tests = Step 6: Association Tests.R

2.8	Pathway Analysis (Differential Gene Analysis) =	Step 7: Pathway Analysis (Differential Gene Analysis).R

2.9	Differential Methylation Analysis	= Step 8: Differential Methylation Analysis.R

2.10	Differential Protein Analysis =	Step 9: Protein Differential Expression Analysis.R

3.11.1	Model 1 – RNA-Seq Model	= Step 10.1: Model 1 RNA-Seq Data.R

3.11.2	Model 2 – Methylation Model	= Step 10.2: Model 2 Methylation Data.R

3.11.3	Model 3 – Integrated Model (RNA-Seq combined with Methylation Data)	= Step 10.3: Model 3 RNA-Seq and Methylation Data.R

3.11.4	Variable Importance	Variable = Importance Model 1.R, Variable Importance Model 2.R, Variable Importance Model 3.R

3.11.5	Box Plot Visualisation	= Variable Importance Boxplot.R

3.11.6	Heatmap Visualisation =	Model 1 Heatmap.R, Model 2 Heatmap.R, Model 3 Heatmap.R

3.11.7	Gene Ontology (GO) Enrichment Analysis = Gene Ontology and Model Retraining Model 1.R, Gene Ontology and Model Retraining Model 2.R, Gene Ontology and Model Retraining Model 3.R

3.11.8	Model Retraining with Important Variables = 	Gene Ontology and Model Retraining Model 1.R, Gene Ontology and Model Retraining Model 2.R, Gene Ontology and Model Retraining Model 3.R

How to Use

Clone the repository:
git clone https://github.com/Laurielou26/Project-1-Dissertation-2024.git
Open the project in RStudio.
Run the scripts in the order specified in the table above.

Data Source

The data used in this project was obtained from cBioPortal for Cancer Genomics.
Data source link : https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018 

Feel free to modify this as needed! 
