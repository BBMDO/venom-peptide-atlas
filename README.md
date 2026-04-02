A Structure-Informed Atlas of Venom-Derived Peptides Reveals the Organization of Chemical Space

Authors:
Thaís Caroline Gonçalves
Eduardo Henrique Toral Cortez
Danilo Trabuco Amaral

------------------------------------------------------------
OVERVIEW
------------------------------------------------------------

This repository contains the scripts used to generate the 
structure-informed atlas of venom-derived peptides described 
in the manuscript.

The workflow integrates sequence-based and structure-based 
descriptors to construct a unified feature space and perform 
multivariate analyses of peptide chemical space.

The pipeline includes:

1. Sequence preprocessing and filtering
2. Physicochemical feature extraction
3. Structural feature extraction from AlphaFold models
4. Data integration and quality control
5. Multivariate analyses (PCA, UMAP, clustering)
6. Comparative analyses (venom vs non-venom datasets)

All scripts are designed for reproducibility and can be 
executed independently or as part of an integrated workflow.

------------------------------------------------------------
FILES DESCRIPTION
------------------------------------------------------------

filter_peptides.py
- Cleans FASTA sequences and removes invalid entries
- Computes physicochemical properties using Biopython:
  molecular weight, pI, GRAVY, charge at pH 7, instability index
- Outputs:
  - Clean FASTA file
  - Table with computed properties
  - Log of discarded sequences

param_alphafold.py
- Extracts structural descriptors from PDB/mmCIF files
- Compatible with AlphaFold models
- Computes:
  - Secondary structure (DSSP)
  - Solvent accessible surface area (FreeSASA)
  - Hydrophobic surface fraction
  - Radius of gyration
  - Contact-based features
- Outputs a feature table for downstream analysis

protein_attribute_analysis_workflow.R
- Main pipeline for atlas construction
- Integrates sequence and structure descriptors
- Performs:
  - Data cleaning and standardization
  - Feature selection and QC
  - Summary tables for classes and taxa
  - Correlation analysis (Spearman)
  - PCA and UMAP projections
  - Clustering (k-means, HDBSCAN)
  - PERMANOVA and statistical tests
- Generates publication-ready figures and tables

venom_vs_nonvenom_feature_comparison.R
- Comparative analysis between venom and non-venom peptides
- Identifies shared feature space
- Performs:
  - Balanced sampling (equal N and length-matched)
  - PCA and UMAP visualization
  - Random Forest classification
  - PERMANOVA testing
  - Univariate statistical comparisons
- Outputs comparative figures and importance metrics

------------------------------------------------------------
INPUT DATA
------------------------------------------------------------

The pipeline requires:

1. FASTA files (for sequence preprocessing)
2. Tabular feature files (.tsv/.csv) for descriptor analysis
3. Structural models (PDB or mmCIF format)

Expected formats:

- Tabular files should contain one peptide per row
- Column names are automatically standardized
- Numeric features are automatically detected and filtered

------------------------------------------------------------
OUTPUT
------------------------------------------------------------

The scripts generate:

- Cleaned datasets
- Feature matrices
- Summary tables (CSV/TSV)
- Multivariate analysis outputs
- Publication-ready figures (PDF)
- Supplementary tables

All outputs are written to user-defined directories.

------------------------------------------------------------
REPRODUCIBILITY
------------------------------------------------------------

All analyses were performed using:

R (>= 4.2)
Python (>= 3.8)

Main R packages:
tidyverse, vegan, uwot, randomForest, ComplexHeatmap

Main Python dependencies:
Biopython, freesasa, numpy, pandas

Random seeds are fixed where applicable to ensure 
reproducibility.

------------------------------------------------------------
NOTES
------------------------------------------------------------

- Structural analyses rely on AlphaFold model confidence.
  Low-confidence models should be filtered prior to use.
- Feature selection includes automatic removal of:
  - invariant variables
  - columns with excessive missing values
- The workflow is modular and can be adapted for other 
  peptide datasets.

------------------------------------------------------------
CONTACT
------------------------------------------------------------

For questions or requests:
Danilo Trabuco Amaral
Universidade Federal do ABC (UFABC)
Brazil
