# Drug Side-Effect and Target Analysis  
Final Project – Group 22  
R for Bio Data Science (DTU)

This repository contains a complete workflow for acquiring, cleaning, merging, and analyzing drug side-effect data. The project integrates two major data sources (SIDER and DrugBank) and investigates how drug targets and molecular properties relate to observed clinical side effects.

---

## Repository Structure

```text
R/
  00_all.qmd            Run entire pipeline
  01_load.qmd           Download SIDER and DrugBank data
  02_clean.qmd          Clean and merge datasets
  03_augment.qmd        Add physicochemical descriptors
  05_analysis_1.qmd     Global analyses (side effects, targets, indications)
  05_analysis_2.qmd     Targeted analyses (hormonal and stress-axis receptors)
  06_analysis_2.qmd     Descriptor distributions and side-effect burden
  99_proj_func.R        Helper functions

data/                   Cleaned data outputs  
results/                Plots and analysis figures
```
---

# Pipeline Summary
## 1. Loading (01_load.qmd)

Downloads SIDER and DrugBank datasets, assigns column names, and exports raw tables. DrugBank files require user credentials.

## 2. Cleaning (02_clean.qmd)

Standardizes SIDER tables (side effects, frequencies, indications) and merges DrugBank vocabulary, structures, and target data. Produces cleaned datasets for both sources.

## 3. Augmentation (03_augment.qmd)

Adds molecular descriptors from SMILES:

Molecular weight

H-bond donors and acceptors

cLogP

Lipinski rule indicator

Merges SIDER and DrugBank into a single dataset (merged_data.tsv).

## 4. Global Analyses (05_analysis_1.qmd)

Examines:

Frequency of side effects

Correlation between reporting frequency and number of reports

Co-occurrence patterns (indications vs. side effects, protein targets vs. side effects)

Normalized heatmaps using average frequency rates

## 5. Targeted Receptor Analyses (05_analysis_2.qmd)

Tests expected biological relationships:

Drugs targeting hormone-related receptors show higher proportions of endocrine side effects.

Drugs targeting stress-axis receptors show higher proportions of metabolic side effects.

Corresponding figures are saved in results/.

## 6. Descriptor and Burden Analysis (06_analysis_2.qmd)

Aggregates side effects per drug and evaluates relationships between:

Molecular descriptors (HBD, HBA, MW, cLogP)

Lipinski compliance

Total number of unique side effects (side-effect burden)

Includes distribution plots and scatter/boxplots.

---

#Reproducibility

Required packages include:
tidyverse, purrr, readr, stringr, curl, rcdk, and RChemMass.

To rerun the full pipeline:
```r
quarto render R/00_all.qmd
```
---

## Contributors:
#### Malte Lau (s224183) - Malteflau
#### Duco Lam (252126) - DucoLam
#### Benedek Dunavolgyi (s243161) - adunavolgyi
#### Manuel Charneca (s253708) - Manuelidk


Group 22
R for Bio Data Science, DTU
