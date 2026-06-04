# CICIV

A Causality-Based Framework for Identifying Drivers of Breast Cancer Progression

## Overview

This repository provides the reference implementation of **CICIV**, a causality-based framework for identifying potential driver genes associated with breast cancer progression.

CICIV is designed to address two major challenges in transcriptomic causal analysis: high-dimensional gene expression data and potential confounding bias. The framework integrates causal feature selection with causal effect estimation to identify genes that may have causal effects on breast cancer progression rather than merely showing statistical correlations.

The source code and data processing scripts in this repository support the main analyses described in the manuscript:

**A Novel Causality-Based Method for Identifying Drivers of Breast Cancer Progression**

## Framework Design

CICIV follows a two-stage causal inference workflow:

1. **Causal feature selection using PC-simple**

   The PC-simple module performs conditional independence tests to reduce the dimensionality of high-dimensional gene expression data. It identifies candidate parent genes and preserves potential causal regulatory relationships among genes.

2. **Causal effect estimation using CIV**

   The conditional instrumental variable representation learning module estimates the absolute treatment effect of selected genes while accounting for potential latent confounding factors. The estimated absolute causal effect is used to rank candidate breast cancer driver genes.

Together, these two stages allow CICIV to move from causal structure discovery to quantitative causal effect estimation in transcriptomic data.

## Repository Structure

```text
CICIV/
├── README.md
├── LICENSE
├── requirements.txt
├── brca_metabric.py
├── example_PC.py
├── pc_simple.py
├── BAMB.py
├── FBEDk.py
├── GSMB.py
├── IAMB.py
├── IAMBnPC.py
├── KIAMB.py
├── LCMB.py
├── MBOR.py
├── STMB.py
├── HITON/
├── IPCMB/
├── MMMB/
├── PCMB/
├── TIE_star/
├── brca_metabric/
├── common/
├── idea/
└── semi_HITON/
```

## Core Modules

The main components of this repository include:

* `pc_simple.py`: implementation of the PC-simple-based causal feature selection process.
* `brca_metabric.py`: script for external validation using the METABRIC breast cancer dataset.
* `example_PC.py`: example script for PC set discovery.
* `BAMB.py`, `FBEDk.py`, `GSMB.py`, `IAMB.py`, `KIAMB.py`, `LCMB.py`, `MBOR.py`, `STMB.py`: causal feature selection and Markov blanket-related comparison methods.
* `HITON/`, `IPCMB/`, `MMMB/`, `PCMB/`, `TIE_star/`, `semi_HITON/`: algorithm modules used for causal discovery or comparison experiments.

## Data Sources

The main analysis is based on breast cancer transcriptomic data from **TCGA-BRCA**. Gene expression data were preprocessed before causal analysis, including normalization, filtering, and transformation steps.

External validation was performed using the independent **METABRIC** breast cancer dataset. The METABRIC analysis was used to evaluate the cross-dataset robustness of the identified breast cancer causal genes.

Due to data size and database usage policies, full raw datasets are not directly included in this repository. Users should download the original datasets from the corresponding public databases and process them according to the scripts and descriptions provided in this repository.

## Requirements

The code is implemented mainly in Python. To install the required Python packages, run:

```bash
pip install -r requirements.txt
```

Recommended environment:

```text
Python >= 3.8
numpy
pandas
scipy
scikit-learn
```

Additional packages may be required depending on the specific script or experiment.

## Usage

To run an example PC-simple analysis:

```bash
python example_PC.py
```

To run the METABRIC validation script:

```bash
python brca_metabric.py
```

Before running the scripts, please make sure that the input data paths are correctly specified in the corresponding files.

## Output

The main outputs include:

* candidate causal parent sets
* causal feature selection results
* external validation results on METABRIC
* intermediate files for downstream causal effect estimation and comparison analysis

The final results can be used to evaluate the causal relevance of candidate breast cancer driver genes.

## Intended Use

This software is intended for research purposes only. The identified causal genes should be interpreted as computationally inferred candidate drivers and should be further validated through biological experiments and clinical evidence.

Automatically generated results should be reviewed carefully before being used for publication, biological interpretation, or downstream experimental design.

## Citation

If you use this code or build upon this repository, please cite the corresponding manuscript:

```text
A Novel Causality-Based Method for Identifying Drivers of Breast Cancer Progression.
```

## License

This project is released under the MIT License.
