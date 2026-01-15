# DEX_neural_models

This repository contains the analysis scripts developed for the study: **"Glucocorticoid Signaling in PSC-Derived Neural Systems to Elucidate Mechanisms of Stress-Induced Psychiatric Vulnerability"** (under review in *Molecular Neurobiology*).

## Scientific Context & Data Source
The datasets analyzed in this project were obtained from previously published studies. This work focuses on the comparative and integrative analysis of existing transcriptomic data to investigate the effects of glucocorticoids (DEX) across diverse human neural models.

**Data Attribution:**
The raw and processed files in the `/data` folder are curated versions of supplemental data from the following publications:
* **Babaniyi et al. (2023)** - [DOI: 10.1016/j.scr.2023.103086]
* **Cruceanu et al. (2022)** - [DOI: 10.1176/appi.ajp.2021.21010095]
* **Dony et al. (2025)** - [DOI: 10.1126/sciadv.adn8631]
* **Krontira et al. (2024)** - [DOI: 10.1016/j.neuron.2024.02.005]

*Note: For original experimental designs and primary data generation, please refer to the studies cited above.*

## Project Overview
This repository provides the Python pipeline used to integrate these diverse datasets, specifically for:
1. Differential gene expression (DGE) comparison.
2. Functional enrichment mapping (Pathways).
3. Pearson correlation across different neural models (Acute vs. Chronic exposure).

## Repository Structure
* **`correlation_matrix.py`**: Computes Pearson correlation across multiple datasets.
* **`heatmap.py`**: Generates a clustered heatmap linking genes to functional pathways.
* **`scatter.py`**: Visualizes acute vs. chronic log2FC correlations.
* **`volcano.py`**: Plots gene expression significance.
* **`/data`**: Curated Excel files from the cited studies.
* **`/results`**: Directory for output plots (images).

## How to Use
1. **Prerequisites**: You need Python installed with the following libraries:
   `pip install pandas seaborn matplotlib openpyxl adjustText scipy`
2. **Setup**: Ensure the Excel files from the primary studies are placed in the `/data` folder using the naming convention expected by the scripts.
3. **Execution**: Run any script from the root directory to generate the corresponding visualization in the `/results` folder.
* **`volcano.py`**: Plots gene expression significance.
* **`/data`**: Curated Excel files from the cited studies.
* **`/results`**: Directory for output plots.
