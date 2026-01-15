# DEX_neural_models

Analysis scripts for transcriptomic and functional enrichment analyses of glucocorticoid (DEX) exposure in human neural models.

## 📌 Project Overview
This repository contains Python scripts developed to analyze differential gene expression (DGE) and functional enrichment data, comparing acute vs. chronic exposure effects in human neural cell lines.

## 📂 Repository Structure
* **`correlation_matrix.py`**: Computes Pearson correlation across multiple datasets.
* **`heatmap.py`**: Generates a clustered heatmap of functional pathways.
* **`scatter.py`**: Visualizes acute vs. chronic log2FC correlations.
* **`volcano.py`**: Plots gene expression significance and fold change.
* **`/data`**: Input Excel files for all analyses.
* **`/results`**: Directory where generated plots are automatically saved.

## 🚀 Getting Started
1. **Requirements**: You need Python installed with the following libraries:
   `pip install pandas seaborn matplotlib openpyxl adjustText scipy`
2. **Execution**: Ensure your Excel data is in the `/data` folder and run any script, for example:
   `python heatmap.py`
