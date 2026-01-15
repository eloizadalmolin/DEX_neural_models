import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# =========================================================
# CONFIGURATION & PATHS
# =========================================================
base_path = "data"
results_path = "results"

# Ensure the results directory exists
if not os.path.exists(results_path):
    os.makedirs(results_path)

# Input files mapping
files = {                        
    "babaniyi_H1": "babaniyih1.xlsx",
    "babaniyi_H9": "babaniyih9.xlsx",
    "cruceanu_dn": "cruceanu_dn.xlsx",
    "dony": "donysemfc.xlsx",
    "krontira_dn": "krontira_dn.xlsx"
}

# Metadata for gene columns and Fold Change columns
meta = {                        
    "babaniyi_H1": ("hgnc_symbol", ["log2FoldChange"]),
    "babaniyi_H9": ("hgnc_symbol", ["log2FoldChange"]),
    "cruceanu_dn": ("gene", ["log2FC"]),
    "dony": ("gene", ["log2FC_Line409b2", "log2FC_LineFOK4"]),
    "krontira_dn": ("gene", ["log2FoldChange"]),
}

dony_models = ["ChP", "Ex.Neurons", "Imm.ChP", "Inh.Neurons", "RGS5Neurons"]

# =========================================================
# PROCESSING FUNCTIONS
# =========================================================
def expand_long_format_and_count(df, gene_col, fc_cols, tag):
    """Filters, cleans, and prepares the data for correlation."""
    if tag == "cruceanu_dn" and "model" in df.columns:
        df = df[df["model"].str.strip().str.upper() == "NEURONS"]
        print(f"📌 Cruceanu filtered for 'Neurons', {df.shape[0]} genes.")
    elif tag == "dony" and "model" in df.columns:
        df = df[df["model"].isin(dony_models)]
        print(f"📌 Dony filtered for models {dony_models}, {df.shape[0]} genes.")

    df = df[[gene_col] + fc_cols].copy()
    df.columns = ["gene"] + fc_cols
    df[fc_cols] = df[fc_cols].apply(pd.to_numeric, errors='coerce')
    df = df[df[fc_cols].notna().any(axis=1)].copy()
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()

    long_frames = []
    counts = {}
    for col in fc_cols:
        col_name = (
            "Dony_409b" if col == "log2FC_Line409b2" else
            "Dony_FOK4" if col == "log2FC_LineFOK4" else
            "Babaniyi_H9" if tag == "babaniyi_H9" else
            "Babaniyi_H1" if tag == "babaniyi_H1" else
            "Krontira" if tag == "krontira_dn" else
            "Cruceanu" if tag == "cruceanu_dn" else tag
        )
        subset = df[["gene", col]].dropna()
        long_frames.append(subset.rename(columns={col: "logFC"}).assign(dataset=col_name))
        counts[col_name] = subset.shape[0]
    return pd.concat(long_frames, ignore_index=True), counts

# =========================================================
# DATA LOADING & MERGING
# =========================================================
frames = []
counts_dict = {}

for tag, fname in files.items():
    gene_col, fc_cols = meta[tag]
    path = os.path.join(base_path, fname)
    
    # Reading Excel files from the /data folder
    df = pd.read_excel(path, engine="openpyxl")
    expanded_df, counts = expand_long_format_and_count(df, gene_col, fc_cols, tag)
    frames.append(expanded_df)
    counts_dict.update(counts)

# Gene counts summary
counts_df = pd.DataFrame(list(counts_dict.items()), columns=["Dataset", "Num_genes"])
print("\n📌 Number of genes considered per dataset/lineage:")
print(counts_df)

# Concatenate all data into a wide matrix
full_df = pd.concat(frames, ignore_index=True)
wide_df = full_df.pivot_table(index='gene', columns='dataset', values='logFC', aggfunc='first')

# =========================================================
# STATISTICAL ANALYSIS: PEARSON CORRELATION
# =========================================================
corr_matrix = wide_df.corr(method='pearson')

print("\n📌 Pearson Correlation Matrix:")
print(corr_matrix)

# =========================================================
# EXPORTING RESULTS
# =========================================================
# Save Correlation Matrix to Excel
corr_out_path = os.path.join(results_path, "acute_vs_chronic_correlation.xlsx")
corr_matrix.to_excel(corr_out_path)
print(f"\n✓ Correlation matrix saved to: {corr_out_path}")

# Save Correlation Heatmap to PNG
heatmap_out_path = os.path.join(results_path, "acute_vs_chronic_heatmap.png")
plt.figure(figsize=(8,6))
sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", square=True)
plt.tight_layout()
plt.savefig(heatmap_out_path, dpi=600)
plt.show()

print(f"✓ Heatmap saved to: {heatmap_out_path}")