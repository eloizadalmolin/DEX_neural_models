import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from adjustText import adjust_text

# =========================================================
# CONFIGURATION & PATHS
# =========================================================
base_path = "data"
results_path = "results"

# Ensure the results directory exists
if not os.path.exists(results_path):
    os.makedirs(results_path)

files = {
    "Babaniyi_H1": "babaniyih1.xlsx",
    "Babaniyi_H9": "babaniyih9.xlsx",
    "Cruceanu": "cruceanu_dn.xlsx",
    "Dony": "donysemfc.xlsx",
    "Krontira": "krontira_dn.xlsx"
}

meta = {
    "Babaniyi_H1": ("hgnc_symbol", ["log2FoldChange"]),
    "Babaniyi_H9": ("hgnc_symbol", ["log2FoldChange"]),
    "Cruceanu": ("gene", ["log2FC"]),
    "Dony": ("gene", ["log2FC_Line409b2", "log2FC_LineFOK4"]),
    "Krontira": ("gene", ["log2FoldChange"])
}

# =========================================================
# DATA PROCESSING FUNCTIONS
# =========================================================
def expand_long_format(df, gene_col, fc_cols, tag):
    """Filters and expands data frames into long format."""
    if tag == "Cruceanu" and "model" in df.columns:
        df = df[df["model"].str.strip().str.upper() == "NEURONS"]
    if tag == "Dony" and "model" in df.columns:
        df = df[df["model"].isin(["ChP", "Ex.Neurons", "Imm.ChP", "Inh.Neurons", "RGS5Neurons"])]
        
    df = df[[gene_col] + fc_cols].copy()
    df.columns = ["gene"] + fc_cols
    df[fc_cols] = df[fc_cols].apply(pd.to_numeric, errors='coerce')
    df = df[df[fc_cols].notna().any(axis=1)].copy()
    df["gene"] = df["gene"].astype(str).str.upper().str.strip()
    
    long_frames = []
    for col in fc_cols:
        col_name = f"{tag}_{col}" if len(fc_cols) > 1 else tag
        long_frames.append(df[["gene", col]].rename(columns={col: "logFC"}).assign(dataset=col_name))
    return pd.concat(long_frames, ignore_index=True)

# =========================================================
# DATA LOADING & MERGING
# =========================================================
frames = []
for tag, fname in files.items():
    gene_col, fc_cols = meta[tag]
    path = os.path.join(base_path, fname)
    df = pd.read_excel(path, engine="openpyxl")
    frames.append(expand_long_format(df, gene_col, fc_cols, tag))

full_df = pd.concat(frames, ignore_index=True)

# Pivot table for gene × dataset matrix
wide_df = full_df.pivot_table(index="gene", columns="dataset", values="logFC", aggfunc="first")

# Define Acute and Chronic groups
acute_cols = ["Babaniyi_H1", "Babaniyi_H9", "Cruceanu"]
chronic_cols = ["Dony_log2FC_Line409b2", "Dony_log2FC_LineFOK4", "Krontira"]

# Select the most extreme log2FC for each group
wide_df["acute"] = wide_df[acute_cols].apply(
    lambda row: row.dropna().loc[abs(row.dropna()).idxmax()] if not row.dropna().empty else np.nan, axis=1)
wide_df["chronic"] = wide_df[chronic_cols].apply(
    lambda row: row.dropna().loc[abs(row.dropna()).idxmax()] if not row.dropna().empty else np.nan, axis=1)

# Valid subset for correlation analysis
subset_df = wide_df[["acute", "chronic"]].dropna()
r, pval = pearsonr(subset_df["acute"], subset_df["chronic"])
print(f"📊 Pearson Correlation (extreme values per gene): r = {r:.3f}, p = {pval:.2e}")

# =========================================================
# ANNOTATION LOGIC
# =========================================================
def top_genes_by_quadrant_distance(df, top_n=4):
    """Finds the most extreme genes by distance from origin in each quadrant."""
    df = df.copy()
    df["dist"] = np.sqrt(df["acute"]**2 + df["chronic"]**2)

    quadrants = {
        "Q1": df[(df["acute"] > 0) & (df["chronic"] > 0)],
        "Q2": df[(df["acute"] < 0) & (df["chronic"] > 0)],
        "Q3": df[(df["acute"] < 0) & (df["chronic"] < 0)],
        "Q4": df[(df["acute"] > 0) & (df["chronic"] < 0)],
    }

    selected = []
    for q_name, qdf in quadrants.items():
        if not qdf.empty:
            top_genes = qdf.nlargest(top_n, "dist")
            selected.append(top_genes)

    return pd.concat(selected).drop(columns="dist")

# Select genes for labeling
labels_df = top_genes_by_quadrant_distance(subset_df, top_n=4)

# Ensure most extreme chronic genes are also labeled
extreme_chronic_genes = subset_df.loc[
    [subset_df["chronic"].idxmax(), subset_df["chronic"].nlargest(2).index[1],
     subset_df["chronic"].idxmin(), subset_df["chronic"].nsmallest(2).index[1]]
]

# =========================================================
# VISUALIZATION
# =========================================================
plt.figure(figsize=(8, 8))
plt.scatter(subset_df["acute"], subset_df["chronic"], alpha=0.6, s=30, edgecolor='black')

# Add reference lines
lims = [subset_df.min().min() - 0.2, subset_df.max().max() + 0.2]
plt.plot(lims, lims, '--', color='grey')
plt.axhline(0, color='grey', ls='--', lw=0.8)
plt.axvline(0, color='grey', ls='--', lw=0.8)
plt.xlim(lims)
plt.ylim(-1, 1) # Specific Y-axis limit for this analysis

plt.xlabel("log₂FC (Acute: Babaniyi_H1, H9, Cruceanu)")
plt.ylabel("log₂FC (Chronic: Dony_409b, FOK4, Krontira)")
plt.grid(ls=":", lw=0.5)

# Annotation process
texts = []
for _, row in labels_df.iterrows():
    texts.append(plt.text(row["acute"], row["chronic"], row.name, fontsize=8, color='black'))

for _, row in extreme_chronic_genes.iterrows():
    texts.append(plt.text(row["acute"], row["chronic"], row.name, fontsize=8, color='black'))

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))

plt.tight_layout()

# Save figure
output_path = os.path.join(results_path, "acute_vs_chronic_scatter.png")
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.show()

print(f"✓ Scatter plot successfully saved to: {output_path}")