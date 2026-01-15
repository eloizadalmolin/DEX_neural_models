import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ─── CONFIGURATION & PATHS ────────────────────────────────
input_file = "data/gene_functional_terms_mapping_w_fc.xlsx"
output_img = "results/clustermap_functional_pathways.png"

# Analysis Parameters
TOP_PATHWAYS = 50
TOP_GENES = 50
COLOR_LIMIT = 3  # V_LIM for the heatmap scale

# ─── DATA LOADING ─────────────────────────────────────────
if not os.path.exists("results"):
    os.makedirs("results")

df = pd.read_excel(input_file, engine="openpyxl")

# Data Cleaning
df["gene"] = df["gene"].astype(str).str.upper()
df["functional_term"] = df["functional_term"].astype(str)
df["term_code"] = df["term_code"].astype(str)

# Final Y-axis label format
df["term_label"] = df["term_code"] + ": " + df["functional_term"]

# ─── TOP PATHWAYS SELECTION ───────────────────────────────
top_terms = (
    df["term_label"]
    .value_counts()
    .head(TOP_PATHWAYS)
    .index
    .tolist()
)

df = df[df["term_label"].isin(top_terms)]

# ─── REMOVE GENE × PATHWAY DUPLICATES (Keep highest |log2FC|) 
df["_abs_"] = df["log2FC"].abs()
df = (
    df.sort_values("_abs_", ascending=False)
      .drop_duplicates(["gene", "term_label"])
      .drop(columns="_abs_")
)

# ─── TOP GENES SELECTION BY |log2FC| ──────────────────────
top_gene_list = (
    df.groupby("gene")["log2FC"]
      .apply(lambda x: x.abs().max())
      .nlargest(TOP_GENES)
      .index
      .tolist()
)

df = df[df["gene"].isin(top_gene_list)]

# ─── MATRIX GENERATION (Terms × Genes) ────────────────────
heatmap_df = df.pivot(
    index="term_label",
    columns="gene",
    values="log2FC"
)

heatmap_df.dropna(axis=0, how="all", inplace=True)
heatmap_df.dropna(axis=1, how="all", inplace=True)

# Fill missing values with 0 for clustering
clustering_data = heatmap_df.fillna(0)

# ─── CLUSTERMAP GENERATION ────────────────────────────────
sns.set(style="white", font_scale=0.9)
cmap = sns.diverging_palette(240, 10, as_cmap=True)

g = sns.clustermap(
    clustering_data,
    cmap=cmap,
    vmin=-COLOR_LIMIT,
    vmax=COLOR_LIMIT,
    center=0,
    figsize=(16, 14),
    linewidths=0.05,
    linecolor="white",
    cbar_kws={"label": "log2FC"},
    xticklabels=True,
    yticklabels=True
)

# ─── VISUAL ADJUSTMENTS ───────────────────────────────────
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

g.ax_heatmap.set_xlabel("Regulated Genes", fontsize=10)
g.ax_heatmap.set_ylabel("Enriched Functional Pathways", fontsize=10)

g.cax.set_ylabel("log2FC", fontsize=7)
g.cax.tick_params(labelsize=6)

plt.subplots_adjust(top=0.94)

# ─── SAVE AND DISPLAY ─────────────────────────────────────
plt.savefig(output_img, dpi=600, bbox_inches="tight")
plt.show()

print(f"✓ Heatmap successfully saved to: {output_img}")
