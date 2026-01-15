import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ─── CONFIG ──────────────────────────────────────────────
input_file = "data/gene_functional_terms_mapping_w_fc.xlsx"
output_img = "results/clustermap_vias_topgenes_com_codigo.png"

TOP_VIAS = 50
TOP_GENES = 50
V_LIM = 3

# ─── LEITURA ─────────────────────────────────────────────
df = pd.read_excel(input_file, engine="openpyxl")

df["gene"] = df["gene"].astype(str).str.upper()
df["functional_term"] = df["functional_term"].astype(str)
df["term_code"] = df["term_code"].astype(str)

# Label final no eixo Y
df["term_label"] = df["term_code"] + ": " + df["functional_term"]

# ─── SELEÇÃO DAS TOP VIAS ────────────────────────────────
top_terms = (
    df["term_label"]
    .value_counts()
    .head(TOP_VIAS)
    .index
    .tolist()
)

df = df[df["term_label"].isin(top_terms)]

# ─── REMOVER DUPLICATAS GENE × VIA (maior |log2FC|) ──────
df["_abs_"] = df["log2FC"].abs()
df = (
    df.sort_values("_abs_", ascending=False)
      .drop_duplicates(["gene", "term_label"])
      .drop(columns="_abs_")
)

# ─── TOP GENES POR |log2FC| ──────────────────────────────
top_gene_list = (
    df.groupby("gene")["log2FC"]
      .apply(lambda x: x.abs().max())
      .nlargest(TOP_GENES)
      .index
      .tolist()
)

df = df[df["gene"].isin(top_gene_list)]

# ─── MATRIZ TERMOS × GENES ───────────────────────────────
heatmap_df = df.pivot(
    index="term_label",
    columns="gene",
    values="log2FC"
)

heatmap_df.dropna(axis=0, how="all", inplace=True)
heatmap_df.dropna(axis=1, how="all", inplace=True)

clustering_data = heatmap_df.fillna(0)

# ─── CLUSTERMAP ──────────────────────────────────────────
sns.set(style="white", font_scale=0.9)
cmap = sns.diverging_palette(240, 10, as_cmap=True)

g = sns.clustermap(
    clustering_data,
    cmap=cmap,
    vmin=-V_LIM,
    vmax=V_LIM,
    center=0,
    figsize=(16, 14),
    linewidths=0.05,
    linecolor="white",
    cbar_kws={"label": "log2FC"},
    xticklabels=True,
    yticklabels=True
)

# ─── AJUSTES VISUAIS ─────────────────────────────────────
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

g.ax_heatmap.set_xlabel("Regulated genes", fontsize=10)
g.ax_heatmap.set_ylabel("Enriched functional pathways", fontsize=10)

g.cax.set_ylabel("log2FC", fontsize=7)
g.cax.tick_params(labelsize=6)

plt.subplots_adjust(top=0.94)

# ─── SALVAR ──────────────────────────────────────────────
plt.savefig(output_img, dpi=600, bbox_inches="tight")
plt.show()

print("✓ Gráfico salvo em:", output_img)
