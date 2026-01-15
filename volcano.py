import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text
import matplotlib.patheffects as pe

# =========================
# PATHS
# =========================
input_file = "data/gene_functional_terms_mapping_w_fc_pval.xlsx"
output_file = "results/volcano_plot_key_genes_top10.png"

# =========================
# LOAD DATA
# =========================
df = pd.read_excel(input_file)

df["gene"] = df["gene"].astype(str).str.upper()
df["p_value"] = df["p_value"].replace(0, 1e-300)
df["neg_log10_pval"] = -np.log10(df["p_value"])

# =========================
# REGULATION
# =========================
df["regulation"] = "Not significant"
df.loc[(df["log2FC"] > 1) & (df["p_value"] < 0.05), "regulation"] = "Up"
df.loc[(df["log2FC"] < -1) & (df["p_value"] < 0.05), "regulation"] = "Down"

# =========================
# Y AXIS VISUAL LIMITS
# =========================
Y_MAX = 300
Y_MIN = 0
MARGIN = 8

df["neg_log10_pval_plot"] = df["neg_log10_pval"].clip(
    lower=Y_MIN + MARGIN,
    upper=Y_MAX - MARGIN
)

# =========================
# PLOT
# =========================
plt.figure(figsize=(10, 8))

colors = {
    "Up": "red",
    "Down": "blue",
    "Not significant": "gray"
}

for reg, color in colors.items():
    subset = df[df["regulation"] == reg]
    plt.scatter(
        subset["log2FC"],
        subset["neg_log10_pval_plot"],
        color=color,
        alpha=0.6,
        s=30,
        edgecolor="black",
        linewidth=0.5,
        label=reg
    )

# =========================
# CUTOFF LINES
# =========================
plt.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=0.8)
plt.axvline(1, color="black", linestyle="--", linewidth=0.8)
plt.axvline(-1, color="black", linestyle="--", linewidth=0.8)

# =========================
# TOP 10 UP / DOWN
# =========================
top_up = (
    df[df["regulation"] == "Up"]
    .sort_values("log2FC", ascending=False)
    .drop_duplicates("gene")
    .head(10)
)

top_down = (
    df[df["regulation"] == "Down"]
    .sort_values("log2FC", ascending=True)
    .drop_duplicates("gene")
    .head(10)
)

# =========================
# KEY CONSISTENT GENES
# =========================
KEY_GENES = [
    "FKBP5", "HIF3A", "RASSF4", "TSC22D3",
    "PIK3R1", "LIFR", "CHST15"
]

key_df = (
    df[df["gene"].isin(KEY_GENES)]
    .loc[
        df[df["gene"].isin(KEY_GENES)]
        .groupby("gene")["log2FC"]
        .apply(lambda x: x.abs().idxmax())
    ]
)

# =========================
# MERGE LABELS
# =========================
label_df = (
    pd.concat([top_up, top_down, key_df])
    .drop_duplicates("gene")
)

# =========================
# ANNOTATIONS
# =========================
texts = []

for _, row in label_df.iterrows():

    x = row["log2FC"]
    y = min(
        max(row["neg_log10_pval"], Y_MIN + MARGIN),
        Y_MAX - MARGIN
    )

    if row["regulation"] == "Up":
        text_color = "red"
    elif row["regulation"] == "Down":
        text_color = "blue"
    else:
        text_color = "gray"

    texts.append(
        plt.text(
            x,
            y,
            row["gene"],
            fontsize=8,
            color=text_color,
            ha="center",
            va="bottom",
            path_effects=[
                pe.withStroke(linewidth=2, foreground="white")
            ]
        )
    )

adjust_text(
    texts,
    arrowprops=dict(arrowstyle="-", color="black", lw=0.6),
    expand_points=(1.6, 1.8),
    expand_text=(1.6, 1.8),
    force_text=1.0,
    force_points=0.6
)

# =========================
# FINAL STYLE
# =========================
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(p-value)")
plt.legend()
plt.grid(ls=":", lw=0.5)
plt.ylim(Y_MIN - MARGIN, Y_MAX + MARGIN)
plt.tight_layout()

# =========================
# SAVE FIGURE (PNG)
# =========================
plt.savefig(
    output_file,
    dpi=600,
    bbox_inches="tight"
)

plt.show()

print("✓ Volcano plot saved as:")
print(output_file)
