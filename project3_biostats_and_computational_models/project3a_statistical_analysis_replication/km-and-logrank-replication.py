"""
Script: tex_signature_survival_GSE14520.py
Author: Juliana Patiño Gallego

Description:
    This script replicates the survival analysis of the 13-gene
    TEX-related exhaustion signature using the GSE14520 hepatocellular
    carcinoma dataset (platform GPL3921).

    The workflow includes:
        1. Loading and parsing the GEO series matrix (expression data)
        2. Reading the GPL3921 annotation file and mapping probe IDs to gene symbols
        3. Extracting the 13-gene TEX-related signature and collapsing probe-level
           expression to gene-level values
        4. Loading an additional clinical metadata file and extracting survival
           time and status
        5. Merging expression and clinical data by GSM sample ID
        6. Computing a TEX risk score using published Cox coefficients
        7. Stratifying patients into High- and Low-risk groups (median split)
        8. Running Kaplan–Meier survival analysis and log-rank test
        9. Plotting KM curves with confidence intervals and number-at-risk table

Inputs:
    - GSE14520_GPL3921_series_matrix.txt.gz
        Gene expression matrix (GEO format)
    - GPL3921.annot
        Platform annotation file mapping probes to gene symbols
    - GSE14520_Extra_Supplement.txt.gz
        Additional clinical metadata with survival information

Signature:
    13-gene TEX-related exhaustion signature:
        HSPD1, UBB, DNAJB4, CALM1, LGALS3,
        BATF, COMMD3, IL7R, FDPS, DRAP1,
        RPS27L, PAPOLA, GPR171

Output:
    - Kaplan–Meier survival plot:
        * km_GSE14520_replicated.png
    - Printed tables describing:
        * Expression matrix dimensions
        * Annotation mappings
        * Missing signature genes (if any)
        * Risk scores and risk group counts
        * Survival data preview and event distribution

Notes / assumptions:
    * GEO series matrices contain metadata lines beginning with "!", which
      are automatically skipped using `comment="!"` in pandas.
    * GPL annotation tables may start at a variable line; this script
      searches for the correct header automatically.
    * Probe-to-gene mapping may not be one-to-one; probe-level data
      are collapsed to gene-level by averaging across probes.
    * Survival status is parsed in a flexible manner to capture variations
      such as "dead", "deceased", "1", etc.
    * Risk groups are defined using the median signature score, following
      the original paper.
    * Kaplan–Meier analysis uses the lifelines library, including
      log-rank testing and number-at-risk tables.

"""

# ==== Library imports ====
import gzip
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines.statistics import logrank_test


# ==== Helper to open plain or gzip-compressed files ====
def smart_open(path, mode="rt"):
    """
    Open a text file that may be plain or gzip-compressed.

    Parameters
        path : str
            Path to the file.
        mode : str, optional
            Open mode (text). Default is "rt".

    Returns
        file object
            A file handle that can be iterated or read as usual.

    Notes
    -----
    * If the file ends with '.gz', it is opened with gzip.
    * Uses UTF-8 decoding and ignores problematic characters.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode=mode, encoding="utf-8", errors="ignore")
    else:
        return open(path, mode=mode, encoding="utf-8", errors="ignore")


# ==== Load expression matrix (GSE14520, platform GPL3921) ====

series_path = "GSE14520_GPL3921_series_matrix.txt.gz"  # adjust path if needed

# The series_matrix table appears after GEO header lines starting with "!".
# Using comment="!" tells pandas to skip those metadata lines automatically.
df_raw = pd.read_csv(
    series_path,
    sep="\t",
    comment="!"
)

# First column corresponds to the probe ID (Affymetrix probes)
df_raw = df_raw.rename(columns={df_raw.columns[0]: "ProbeID"})
df_raw = df_raw.set_index("ProbeID")

# Convert expression values to numeric (non-numeric entries become NaN)
df_raw = df_raw.apply(pd.to_numeric, errors="coerce")

# ==== Load GPL3921 annotation and build mapping: gene_symbol -> list of probe IDs ====

gpl_path = "GPL3921.annot"

# ---- Find the header line where the actual table starts ----
header_line_idx = None
with smart_open(gpl_path, "rt") as f:
    for i, line in enumerate(f):
        line_stripped = line.strip()
        # GEO platform files typically start the table with "ID" or "ID_REF"
        if line_stripped.startswith("ID\t") or line_stripped.startswith("ID_REF\t"):
            header_line_idx = i
            break

if header_line_idx is None:
    raise RuntimeError("Could not find a header line (ID/ID_REF) in GPL3921 annotation file.")

# ---- Read the annotation table from that header line ----
gpl = pd.read_csv(
    gpl_path,
    sep="\t",
    header=header_line_idx,
    dtype=str,
    compression="gzip" if gpl_path.endswith(".gz") else "infer"
)

# Normalize column names: lowercase and replace spaces with underscores
gpl.columns = gpl.columns.str.lower().str.replace(" ", "_")

def get_col(df, candidates):
    """
    Return the first column name in `df` that matches any
    of the strings in `candidates`.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame whose columns will be searched.
    candidates : list of str
        Possible column names (already normalized).

    Returns
    -------
    str
        The first matching column name.

    Raises
    ------
    KeyError
        If none of the candidate names are present.
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of the candidate columns found: {candidates}")


# Typical columns in a GPL file:
#  * "ID"  (or "ID_REF") for probe ID
#  * "Gene Symbol" for gene symbol
probe_col = get_col(gpl, ["id", "id_ref", "probe_id"])
symbol_col = get_col(gpl, ["gene_symbol", "genesymbol", "symbol"])

# ---- Build gene_symbol -> list of probe IDs ----
gpl_sub = gpl[[symbol_col, probe_col]].dropna()
gpl_sub = gpl_sub[gpl_sub[symbol_col] != ""]

symbol_to_probe = (
    gpl_sub
    .groupby(symbol_col)[probe_col]
    .apply(list)
)


# ==== 13-gene TEX-related signature from the paper ====

signature_genes = [
    "HSPD1", "UBB", "DNAJB4", "CALM1", "LGALS3",
    "BATF", "COMMD3", "IL7R", "FDPS", "DRAP1",
    "RPS27L", "PAPOLA", "GPR171"
]

# Expand to all probe IDs associated with these 13 genes
signature_probes = []
for gene in signature_genes:
    if gene in symbol_to_probe.index:
        signature_probes.extend(symbol_to_probe[gene])
    else:
        # Some genes may not have probes on this platform
        print(f"WARNING: gene {gene} not found in GPL3921 annotation.")

# Subset the expression matrix to those probes only
expr_probes = df_raw.loc[df_raw.index.intersection(signature_probes)]


# ==== Collapse probes to gene-level expression (mean per gene) ====

# Build a reverse mapping: probe ID -> gene symbol
probe_to_gene = {}
for gene in signature_genes:
    if gene in symbol_to_probe.index:
        for probe in symbol_to_probe[gene]:
            probe_to_gene[probe] = gene

# Rename index from probe IDs to gene symbols and average
# expression if multiple probes map to the same gene.
expr_genes = (
    expr_probes
    .rename(index=probe_to_gene)
    .groupby(level=0)
    .mean()
)

# Transpose to get samples as rows and genes as columns
expr_genes = expr_genes.transpose()


# ==== Load extra clinical file (GSE14520_Extra_Supplement) and extract survival time + status ====

extra_path = "GSE14520_Extra_Supplement.txt.gz"  # adjust path if needed

extra = pd.read_csv(
    extra_path,
    sep="\t",
    dtype=str,
    compression="gzip"
)

def find_col_by_substring(df, substrings):
    """
    Locate the first column whose name contains any of the
    given substrings (case-insensitive).

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame whose columns will be searched.
    substrings : list of str
        List of substrings to search for.

    Returns
    -------
    str
        Name of the first matching column.

    Raises
    ------
    KeyError
        If no column name contains any of the substrings.
    """
    cols = df.columns
    lower = cols.str.lower()
    for sub in substrings:
        sub = sub.lower()
        mask = lower.str.contains(sub)
        if mask.any():
            return cols[mask][0]
    raise KeyError(f"No column containing any of {substrings}")


# ---- Identify the sample ID column (GSM accession, etc.) ----
id_col = find_col_by_substring(extra, ["geo_accession", "gsm", "sample", "id"])

# Use that column as index so it matches GSM IDs in the expression matrix
extra = extra.set_index(extra[id_col])

# ---- Select survival time and status columns (known names for this file) ----
time_col = "Survival months"
status_col = "Survival status"

time_raw = extra[time_col]
status_raw = extra[status_col]

# Extract numeric values from "Survival months" (e.g. "54.3 months" -> 54.3)
time_clean = time_raw.str.extract(r"(\d+\.?\d*)", expand=False)
time = pd.to_numeric(time_clean, errors="coerce")

# Map status to event indicator:
# 1 = death, 0 = alive / censored
status_clean = status_raw.str.lower()
event = status_clean.isin(["dead", "deceased", "1"]).astype(int)

clinical_surv = pd.DataFrame(
    {"time": time, "event": event},
    index=extra.index
)


# ==== Merge expression (13 genes) with clinical survival data ====

# expr_genes index: GSM IDs from series matrix
# clinical_surv index: IDs from extra (hopefully same GSM IDs)
merged = expr_genes.join(clinical_surv, how="inner")

# Drop samples with missing survival information
merged = merged.dropna(subset=["time", "event"])

# Ensure that 'time' and 'event' are numeric (safety step)
merged["time"] = pd.to_numeric(merged["time"], errors="coerce")
merged["event"] = pd.to_numeric(merged["event"], errors="coerce")

# If any rows became NaN after numeric conversion, drop them
merged = merged.dropna(subset=["time", "event"])


# ==== Compute TEX-related 13-gene risk score (as in the paper) ====

# Cox coefficients for the 13-gene TEX-related signature
coeffs = {
    "HSPD1": 0.32788,
    "UBB": -0.72135,
    "DNAJB4": 0.33201,
    "CALM1": 0.38101,
    "LGALS3": 0.33354,
    "BATF": 0.23593,
    "COMMD3": 0.49528,
    "IL7R": -0.25727,
    "FDPS": 0.31791,
    "DRAP1": 0.28689,
    "RPS27L": 0.25298,
    "PAPOLA": 0.24165,
    "GPR171": -0.57949
}

# Check if any genes in the signature are missing from the merged matrix
missing_genes = [g for g in coeffs if g not in merged.columns]
if missing_genes:
    print("\nWARNING: these genes are missing in merged expression:", missing_genes)

# Compute risk score as a linear combination:
# risk_score = sum( expression_gene_i * coef_i )
merged["risk_score"] = merged.apply(
    lambda row: sum(row[g] * coeffs[g] for g in coeffs if g in merged.columns),
    axis=1
)


# ==== Define high vs low risk groups (median split) ====

# Median cut-off (simplest and symmetric choice)
median_score = merged["risk_score"].median()
merged["risk_group"] = np.where(merged["risk_score"] >= median_score, "High", "Low")

# Convert survival time from months to years for more intuitive X-axis
merged["time_years"] = merged["time"] / 12.0


# ==== Kaplan–Meier analysis + number-at-risk table ====

# Global plot aesthetics: classic "journal" style
custom_params = {'axes.spines.right': False, 'axes.spines.top': False}
sns.set_theme(style='ticks', rc=custom_params)

kmf_low = KaplanMeierFitter()
kmf_high = KaplanMeierFitter()

# Subsets for each risk group
low = merged[merged["risk_group"] == "Low"]
high = merged[merged["risk_group"] == "High"]

# --- Fit Kaplan–Meier curves for each group ---
kmf_low.fit(
    durations=low["time_years"],
    event_observed=low["event"],
    label="Low-risk"
)
kmf_high.fit(
    durations=high["time_years"],
    event_observed=high["event"],
    label="High-risk"
)

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Plot KM curves with confidence intervals
kmf_low.plot(
    ax=ax, ci_show=True, ci_alpha=0.2,
    lw=2, color="#1f77b4"
)
kmf_high.plot(
    ax=ax, ci_show=True, ci_alpha=0.2,
    lw=2, color="#d62728"
)

# Axis labels and title
ax.set_title("GSE14520 – Kaplan–Meier by 13-gene TEX-related risk score")
ax.set_xlabel("Time (years)")
ax.set_ylabel("Survival probability")

# Y-axis from 0 to 1 with ticks every 0.25
ax.set_ylim(0, 1)
ax.set_yticks(np.arange(0, 1.01, 0.25))

# X-axis ticks from 0 to 5 years (as in the original paper)
ax.set_xticks(np.arange(0, 6, 1))

# --- Log-rank test comparing High vs Low risk groups ---
lr = logrank_test(
    low["time_years"], high["time_years"],
    event_observed_A=low["event"],
    event_observed_B=high["event"]
)

# Annotate p-value (4 decimal places) inside the plot
ax.text(
    0.05, 0.15,
    f"Log-rank p = {lr.p_value:.4f}",
    transform=ax.transAxes,
    fontsize=11,
    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7)
)

# --- Add number-at-risk table under the KM plot ---
# Only show the "At risk" row for each group,
# aligned with the chosen X-axis ticks.
add_at_risk_counts(
    kmf_high, kmf_low,
    ax=ax,
    rows_to_show=["At risk"],
    xticks=np.arange(0, 6, 1)  # 0, 1, 2, 3, 4, 5 years
)

plt.tight_layout()

fig = plt.gcf()
fig.savefig("km_GSE14520_replicated.png", format="png", bbox_inches="tight", dpi=300)

plt.show()
