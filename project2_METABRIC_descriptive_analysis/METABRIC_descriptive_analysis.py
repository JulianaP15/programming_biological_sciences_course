"""
Script: METABRIC_descriptive_analysis.py
Author: Juliana Patiño Gallego
Description:
    Exploratory analysis of the METABRIC breast cancer dataset.
    Generates visual summaries relating molecular subtype classifications
    (PAM50 and 3-gene classifier) with patient age, survival, and therapies.

Analyses performed:
    1) Age distribution histogram for all patients and those who died of disease.
    2) Violin plots of overall survival by PAM50 and 3-gene classifications.
    3) Bar plots: proportion of patients receiving each therapy per subtype.
    4) Heatmaps: distribution of the number of therapies per subtype.
    5) Scatterplot: mutation_count vs tumor_size by tumor stage.

Notes / assumptions:
    * Input file must be the METABRIC_RNA_Mutation.csv dataset.
    * Some columns contain mixed datatypes and are cast to 'object' to avoid parsing errors.
    * Therapy columns are binary (1 = received, 0 = not received).
    * Output: Figures are saved as PDF and PNG under the 'figures' folder.
"""

# ==== Library imports ====
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# ==== Functions ====

def save_figure(fig, name, output_dir):
    """
    Save a Matplotlib figure as .pdf and .png

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure object to be saved.
    name : str
        Base name of the output file (without extension).
    output_dir : str, optional
        Path to the output directory where figures will be stored.

    Notes
    -----
    * Uses tight bounding boxes to avoid excessive white margins.
    * Saves vector PDF (for print) and 300 DPI PNG (for web).
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure the output directory exists

    # Define output paths
    path_pdf = os.path.join(output_dir, f"{name}.pdf")
    path_png = os.path.join(output_dir, f"{name}.png")

    # Save both formats
    fig.savefig(path_pdf, format="pdf", bbox_inches="tight")
    fig.savefig(path_png, format="png", bbox_inches="tight", dpi=300)

# ==== Configuration ====
# Input file: METABRIC dataset
filename = 'METABRIC_RNA_Mutation.csv'  # expected in working directory
# Output folder for saved figures (files with the same name will be overwritten)
output_dir = "figures"

# ==== Global aesthetics ====
# Style: light grid, no top/right spines
custom_params = {'axes.spines.right': False, 'axes.spines.top': False}
sns.set_theme(style='ticks', rc=custom_params)

# Base palette
palette_base = sns.color_palette('Set2')  # used for therapy bars

# Consistent color mapping for the analysis
palette_dict = {'All Patients': 'lightblue',
                'Died of Disease': 'salmon',
                'Chemotherapy': palette_base[0],
                'Hormone Therapy': palette_base[1],
                'Radiotherapy': palette_base[2],
                'Tumor Stage': 'viridis',  # colormap for scatter hue
                }

# Define consistent subtype orders from worst to best prognosis
pam50_order = ['Basal','claudin-low','Her2','LumB','LumA','Normal','NC']
three_gene_order = ['ER-/HER2-','ER+/HER2- High Prolif','HER2+','ER+/HER2- Low Prolif']
# NOTE: Ensure these labels exactly match the dataset values (case/spelling).

# ==== 0) Set and Subset used across plots ====

# Load CSV with explicit dtype casting for mixed-type columns
df = pd.read_csv(filename,
                 dtype={678: object, 688: object, 690: object, 692: object})

# Filter and save subset (patients who died of breast cancer)
selection_death = df['death_from_cancer'] == 'Died of Disease'
df_death = df[selection_death]

# ==== 1) Histogram: Age at diagnosis ====
plt.figure(figsize=(9.5, 6.2))

# Distributions
sns.histplot(data=df,
             x='age_at_diagnosis',
             binwidth=2,
             label='All Patients',
             color=palette_dict['All Patients']
             )
sns.histplot(data=df_death,
             x='age_at_diagnosis',
             binwidth=2,
             label='Died of Disease',
             color=palette_dict['Died of Disease']
             )

# Medians
median_all = df['age_at_diagnosis'].median()
median_death = df_death['age_at_diagnosis'].median()

# Median lines
plt.axvline(median_all, color='blue', linestyle='--', linewidth=2, alpha=0.8)
plt.axvline(median_death, color='red', linestyle='--', linewidth=2, alpha=0.8)

# Annotate each line with the median value
plt.text(median_all + 0.5, plt.ylim()[1]*0.9, f'Median (All) = {median_all:.1f}', color='blue')
plt.text(median_death + 0.5, plt.ylim()[1]*0.85, f'Median (Died of Disease) = {median_death:.1f}', color='red')

# Labels and formatting
plt.title('Age Distribution at Diagnosis', pad=15)
plt.xlabel('Age at Diagnosis (years)')
plt.ylabel('Number of Patients')
plt.legend(frameon=False,
           bbox_to_anchor=(1.2, 1),
           )
plt.tight_layout()

# Save figure
fig = plt.gcf()  # get the current active figure
save_figure(fig, '01_age_distribution',output_dir)

plt.show()


# ==== 2) Violin plots: Survival by classification ====
fig, axes = plt.subplots(1,2, figsize=(15, 6.2), sharey=True)

# -- Left panel: PAM50
ax = axes[0]

sns.violinplot(data=df,
               x='pam50_+_claudin-low_subtype',
               y='overall_survival_months',
               order=pam50_order,
               split=True,
               inner='quart',
               color=palette_dict['All Patients'],
               ax=ax
               )
sns.violinplot(data=df_death,
               x='pam50_+_claudin-low_subtype',
               y='overall_survival_months',
               order= pam50_order,
               split=True,
               inner='quart',
               color=palette_dict['Died of Disease'],
               ax=ax
               )

# Medians overlay
medians_all_pam50 = df.groupby('pam50_+_claudin-low_subtype')['overall_survival_months'].median()
medians_death_pam50 = df_death.groupby('pam50_+_claudin-low_subtype')['overall_survival_months'].median()

for i, subtype in enumerate(pam50_order):
    if subtype in medians_all_pam50.index:
        y_val = medians_all_pam50[subtype]
        # Line across the violin (blue = All)
        ax.hlines(y=y_val,
                  xmin=i-0.25,
                  xmax=i+0.25,
                  colors='blue',
                  linestyles='--',
                  linewidth=2)
        # Annotate numeric value slightly above the line
        ax.text(i,
                y_val + 2,
                f'{y_val:.1f}',
                color='blue',
                ha='center',
                va='bottom',
                fontsize=9)
    if subtype in medians_death_pam50.index:
        y_val = medians_death_pam50[subtype]
        # Line across the violin (red = Death)
        ax.hlines(y=y_val,
                  xmin=i-0.25,
                  xmax=i+0.25,
                  colors='red',
                  linestyles='--',
                  linewidth=2)
        ax.text(i,
                y_val + 2,
                f'{y_val:.1f}',
                color='red',
                ha='center',
                va='bottom',
                fontsize=9)

# Formatting
ax.set_title('Overall Survival by PAM50 Classification', pad=12)
ax.set_xlabel('PAM50 Subtype', labelpad=12)
ax.set_ylabel('Overall Survival (months)')
ax.tick_params(axis='x', rotation=20)

# -- Right panel: 3-gene
ax = axes[1]

sns.violinplot(data=df,
               x='3-gene_classifier_subtype',
               y='overall_survival_months',
               order=three_gene_order,
               split=True,
               inner='quart',
               color=palette_dict["All Patients"],
               ax=ax
               )
sns.violinplot(data=df_death,
               x='3-gene_classifier_subtype',
               y='overall_survival_months',
               order=three_gene_order,
               split=True,
               inner='quart',
               color=palette_dict["Died of Disease"],
               ax=ax
               )

# Medians overlay
medians_all_3g   = df.groupby('3-gene_classifier_subtype')['overall_survival_months'].median()
medians_death_3g = df_death.groupby('3-gene_classifier_subtype')['overall_survival_months'].median()

for i, subtype in enumerate(three_gene_order):
    if subtype in medians_all_3g.index:
        y_val = medians_all_3g[subtype]
        # Line across the violin (blue = All)
        ax.hlines(y=y_val,
                  xmin=i-0.25,
                  xmax=i+0.25,
                  colors='blue',
                  linestyles='--',
                  linewidth=2)
        ax.text(i,
                y_val + 2,
                f'{y_val:.1f}',
                color='blue',
                ha='center',
                va='bottom',
                fontsize=9)
    if subtype in medians_death_3g.index:
        y_val = medians_death_3g[subtype]
        # Line across the violin (red = Death)
        ax.hlines(y=y_val,
                  xmin=i-0.25,
                  xmax=i+0.25,
                  colors='red',
                  linestyles='--',
                  linewidth=2)
        ax.text(i,
                y_val + 2,
                f'{y_val:.1f}',
                color='red',
                ha='center',
                va='bottom',
                fontsize=9)

# Formatting
ax.set_title('Overall Survival by 3-Gene Classification', pad=12)
ax.set_xlabel('3-Gene Subtype')
ax.tick_params(axis='x', rotation=20)
ax.tick_params(axis='y', labelleft=True)  # display Y tick numbers on right panel
ax.set_ylabel("")  # remove duplicate Y-axis label

# -- Global labels (general legend for both panels) ---
legend_elements = [Patch(facecolor=palette_dict['All Patients'],
                         edgecolor='none',
                         alpha=0.55,
                         label='All Patients'),
                   Patch(facecolor=palette_dict['Died of Disease'],
                         edgecolor='none',
                         alpha=0.65,
                         label='Died of Disease'),
                   Line2D([0], [0],
                          color='blue',
                          linestyle='--',
                          lw=2,
                          label='Median (All Patients)'),
                   Line2D([0], [0],
                          color='red',
                          linestyle='--',
                          lw=2,
                          label='Median (Died of Disease)')
                   ]

fig.legend(handles=legend_elements,
           loc='lower center',
           ncol=2,
           frameon=False,
           bbox_to_anchor=(0.5, 0.02),
           bbox_transform=fig.transFigure
           )

# -- Layout
plt.subplots_adjust(left=0.08, right=0.98,
                    top=0.85, bottom=0.25,
                    wspace=0.18
                    )

# Save figure
fig = plt.gcf()  # get the current active figure
save_figure(fig, '02_violin_survival',output_dir)

plt.show()


# ==== 3) Bar plots: Therapy proportions per classification ====
def cancer_classification_therapy(classification, therapies):
    """
    Calculate therapy percentages per molecular classification group.

    Parameters:
        classification (str): Column name for subtype classification.
        therapies (list[str]): Therapy column names (binary).

    Returns:
        pd.DataFrame: Percentages of patients receiving each therapy.
        - Returns a wide DataFrame with one row per classification level and one column per therapy.
        - Values are percentages in [0, 100].
    """
    df_therapy = df.copy()
    df_therapy[classification] = df_therapy[classification].fillna('Unknown')
    # Calculate mean (proportion) and convert to percent
    percent_by_therapy = df_therapy.groupby(classification)[therapies].mean()
    percent_by_therapy = percent_by_therapy * 100
    percent_by_therapy = percent_by_therapy.reset_index()
    return percent_by_therapy

fig, axes = plt.subplots(1, 2, figsize=(16, 8.5), sharey=True)

# -- Left: PAM50
ax = axes[0]
classification = 'pam50_+_claudin-low_subtype'
therapies = ['chemotherapy', 'hormone_therapy', 'radio_therapy']

percent_by_therapy_pam50 = cancer_classification_therapy(classification, therapies)
# transform to long format
df_therapy_plot_pam50 = percent_by_therapy_pam50.melt(id_vars=classification, var_name='Therapy', value_name='Percent')

sns.barplot(data=df_therapy_plot_pam50,
            x=classification,
            y='Percent',
            hue='Therapy',
            order=pam50_order,
            palette=[palette_dict['Chemotherapy'],
                     palette_dict['Hormone Therapy'],
                     palette_dict['Radiotherapy']],
            ax=ax
            )

# Formatting
ax.set_title('Therapy Use by PAM50 Classification', pad=12)
ax.set_xlabel('PAM50 Subtype', labelpad=17)
ax.set_ylabel('Patients Receiving Therapy (%)', labelpad=10)
ax.set_ylim(0, 100)
ax.tick_params(axis='x', rotation=20)

# Remove per-axes legend (will add a global one later)
leg = ax.get_legend()
if leg: leg.remove()

# -- Right: 3-gene
ax = axes[1]
classification = '3-gene_classifier_subtype'

percent_by_therapy_3g = cancer_classification_therapy(classification, therapies)
# transform to long format
df_therapy_plot_3g = percent_by_therapy_3g.melt(id_vars=classification, var_name='Therapy', value_name='Percent')

sns.barplot(data=df_therapy_plot_3g,
            x=classification,
            y='Percent',
            hue='Therapy',
            order=three_gene_order,
            palette=[palette_dict['Chemotherapy'],
                     palette_dict['Hormone Therapy'],
                     palette_dict['Radiotherapy']],
            ax=ax
            )

# Formatting
ax.set_title('Therapy Use by 3-Gene Classifier', pad=12)
ax.set_xlabel('3-Gene Subtype')
ax.tick_params(axis='x', rotation=20)
ax.tick_params(axis='y', labelleft=True)  # display Y tick numbers on right panel
ax.set_ylabel("")  # remove duplicate Y-axis label
ax.set_ylim(0, 100)

# Remove per-axes legend (will add a global one later)
leg = ax.get_legend()
if leg: leg.remove()

# -- Shared legend for both barplots
handles, labels = axes[1].get_legend_handles_labels()
fig.legend(handles,
           labels,
           ncol=3,
           frameon=False,
           loc='lower center',
           bbox_to_anchor=(0.5, 0.02))

# -- Layout
plt.subplots_adjust(left=0.07, right=0.99,
                    top=0.90, bottom=0.22,
                    wspace=0.20)

# Save figure
fig = plt.gcf()  # get the current active figure
save_figure(fig, '03_bar_therapy',output_dir)

plt.show()


# ==== 4) Heatmaps: Number of therapies per classification ====
def cancer_n_therapy(classification, therapies):
    """
    Calculate the percentage of patients receiving 0, 1, 2, or 3 therapies
    within each molecular classification group.

    This function sums the selected therapy columns per patient to count
    how many therapies each received. Then it groups the data by the chosen
    classification (e.g., PAM50 or 3-gene subtype) and calculates the
    percentage of patients in each category of therapy count.

    Parameters:
        classification (str): Column name for subtype classification.
        therapies (list[str]): Therapy column names (binary).

    Returns:
        pd.Series: A Series with a MultiIndex:
        (classification value, number of therapies) → percent of patients

    Notes:
    - Use `.unstack(fill_value=0)` on the output to get a DataFrame
      suitable for heatmap visualization.
    """
    df_n_therapies = df.copy()
    df_n_therapies['therapies count'] = df_n_therapies[therapies].sum(axis=1)

    proportion_by_n_therapy = df_n_therapies.groupby(classification)['therapies count'].value_counts(normalize=True)
    percent_by_n_therapy = proportion_by_n_therapy * 100
    return percent_by_n_therapy

fig, axes = plt.subplots(1, 2, figsize=(15, 9), sharey=False)

# -- Left: PAM50
ax = axes[0]
classification = 'pam50_+_claudin-low_subtype'
therapies = ['chemotherapy','hormone_therapy','radio_therapy']

percent_by_n_therapy_pam50 = cancer_n_therapy(classification,therapies)
# transform to wide format
df_n_therapies_plot_pam50 = percent_by_n_therapy_pam50.unstack(fill_value=0)

# Reorder rows according to PAM50 biological order
df_n_therapies_plot_pam50 = df_n_therapies_plot_pam50.reindex(index=pam50_order)

sns.heatmap(data=df_n_therapies_plot_pam50,
            annot=True,
            fmt='.1f',
            cmap='YlOrBr',
            cbar=False,
            ax=ax)

# Formatting
ax.set_title('Number of Therapies per Patient (PAM50)', pad=12)
ax.set_xlabel('Number of Therapies', labelpad=8)
ax.set_ylabel('PAM50 Subtype', labelpad=8)

# -- Right: 3-gene
ax = axes[1]
classification = '3-gene_classifier_subtype'

percent_by_n_therapy_3g = cancer_n_therapy(classification,therapies)
# transform to wide format
df_n_therapies_plot_3g = percent_by_n_therapy_3g.unstack(fill_value=0)

# Reorder rows according to 3-gene biological order
df_n_therapies_plot_3g = df_n_therapies_plot_3g.reindex(index=three_gene_order)

sns.heatmap(data=df_n_therapies_plot_3g,
            annot=True,
            fmt='.1f',
            cmap='YlOrBr',
            cbar=False,
            ax=ax)

# Formatting
ax.set_title('Number of Therapies per Patient (3-Gene Classifier)', pad=12)
ax.set_xlabel('Number of Therapies', labelpad=8)
ax.set_ylabel('3-Gene Subtype', labelpad=8)

# -- Shared colorbar (horizontal)
cbar = fig.colorbar(axes[1].collections[0],
                    ax=axes,
                    orientation='horizontal',
                    fraction=0.05,
                    pad=0.05,
                    label='Patients (%)')

# -- Layout
for ax in axes:
    ax.tick_params(axis='x', rotation=0)
    ax.tick_params(axis='y', labelrotation=0)
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=9)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=9)

plt.subplots_adjust(left=0.1, right=0.97,
                    top=0.85, bottom=0.25,
                    wspace=0.45
                    )
plt.suptitle('Distribution of Therapy Counts by Molecular Classification', y=0.95, fontsize=13)

# Save figure
fig = plt.gcf()  # get the current active figure
save_figure(fig, '04_heatmap_therapies',output_dir)

plt.show()

# ==== 5) Scatterplot: mutation_count vs tumor size ====
plt.figure(figsize=(9.2, 6.2))

sns.scatterplot(data=df,
                x='mutation_count',
                y='tumor_size',
                hue='tumor_stage',
                palette=palette_dict['Tumor Stage'],
                alpha=0.6,
                s=50)

# Labels and Formatting
plt.title("Mutation Load vs Tumor Size by Stage", pad=12)
plt.xlabel("Mutation Count", labelpad=8)
plt.ylabel("Tumor Size (mm)", labelpad=8)
plt.legend(title="Tumor Stage",
           frameon=False,
           bbox_to_anchor=(1, 1),
           )
plt.grid(alpha=0.2)
plt.tight_layout()

# Save figure
fig = plt.gcf()  # get the current active figure
save_figure(fig, '05_scatter_mutations',output_dir)

plt.show()
