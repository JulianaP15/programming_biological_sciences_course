# Project 2 – METABRIC descriptive analysis

## Description

This project performs **exploratory analysis** on the METABRIC breast cancer dataset to relate **molecular subtype classifications** (PAM50 and 3-gene classifier) with **patient age, survival, and therapies.**
It generates a compact set of **publication-ready figures** (PDF + PNG) for quick inspection and reporting.

---

## Features

- **Age distribution** (histogram) for all patients vs. those who **died of disease**, including median overlays.
- **Overall survival** by **PAM50** and **3-gene** subtypes (violin plots) with **median annotations**.
- **Therapy usage** per subtype (bar plots): chemotherapy, hormone therapy, and radiotherapy.
- **Number of therapies** per subtype (heatmaps) showing % of patients with 0/1/2/3 therapies.
- **Mutation load vs tumor size** (scatter) colored by tumor stage.
- **Automatic figure saving** (`.pdf` vector + `.png` @ 300 DPI) via a helper `save_figure()`.
- **Beginner-friendly code**, fully commented for learning purposes.

---

## Requirements

- **Python 3.9** or higher  
- Libraries: **pandas**, **seaborn**, **matplotlib**

Install dependencies:
```bash
pip install pandas seaborn matplotlib
```

---

## Project Files

- `METABRIC_descriptive_analysis.py`: Main script with all analyses and figure exports.
- `METABRIC_RNA_Mutation.csv`: Input dataset (expected in the working directory).
- `figures/`: Output folder created automatically to store exported plots (PDF + PNG).
- `README.md`: Documentation and tutorial.

---

## Input File

The script expects a METABRIC-like CSV named: `METABRIC_RNA_Mutation.csv`
 
- Note: Certain columns are explicitly cast to object (by index) to avoid mixed-type parsing issues. If your CSV schema changes, update the indices passed to pd.read_csv(dtype=...)
---

## Usage

1. Run the script from a terminal:
```bash
python METABRIC_descriptive_analysis.py
```
3. Outputs (saved to `figures/`):
- `01_age_distribution.[pdf|png]`
- `02_violin_survival.[pdf|png]`
- `03_bar_therapy.[pdf|png]`
- `04_heatmap_therapies.[pdf|png]`
- `05_scatter_mutations.[pdf|png]`
  
You can change the output folder by editing:
```bash
output_dir = "figures"
```

---

## Code Explanation

Although everything is in a single script, each block is scoped and documented:

- **Global aesthetics & palettes**
  - Sets a light grid style and uses a consistent color scheme:
  - Histograms: `lightblue` vs `salmon`
  - Therapy bars: `seaborn Set2`
  - Heatmaps: `YlOrBr`
  - Scatter hue: `viridis`
- Subtype orders (`pam50_order`, `three_gene_order`) are defined and applied for consistent plots.
- `save_figure(fig, name, output_dir)` Saves each figure as vector PDF and 300 DPI PNG with tight bounding boxes.
- **Histogram (Age):** Two overlaid distributions: all patients vs. died of disease + median lines and labels.
- **Violin plots (Survival by subtype):** Two-panel figure (PAM50 and 3-gene). Each panel layers All patients and Died of disease and overlays per-subtype medians.
- `cancer_classification_therapy(classification, therapies)`: Computes % of patients receiving each therapy per subtype (wide DataFrame).
Used to build bar plots for PAM50 and 3-gene.
- `cancer_n_therapy(classification, therapies)`: Computes % distribution of the number of therapies (0/1/2/3) within each subtype (returns a Series with MultiIndex).
This is reshaped with `.unstack(fill_value=0)` to plot the heatmaps.
- **Scatter (Mutation load vs. tumor size):** Colored by tumor stage; legend placed outside to avoid overlap.

---

## Output Files

Each figure is exported twice: PDF (vector, ideal for print) and PNG (300 DPI, web-ready).
Filenames are prefixed with a numeric order for easy reference in reports:

- `01_age_distribution`
- `02_violin_survival`
- `03_bar_therapy`
- `04_heatmap_therapies`
- `05_scatter_mutations`

---

## Troubleshooting & Notes

- Seaborn warnings with `split=True`:
Some versions warn when `split=True` is used without `hue`. This script keeps it for visual layering since it works in the current environment.

- Category labels:
Ensure the values in your dataset match the strings in `pam50_order` and `three_gene_order` (case/spelling).
If your data use different labels, update those lists.

- CSV schema changes:
The `dtype={...}` argument in `pd.read_csv` uses **column indices**. If the supplier changes the CSV, adjust those indices to prevent parsing issues.

---

## License

This project is distributed under the [MIT License](../LICENSE).

---

## Author

Juliana Patiño Gallego
jpatinoga@unal.edu.co
