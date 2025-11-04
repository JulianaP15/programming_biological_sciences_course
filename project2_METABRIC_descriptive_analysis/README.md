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

Although the script does not define explicit functions, each block handles a specific task:

- File reading
  - Opens the summary file line by line.
  - Identifies sections corresponding to each SNP.

- Field extraction
  - `rs_id`: detected in lines starting with rs.
  - `variant_type`, `allele_ref`, `allele_alt`: parsed from the following line using : and >.
  - `functional_consequence`: extracted from lines containing Functional Consequence:.
  - `clinical_significance`: usually found in the next line; marked as not-in-summary if absent.

- Table building
  - Creates a `DataFrame` with the extracted fields.
  - Facilitates exporting and analysis.

- Saving results
  - Exports results to a TSV file (tab-separated).

---

## Output Files

This script will create a file that contains:
- One row per variant.
- The following collumns: rs_id, variant_type, allele_ref, allele_alt, functional_consequence, and clinical_significance.

For demonstration, this repository includes an example output file generated from the example input: [dbsnp_summary.tsv](examples/dbsnp_summary.tsv)
``` bash
rs_id    variant_type    allele_ref    allele_alt    functional_consequence                     clinical_significance
rs12516  SNV             G             A,C,T         3_prime_UTR_variant,non_coding_transcript  benign
```
This file shows the expected structure and can be used to verify that the script runs correctly.

When analyzing your own data, update the filename in the script if necessary:
```bash
output_file = 'your_output_file.tsv'
```

---

## Error Handling

- If a field is missing, it is marked as not-in-summary.
- The script depends on the current dbSNP summary format; major changes may require parser updates.

---

## License

This project is distributed under the [MIT License](../LICENSE).

---

## Author

Juliana Patiño Gallego
jpatinoga@unal.edu.co