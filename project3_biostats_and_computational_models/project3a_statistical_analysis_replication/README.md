# Project 3a – Kaplan-Meier and Log-Rank Statistical Analysis Replication

## Description

This project replicates the Kaplan–Meier survival analysis of the **13-gene TEX-related T-cell exhaustion signature** using the **GSE14520 hepatocellular carcinoma (HCC)** dataset (platform GPL3921).

It integrates expression data, probe annotation, and clinical metadata to:

- Build gene-level expression matrices
- Compute a TEX signature–based risk score
- Stratify patients into High vs Low risk groups
- Perform survival analysis with Kaplan–Meier curves and log-rank testing

The pipeline reproduces the methodology reported in the original TEX signature publication.

📄 DOI: [10.21037/tcr-24-650](https://doi.org/10.21037/tcr-24-650)

---

## Features

- Automatic loading and cleaning of GEO expression matrix
- Parsing GPL3921 annotation to map probes → genes
- Extraction of the 13-gene TEX-related exhaustion signature
- Collapse of probe-level expression to gene-level values
- Integration of expression and clinical survival metadata
- Risk score computation using published Cox coefficients
- High/Low risk stratification (median split)
- Kaplan–Meier survival analysis with:
  - Confidence intervals
  - Log-rank test
  - Number-at-risk table
- Journal-style plots exported as PNG

---

## Project Files

- `tex_signature_survival_GSE14520.py`  
  Main script for expression processing, annotation mapping, risk scoring, and KM survival analysis.

- `GSE14520_GPL3921_series_matrix.txt.gz`  
  GEO expression matrix.

- `GPL3921.annot`  
  Probe-level annotation file for GPL3921 platform.

- `GSE14520_Extra_Supplement.txt.gz`  
  Supplementary clinical metadata with survival time and status.

- `km_GSE14520_replicated.png`  
  Output Kaplan–Meier survival figure.

---

## 13-Gene TEX-Related Signature

The analysis uses the TEX-related exhaustion signature reported in the paper.

Included genes: `HSPD1`, `UBB`, `DNAJB4`, `CALM1`, `LGALS3`,
`BATF`, `COMMD3`, `IL7R`, `FDPS`, `DRAP1`,
`RPS27L`, `PAPOLA`, `GPR171`

Each gene has an associated Cox regression coefficient, and the risk score is computed as:

$$
\text{Risk Score} = \sum_{i=1}^{13} (\text{Expression}_i \times \beta_i)
$$

---

## Requirements

- **Python 3.8+**
- Python libraries:
  - numpy
  - pandas
  - seaborn
  - matplotlib
  - lifelines

Install dependencies:

```bash
pip install numpy pandas seaborn matplotlib lifelines
```

## Usage

Run the full pipeline with:

```bash
python tex_signature_survival_GSE14520.py
```
The script will:

1. Load and clean the GEO series matrix  
2. Parse GPL3921 annotation and map probes to gene symbols  
3. Extract the 13-gene TEX signature  
4. Collapse probe-level to gene-level expression  
5. Integrate clinical survival metadata  
6. Compute TEX risk scores using Cox coefficients  
7. Split samples into High/Low risk groups (median split)  
8. Perform Kaplan–Meier survival analysis with log-rank test  
9. Save the final KM plot as `km_GSE14520_replicated.png`

---

## Output Files

### Figures

Stored in the working directory:

- **`km_GSE14520_replicated.png`**  
  Contains:
  - High-risk vs Low-risk survival curves  
  - 95% confidence intervals  
  - Number-at-risk table  
  - Log-rank p-value annotation  
  - Journal-style axes and layout  

### Data (During Execution)

- Gene-level expression matrix (samples × 13 genes)  
- Clinical survival dataframe (time, event)  
- TEX signature risk score per patient  
- High/Low risk group assignments  

---

## Notes

- GEO metadata lines beginning with `!` are skipped via `comment="!"`.  
- GPL3921 annotation file may start at different lines; the script **detects the correct header automatically**.  
- Multiple probes per gene are collapsed using the **mean**, which is standard for microarray summarization.  
- Survival status parsing accommodates multiple formats such as `dead`, `deceased`, or `"1"`.  
- Median risk score is used to define High vs Low risk groups, following the methodology in the referenced study.  
- Kaplan–Meier curves, CIs, and log-rank tests are performed using the **lifelines** package.  

---

## License

This project is distributed under the **MIT License**.

---

## Author

Juliana Patiño Gallego  
jpatinoga@unal.edu.co
