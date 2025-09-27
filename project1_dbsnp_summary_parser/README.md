# Project 1 – dbSNP Summary Parser

## Description
This project focuses on **reading and processing dbSNP summary files** to extract variant information (rsID, variant type, alleles, functional consequence and clinical significance) and convert them into a **tabular format** (TSV).

---

## Features

- **Extracts relevant fields**:
  - `rs_id`
  - `variant_type`
  - `allele_ref`
  - `allele_alt`
  - `functional_consequence`
  - `clinical_significance`
- **Tabular output** (`.tsv`) for downstream analysis.
- **Handles missing values**: marks clinical significance as `not-in-summary` when absent.
- **Beginner-friendly code**, fully commented for learning purposes.

---

## Requirements

- **Python 3.8** or higher  
- Library **pandas**  

Install dependencies:
```bash
pip install pandas
```

---

## Project Files

- `dbsnp_summary_parser.py`: Main script that performs the parsing.
- `examples/snp_result.txt`: Example input file (dbSNP summary text).
- `examples/dbsnp_summary.tsv`: Example output file generated from the example input.
- `README.md`: Documentation and tutorial.

---

## Input File

The script is designed to work with dbSNP summary files (.txt), which can be downloaded directly from NCBI dbSNP

To get your own input file:

1. Go to dbSNP and search for a gene or variant (e.g., BRCA1).
2. Use Send to → File → Summary to download the summary text file.
3. Save it in the same folder as the script (or update the filename in the code).

For convenience, this repository also includes an example file in the examples/ folder: snp_result.txt(./examples/snp_result.txt)


This example allows you to test the script immediately without downloading from NCBI.
