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

For convenience, this repository also includes an example file in the examples/ folder: [snp_result.txt](./examples/snp_result.txt)

This example allows you to test the script immediately without downloading from NCBI.

---

## Usage

1. If you want to test quickly, use the example file included in `examples/snp_result.txt`.
2. Run the [script](./dbsnp_summary_parser.py) in the terminal:
```bash
python dbsnp_summary_parser.py
```
3. The output will be saved as: `dbsnp_summary.tsv`
  
To analyze your own data, replace the example file with another dbSNP summary file and update the filename in the script if necessary:
```bash
dbsnp_summary_file = 'your_input_file.txt'
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
