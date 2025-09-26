"""
Script: dbsnp_summary_parser.py
Author: Juliana Patiño Gallego
Description:
    Simple parser for NCBI dbSNP "summary" .txt files.
    It extracts a subset of fields and saves them into a tabular format (TSV).

Fields extracted:
    - rs_id
    - variant_type
    - allele_ref
    - allele_alt
    - functional_consequence
    - clinical_significance

Notes / assumptions:
    * This script is designed for summary text files downloaded from dbSNP.
    * It expects that lines follow the same layout as seen in the example.
    * Example of such a file: 'snp_result.txt', as downloaded directly
      from the dbSNP web interface
    * Parsing relies on string splitting, so format changes may break it.
    * Clinical significance will be marked as 'not-in-summary' if not present.
"""

# ==== Library imports ====
import pandas as pd

# ==== Configuration ====
# Input file: dbSNP summary in plain text format
dbsnp_summary_file = 'snp_result.txt'

# Lists to store extracted values for each SNP
rs_id = []
variant_type = []
allele_ref = []
allele_alt = []
functional_consequence = []
clinical_significance = []

# Open the input file and process line by line
with open(dbsnp_summary_file) as file:
    for line in file:
        # ---- rsID ----
        if 'rs' in line:
            temp = line.split(' ')
            rs_id.append(temp[1])  # take the second element after split
            if len(temp) > 5:  # sometimes rsID is found further in the line when it has merged
                rs_id[-1] = temp[5]

            # Next line usually contains "variant_type:reference_allele>alternative_allele(s)"

            # ---- variant type ----
            new_line = next(file).strip()
            temp = new_line.split(':')
            variant_type.append(temp[0])

            # ---- reference and alternative alleles ----
            temp2 = temp[1].split('>')
            allele_ref.append(temp2[0])
            allele_alt.append(temp2[1])

        # ---- functional consequence ----
        if 'Functional Consequence: ' in line:
            new_line = line.replace('Functional Consequence: ', '')
            functional_consequence.append(new_line)

            # ---- clinical significance ----
            # Clinical significance line usually follows right after functional consequence
            new_line = next(file).strip()
            if 'Clinical significance: ' in new_line:
                temp = new_line.replace('Clinical significance: ', '')
            else:
                temp = 'not-in-summary'  # mark as missing
            clinical_significance.append(temp)

# Build DataFrame with extracted data
df = pd.DataFrame({
    'rs_id': rs_id,
    'variant_type': variant_type,
    'allele_ref': allele_ref,
    'allele_alt': allele_alt,
    'functional_consequence': functional_consequence,
    'clinical_significance': clinical_significance,
})

# Save DataFrame to TSV file (tab-separated)
df.to_csv('dbsnp_summary.tsv', sep='\t', index=False)



# Print quick length check for all lists
print(len(rs_id),
      len(variant_type),
      len(allele_ref),
      len(allele_alt),
      len(functional_consequence),
      len(clinical_significance)
      )



