# PPI-Network-Analysis-Prioritizing-Drug-Targets

This Python script is designed for analyzing summary statistics data obtained from Genome Wide Association studies (GWAS). It offers functionalities for processing data, conducting gene-level analyses, and preparing input for further downstream analysis with SigMod. Below is a detailed description of the functionalities and usage instructions for the script.

## Table of Contents

- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Additional Notes](#additional-notes)

## Introduction

GWAS generate summary statistics files containing information about the genetic variants such as SNPs (Single Nucleotide Polymorphisms), their genomic locations, and associated statistical metrics like p-values. Analyzing such data can provide insights into the genetic basis of various traits and diseases.

GWAS explore the genetic causes ofcomplex diseases. However, classical approaches ignore the biological context of the genetic variants and genes under study. To address this shortcoming, one can use biological networks, which model functional relationships, tosearch for functionally related susceptibility loci. This project tries to analyse the GWAS results with SigMod.

SigMod is a novel and efficient method integrating GWAS results and gene network to identify a strongly interconnected gene module enriched in high association signals.

This Python script automates the analysis process by performing the following tasks:

1. **Data Preprocessing**: Functions for inferring delimiter types, standardizing column names, and extracting relevant information from input files.
2. **Gene Annotation**: Methods to annotate genetic variants with gene information and prepare data for downstream analysis.
3. **Protein Association**: Automation of filtering genes associated with proteins, linking genetic variants to protein-related data.
4. **FastCGP Computation**: Execution of FastCGP analysis, facilitating gene-level p-value computation for genetic variants.
5. **SIGMOD Preparation and Analysis**: Functions for preparing input data and executing Sigmod analysis, enabling the study of gene interactions and associations.


## Requirements

- Python 3.x
- Required Python libraries: `numpy`, `pandas`, `gprofiler`, `tqdm`

## Usage

1. **Prepare Input Data**: Ensure that your summary statistics files are stored in a directory (`data/summary_stats` by default).
2. **Install Dependencies**: Install the required Python libraries using `pip install -r requirements.txt`.
3. **Run the Script**: Execute the script `main.py`.
4. **Monitor Progress**: The script will provide updates on the progress of each analysis step.
5. **Review Output**: Check the output directory (`data/output` by default) for the processed data and analysis results.

## Input

The script expects summary statistics files in CSV, TSV, or TXT format. Each file should contain columns representing variant ID (`variant_id`), base pair location (`base_pair_location`), chromosome (`chromosome`), and p-value (`p_value`). Optionally, the file may include additional columns such as SNP ID (`rs_id`).

## Output

The script generates the following output:

- Processed data files:
  - `snp2gene.csv`: Mapping of SNPs to associated genes.
  - `snp_chr_pos_p.csv`: SNP information including chromosome, position, and p-value.
- Analysis results:
  - Gene-level p-values computed using FastCGP.
  - SIGMOD input files (`sigmod_edges.tsv` and `sigmod_nodes.tsv`) prepared for further analysis.

## Additional Notes

- Ensure that the necessary input files (`gene_with_protein_products.tsv`) are available in the specified locations.
- Customize file paths and settings in the script as per your requirements.
- For detailed information on the algorithm implementations and methodologies used, refer to the respective documentation of external tools like g:Profiler and SIGMOD.
