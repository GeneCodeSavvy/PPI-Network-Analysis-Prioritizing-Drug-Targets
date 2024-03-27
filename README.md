# PPI-Network-Analysis-Prioritizing-Drug-Targets

This Python script is designed for analyzing summary statistics data obtained from genetic studies. It offers functionalities for processing data, conducting gene-level analyses, and preparing input for further downstream analysis using SIGMOD (Significance Analysis of Gene MODules). Below is a detailed description of the functionalities and usage instructions for the script.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Requirements](#requirements)
- [Usage](#usage)
- [Input](#input)
- [Output](#output)
- [Additional Notes](#additional-notes)

## Introduction

Genetic studies often generate summary statistics files containing information about genetic variants such as SNPs (Single Nucleotide Polymorphisms), their genomic locations, and associated statistical metrics like p-values. Analyzing such data can provide insights into the genetic basis of various traits and diseases.

This Python script automates the analysis process by performing the following tasks:

1. **Data Preprocessing**: Standardizes column names, selects relevant columns, and prepares data in BED format.
2. **Gene Annotation**: Utilizes the g:Profiler tool to annotate SNPs with associated genes.
3. **Protein Association**: Identifies genes associated with protein products.
4. **FastCGP Computation**: Runs the FastCGP algorithm for gene-level p-value computation.
5. **SIGMOD Preparation and Analysis**: Prepares input data for SIGMOD analysis and executes the SIGMOD algorithm for identifying significant gene modules.

## Features

- Supports CSV, TSV, and TXT file formats for summary statistics data.
- Batch processing for efficient annotation of large datasets.
- Automated handling of missing or incomplete data.
- Integration with external tools like g:Profiler and SIGMOD.
- Progress tracking using tqdm for better user experience.

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
