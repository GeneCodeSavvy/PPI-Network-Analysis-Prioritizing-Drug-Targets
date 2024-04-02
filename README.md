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
- Cytoscape Software
- Required Cytoscape Plugins: `CytoHubba`, `Stringdb`

## Input

The script expects summary statistics files in CSV, TSV, or TXT format. Each file should contain columns representing variant ID (`variant_id`), base pair location (`base_pair_location`), chromosome (`chromosome`), and p-value (`p_value`). Optionally, the file may include additional columns.
The script also expects `gene_with_protein_products.tsv` in the `\data\` folder.

## Usage

### Part I- Processing GWAS summary statistics for analyses.

1. **Prepare Input Data**: Ensure that your summary statistics files are stored in the directory `data/summary_stats/` .
2. **Run the Script**: Execute the script `main.py`.
3. **Monitor Progress**: The script will provide updates on the progress of each analysis step.
4. **Review Output**: Check the output directory for the processed data and analysis results.

### Part II- Analysis of the network
1. **Prepare Input Data**: Use `computed_genes.tsv` to generate a Protein Protein Interaction Network in Cytoscape using Stringdb plugin
2. **Save Network files**: Save the edge table and node table as `network_edges.csv` and `network_nodes.csv` respectively in the `data/summary_stat/*`
3. **Run the Script**: Execute the script `main.py`
4. **Review Output**: Check the output directory for the results.


## Output

### Part I - Processing GWAS summary statistics for analyses
The script generates the following output:
A folder with the name of the summary statistic file is created
- Processed data files in the new folder:
  - `snp2gene.csv`: Mapping of SNPs to associated genes.
  - `snp_chr_pos_p.csv`: SNP information including chromosome, position, and p-value.
  - `computed_gene_p_value.tsv`: Genes and their respective calculated p-value
  - `computed_genes.tsv`: List of Genes from the `computed_gene_p_value.tsv` file

### Part II - Analysis of the network 
- Analysis results:
  - Gene-level p-values computed using FastCGP.
  - SIGMOD input files (`sigmod_edges.tsv` and `sigmod_nodes.tsv`) prepared for further analysis.

## Additional Notes

- Customize file paths and settings in the script as per your requirements.
- For detailed information on the algorithm implementations and methodologies used, refer to the respective documentation of external tools like g:Profiler, SigMod and fastCGP.
