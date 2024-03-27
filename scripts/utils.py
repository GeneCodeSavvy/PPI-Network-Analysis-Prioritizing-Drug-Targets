from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from gprofiler import GProfiler
from tqdm.auto import tqdm
import ast

def infer_delimiter(file_path):
    """Infers the delimiter based on the file extension."""
    if file_path.suffix == ".csv":
        return ","
    elif file_path.suffix == ".tsv":
        return "\t"
    elif file_path.suffix == ".txt":
        return " "
    else:
        raise ValueError("Unsupported file extension. Please provide a CSV, TSV, or TXT file.")
    
def rename_col_summary_stats(summary_stats_file):
    # Standardize column names, select relevant columns
    summary_stats_file = summary_stats_file.rename(columns={
        'variant_id': 'SNP',
        'base_pair_location': 'BP',
        'chromosome': 'CHR',
        'p_value': 'P',
        'rs_id': 'SNP'
    }, errors='ignore')
    summary_stats_file = summary_stats_file[["SNP", "BP", "CHR", "P"]]
    return summary_stats_file

def batched_gprofiler(query, batch_size=2):
    gp = GProfiler(return_dataframe=True)
    results = []
    
    # Split query list into sublists
    sublist = [x.tolist() for x in np.array_split(query, batch_size)]
    
    # Initialize tqdm for progress tracking
    pbar = tqdm(total=len(sublist))
    
    for batch in sublist:
        df = gp.snpense(query=batch)
        results.append(df)
        pbar.update(1)  # Update progress bar
        
    pbar.close()  # Close progress bar
    
    return pd.concat(results, ignore_index=True)

def extract_first_element(df, columns):
    """
    Extracts the first element from lists in a specified column of a DataFrame.
    
    Parameters:
        df (DataFrame): The DataFrame containing the column with lists.
        column_name (str): The name of the column containing lists.
    
    Returns:
        DataFrame: A new DataFrame with the first element extracted from each list.
    """

    # Function to convert string representation of list to list
    def str_to_list(column):
        return column.apply(ast.literal_eval)
    
    # Define a function to extract the first element from a list
    def extract_first(lst):
        if isinstance(lst, list) and len(lst) > 0:
            return lst[0]
        else:
            return None
    
    # Apply the function to the specified column
    for col in columns:
        #df[col] = str_to_list(df[col])
        df[col] = df[col].apply(extract_first)
    
    return df

def rsids2gene(all_rsids, gene_assc_rsids):
    all_rsids = all_rsids.merge(gene_assc_rsids[['rs_id', 'ensgs', 'gene_names']], left_on='SNP', right_on='rs_id', how='left')
    all_rsids.drop("rs_id", axis=1, inplace=True)

    return all_rsids

def genes_with_proteins(genes_with_proteins, genes_associated_snps):
    #protiens = pd.read_csv('/Users/harshsharma/Desktop/[1] Project - NAFLD/[1] Summary Statistics/[1] Input/gene_with_protein_products.tsv', sep='\t')
    genes_with_proteins.rename(columns={"ensembl_gene_id":"ensgs"}, inplace=True)
    genes_with_proteins.dropna(axis=0, inplace=True, subset=["ensgs"])
    filtered_df = pd.merge(genes_associated_snps, genes_with_proteins[["ensgs"]], on="ensgs", how='inner')

    return filtered_df

def run_fastCGP(snp2gene_file, snp_chr_pos_p_file):
    # Define the Rscript command
    command = ["Rscript", "scripts/gene_level_p/run_fastCGP.R", snp2gene_file, snp_chr_pos_p_file]
    
    # Run the Rscript command
    try:
        subprocess.run(command, check=True)
        print("fastCGP computation completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error:", e)
        print("fastCGP computation failed.")

def sigmod_input(edgefile, nodefile, genefile):
    edges_df = pd.read_csv(edgefile)
    nodes_df = pd.read_csv(nodefile)
    gene_df = pd.read_csv(genefile, sep="\t")

    sigmod_nodes = gene_df.merge(nodes_df, left_on='gene', right_on='query term', how='inner')[["gene", "p"]]
    sigmod_edges = pd.DataFrame()
    
    # Extract gene names
    sigmod_edges['gene1'] = edges_df['name'].str.extract(r'(\d+\.\w+) \(pp\)')
    sigmod_edges['gene2'] = edges_df['name'].str.extract(r'(\d+\.\w+)$')

    # Merge with nodes_df to replace names with query terms
    sigmod_edges = sigmod_edges.merge(nodes_df[['name', 'query term']], left_on='gene1', right_on='name', how='left')
    sigmod_edges= sigmod_edges.merge(nodes_df[['name', 'query term']], left_on='gene2', right_on='name', how='left')

    sigmod_edges = sigmod_edges[["query term_x", "query term_y"]]

    # Rename columns (optional)
    sigmod_edges.rename(columns={'query term_x': 'gene1', 'query term_y': 'gene2'}, inplace=True)

    return sigmod_edges, sigmod_nodes

def run_sigmod(sigmod_node, sigmod_edge):
    # Define the Rscript
    command = ["Rscript", "scripts/Sigmod/Analysis.R", sigmod_node, sigmod_edge]

    # Run the Rscript command
    try:
        subprocess.run(command, check=True)
        print("Sigmod computation completed successfully.")
    except subprocess.CalledProcessError as e:
        print("Error:", e)
        print("Sigmod computation failed.")

