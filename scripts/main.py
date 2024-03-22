import utils
from pathlib import Path
import numpy as np
import pandas as pd

def main(paper:list) :
    output_dir = paper[2]

    if (output_dir/"network_edges.csv").exists() and (output_dir/"network_nodes.csv").exists() and (output_dir/"computed_gene_p_values.tsv").exists():
        sigmod_edges, sigmod_nodes = utils.sigmod_input(edgefile=str(output_dir/"network_edges.csv"), 
                                     nodefile=str(output_dir/"network_nodes.csv"), 
                                     genefile=str(output_dir/"computed_gene_p_values.tsv"))
        edges = output_dir/"sigmod_edges.tsv"
        nodes = output_dir/"sigmod_nodes.tsv"
        sigmod_edges.to_csv(edges, index=False)
        sigmod_nodes.to_csv(nodes, index=False)
        utils.run_sigmod(sigmod_node=nodes, sigmod_edge=edges)
    
    else:
        df = pd.read_csv(paper[0], sep=paper[1])
        df = utils.rename_col_summary_stats(df)
        query = df["SNP"].tolist()
        genes = utils.batched_gprofiler(query, batch_size=30)
        genes = utils.extract_first_element(df= genes, columns=["ensgs", "gene_names"])
        df = utils.rsids2gene(df, genes)
        merged_df = utils.genes_with_proteins(protein, df)
        snp2gene = merged_df[["ensgs", "SNP"]]
        snp_chr_pos_p = df[["SNP", "CHR", "BP", "P"]]
        snp2gene_file = str(output_dir/"snp2gene.csv")
        snp_chr_pos_p_file = str(output_dir/"snp_chr_pos_p.csv")
        snp2gene.to_csv(snp2gene_file, index=False, header=["gene", "SNP"])
        snp_chr_pos_p.to_csv(snp_chr_pos_p_file, index=False, header=["SNP", "chr", "pos", "p"])
        utils.run_fastCGP(snp2gene_file=snp2gene_file, snp_chr_pos_p_file=snp_chr_pos_p_file)
    

if __name__ == "__main__":
    summary_stats_file : Path = Path("data/summary_stats")
    output_file : Path = Path("data/output")
    protein = pd.read_csv("data/gene_with_protein_products.tsv", sep='\t')

    summary_stats :list = [] # Paths of the summary stats file
    sepr :list = [] # Seperators used for the summary stats file
    output :list = [] # Paths of the folder of analyses of each file

    """ Creating a directory for summary statistic file to store analyses """
    for file_path in summary_stats_file.iterdir(): # Iterate over the contents of the summary stats directory
        if file_path.is_file(): # Check if the path is a file (not a directory)
            summary_stats.append(file_path)
            sepr.append(utils.infer_delimiter(file_path))
            folder_path = output_file/file_path.stem
            output.append(folder_path)
            if not folder_path.exists(): # Check if the directory is already created
                folder_path.mkdir() # Make directory for summary stat file

    paths :tuple = tuple(zip(summary_stats, sepr, output))

    results = list(map(main, paths))