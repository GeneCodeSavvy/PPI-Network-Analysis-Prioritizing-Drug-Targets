import utils
from pathlib import Path
import pandas as pd
from typing import List, Tuple

def main(paper: Tuple[Path, str, Path]) -> None:
    # Extracting output directory from the input list
    output_dir: Path = paper[2]

    # Checking if required files exist in the output directory
    if (output_dir/"network_edges.csv").exists() and \
       (output_dir/"network_nodes.csv").exists() and \
       (output_dir/"computed_gene_p_values.tsv").exists():
        
        # If all required files exist, process them for SIGMOD
        sigmod_edges, sigmod_nodes = utils.sigmod_input(edgefile=str(output_dir/"network_edges.csv"),
                                                        nodefile=str(output_dir/"network_nodes.csv"),
                                                        genefile=str(output_dir/"computed_gene_p_values.tsv"))
        
        # Saving processed SIGMOD edges and nodes
        edges: str = str(output_dir/"sigmod_edges.tsv")
        nodes: str = str(output_dir/"sigmod_nodes.tsv")
        sigmod_edges.to_csv(edges, index=False, sep="\t")
        sigmod_nodes.to_csv(nodes, index=False, sep="\t")
        
        # Running SIGMOD
        utils.run_sigmod(sigmod_node=nodes, sigmod_edge=edges)
    
    else:
        # If required files are missing, process summary statistics
        df: pd.DataFrame = pd.read_csv(paper[0], sep=paper[1])
        df = utils.rename_col_summary_stats(df)
        query: List[str] = df["SNP"].tolist()
        genes = utils.batched_gprofiler(query, batch_size=30)
        genes = utils.extract_first_element(df=genes, columns=["ensgs", "gene_names"])
        df = utils.rsids2gene(df, genes)
        merged_df = utils.genes_with_proteins(protein, df)
        snp2gene = merged_df[["ensgs", "SNP"]]
        snp_chr_pos_p = df[["SNP", "CHR", "BP", "P"]]
        snp2gene_file: str = str(output_dir/"snp2gene.csv")
        snp_chr_pos_p_file: str = str(output_dir/"snp_chr_pos_p.csv")
        
        # Saving processed data
        snp2gene.to_csv(snp2gene_file, index=False, header=["gene", "SNP"],)
        snp_chr_pos_p.to_csv(snp_chr_pos_p_file, index=False, header=["SNP", "chr", "pos", "p"])
        print("files saved")
        
         # Running fastCGP
        utils.run_fastCGP(snp2gene_file=snp2gene_file, snp_chr_pos_p_file=snp_chr_pos_p_file)

if __name__ == "__main__":
    # Paths for input and output files
    summary_stats_folder: Path = Path("data/summary_stats")
    output_file: Path = Path("data/output")
    protein: pd.DataFrame = pd.read_csv("data/gene_with_protein_products.tsv", sep='\t')

    summary_stats: List[Path] = []  # Paths of the summary stats files
    separators: List[str] = []      # Separators used for the summary stats files
    output: List[Path] = []          # Paths of the folder of analyses of each file

    # Creating a directory for each summary statistics file to store analyses
    for file_path in summary_stats_folder.iterdir():
        if file_path.is_file():
            summary_stats.append(file_path)
            separators.append(utils.infer_delimiter(file_path))
            folder_path: Path = output_file / file_path.stem
            output.append(folder_path)
            if not folder_path.exists():
                folder_path.mkdir()

    paths: List[Tuple[Path, str, Path]] = list(zip(summary_stats, separators, output))

    # Executing main function for each summary statistics file
    results: List[None] = list(map(main, paths))
