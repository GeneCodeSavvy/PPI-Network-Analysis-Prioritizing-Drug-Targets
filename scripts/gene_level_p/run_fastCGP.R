# Define command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) > 3) {
  cat("Usage: Rscript run_fastCGP.R <snp2gene_file> <snp_chr_pos_p_file> <genes2compute_file>\n")
  quit(status = 1)
}

# Assign command-line arguments to variables
snp2gene_file <- args[1]
snp_chr_pos_p_file <- args[2]
if (length(args) == 3) {
  genes2compute_file <- args[3]
}

# Load necessary functions
source("scripts/gene_level_p/fastCGP_function.R")


# Run the fastCGP function
output_direc <- dirname(snp2gene_file)
result <- fastCGP(snp2gene_file, snp_chr_pos_p_file, output_dir = output_direc)

# Print the result
cat("Computed gene p-values saved to:", output_direc, "\n")
