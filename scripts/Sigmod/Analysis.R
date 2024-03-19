# Load necessary functions/scripts
source('scripts/Sigmod/construct_scored_network.R')
source('scripts/Sigmod/construct_string_net.R')
source('scripts/Sigmod/convert_pvalues_2_scores.R')
source('scripts/Sigmod/remove_isolated_nodes.R')
source('scripts/Sigmod/SigMod_core.R')
source('scripts/Sigmod/selection_path.R')
source('scripts/Sigmod/selection_path_fast.R')
source('scripts/Sigmod/SigMod_bisection.R')
source('scripts/Sigmod/SigMod_v00.R')
source('scripts/Sigmod/determine_selection_fast.R')

# Get command-line arguments for output directory, node list, and edge list
args <- commandArgs(trailingOnly = TRUE)
node_list <- args[1]
edge_list <- args[2]
output_direc <- file.path(dirname(node_list), 'SubNetwork')
if (!dir.exists(output_direc)) {
  dir.create(output_direc)
}

#node_list <- "data/output/AnsteeQM/sigmod_nodes.tsv"
#edge_list <- "data/output/AnsteeQM/sigmod_edges.tsv"

# Read gene information and convert p-values to scores
gene_ps = read.table(node_list, header=TRUE, sep='\t')
gene_scores = p2score(gene_ps, output = output_direc)

# Read network data and construct the scored network
network_data = read.table(edge_list, header=TRUE, sep='\t')
net = construct_scored_net(network_data, 
                           interaction_indices=c(1,2), 
                           weight_index=NULL, 
                           gene_scores, 
                           output = output_direc)

# Run the SigMod bisection algorithm on the network
res_info <- SigMod_bisection(net, lambda_max = 1, nmax = 300, maxjump = 30, output_dir = output_direc)

# Save the network and results
save(net, res_info, file = file.path(output_direc,'/sigmod.RData'))