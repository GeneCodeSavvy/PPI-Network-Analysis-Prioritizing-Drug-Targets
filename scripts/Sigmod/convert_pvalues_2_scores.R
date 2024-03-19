	##	this function used the bisection algorithm to find the optimal module
	
#	}

	p2score <- function(gene_ps, output) {
	  gene_ps$score = qnorm(1 - gene_ps$p) ## convert p-values to upper quantile of normal distribution
	  
	  neg_inf_index = which(is.infinite(gene_ps$score) & (0 > gene_ps$score)) ## p-values converted to -inf
	  pos_inf_index = which(is.infinite(gene_ps$score) & (0 < gene_ps$score)) ## p-values converted to +inf
	  
	  neg_inf_index_len = length(neg_inf_index)
	  pos_inf_index_len = length(pos_inf_index)
	  
	  if (neg_inf_index_len > 0) warning('Warning:\n The p-value of ', neg_inf_index_len, ' genes are converted to a score of negative infinity. We reset their scores as -9 for computation consideration.\n')
	  if (pos_inf_index_len > 0) warning('Warning:\n The p-value of ', pos_inf_index_len, ' genes are converted to a score of positive infinity. We reset their scores as 9 for computation consideration.\n')
	  
	  gene_ps$score[neg_inf_index] = -9
	  gene_ps$score[pos_inf_index] = +9
	  
	  if (length(neg_inf_index) > 0 || length(pos_inf_index) > 0) {
	    gene_ps$score = ifelse(is.infinite(gene_ps$score), ifelse(gene_ps$score < 0, -9, 9), gene_ps$score)
	    warning(paste0('Warning:\n The p-value of ', length(neg_inf_index) + length(pos_inf_index), ' genes are converted to -9 or 9 for computation consideration.\n'))
	  }
	  
	  filename = file.path(output, '/gene_scores.tsv')
	  write.table(gene_ps, file=filename, sep='\t', quote=FALSE, row.names=FALSE)
	  gene_scores = gene_ps
	  return(gene_scores)
	}
	