
######################################################################################
### Annotate Ensembl Gene IDs
######################################################################################

compute_gene_names = function(x,
			      y,
			      test,
			      ...
			      ) {

	temp.ann.data = y %>% dplyr::rename(ensembl_id=gene_id)

	if (test=="DTU") {
		
		temp_df = x[order(x$adj_pvalue), ]
		temp_sep = temp_df %>% 
			tidyr::separate(gene_id, into=c("ensembl"), sep="_", remove=FALSE) %>%
			tidyr::separate(ensembl, into=c("ensembl_id"), remove=FALSE)
		temp.dtu.data = dplyr::left_join(temp_sep, temp.ann.data, by="ensembl_id") %>% 
			rename(symbol=gene_name)
  
  		return(temp.dtu.data)

  	} else if (test=="LFC") {
		
		temp_df = x[order(x$gene_padj), ]
		temp_sep = temp_df %>% 
			tidyr::separate(gene_id, into=c("ensembl"), sep="_", remove=FALSE) %>%
			tidyr::separate(ensembl, into=c("ensembl_id"), remove=FALSE)
		temp.lfc.data = dplyr::left_join(temp_sep, temp.ann.data, by="ensembl_id") %>% 
			rename(symbol=gene_name)

		return(temp.lfc.data)
	}
}


