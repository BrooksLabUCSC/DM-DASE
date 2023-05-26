
######################################################################################
### DRIMSEQ BATCH CORRECTION DTU MAIN FUNCTION
######################################################################################

dirichletMulti_diffUsage_batch_analysis = function(x, 
                                                   y, 
                                                   ann_data,
                                                   gtf_data,
                                                   release,
                                                   input_format,
                                                   select,
                                                   event_types,
                                                   samples_total, 
                                                   samples_alt, 
                                                   condition1,
                                                   condition2,
                                                   condition3,
                                                   layout,
                                                   delta_psi,
                                                   speed,
                                                   threads,
                                                   output
                                                   ) {
  
  cat("\nChecking input data format for BATCH DTU analysis...")
  x = as.data.frame(compute_data_format(x, input_format, select, event_types, output))
  x = filter_funct_as_events(x, gtf_data, release, output)
  y = as.data.frame(y) %>% dplyr::rename("library_layout"=layout)
  check_inputs(x, y)
  
  data_object2 = DRIMSeq::dmDSdata(counts=x, samples=y)
  data_object2_fil = dmFilter(data_object2,
                              min_samps_gene_expr=as.numeric(samples_total), 
                              min_samps_feature_expr=as.numeric(samples_alt),
                              min_gene_expr=10,
                              min_feature_expr=0)
  
  design_full2 = model.matrix(~group+library_layout, 
                              data=DRIMSeq::samples(data_object2_fil)
                              )
  
  set.seed(123)
  precision_est_obj = dmPrecision(data_object2_fil, 
                                  design=design_full2)
  
  object_fit = dmFit(precision_est_obj,
                      design=design_full2,
                      verbose=1)
  
  cat("\nComputing Dirichlet multinomial test...")
  test_object = dmTest(object_fit, 
                       coef=paste0("group", condition1), 
                       verbose=1)
  
  test_object_res = DRIMSeq::results(test_object)
  batch_corr.dtu.gene = annotate_transcripts(test_object_res, gtf_data, 
                                             ann_data, release)

  data.table::fwrite(batch_corr.dtu.gene, 
                     paste0(output, "batch_corrected_drimseq_gene_diff_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  cat("\nComputing Dirichlet multinomial DU for transcripts...")
  test_object_txs = as.data.frame(DRIMSeq::results(test_object,
                                                   level="feature"))
  batch_corr.dtu.tx = annotate_transcripts(test_object_txs, gtf_data, 
                                           ann_data, release)
  
  data.table::fwrite(batch_corr.dtu.tx, 
                     paste0(output, "batch_corrected_drimseq_tx_diff_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


  #############################################################################
  ### Two-stage test
  #############################################################################
  
  cat("\nComputing Two-stage test analysis...\n")
  pval_screen = DRIMSeq::results(test_object)$pvalue
  names(pval_screen) = as.character(DRIMSeq::results(test_object)$`gene_id`)
  
  pval_conf = matrix(DRIMSeq::results(test_object, 
                                        level="feature")$pvalue, ncol=1)
  rownames(pval_conf) = DRIMSeq::results(test_object, 
                                          level="feature")$feature_id
  
  transcript2gene = DRIMSeq::results(test_object, 
                                     level="feature")[, c("feature_id", "gene_id")]

  stagertx_object = stageRTx(pScreen=pval_screen,
                             pConfirmation=pval_conf,
                             pScreenAdjusted=FALSE, 
                             tx2gene=transcript2gene)
  stagertx_adj_object = stageR::stageWiseAdjustment(stagertx_object,
                                                    method="dtu",
                                                    alpha=0.05, 
                                                    allowNA=TRUE)

  stagewise.padj = stageR::getAdjustedPValues(object=stagertx_adj_object, 
                                              order=FALSE, 
                                              onlySignificantGenes=FALSE)
  check_output2(x, y, stagewise.padj)
  
  data.table::fwrite(stagewise.padj, 
                     paste0(output, "batch_stage_wise_adjpval_diff_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  #############################################################################
  ### Calculate mean proportions, LFC and Delta PSI values
  #############################################################################
  
  cat("\nComputing LFC values from fitted proportions...\n")
  counts = as.data.frame(counts(data_object2))
  props = proportions(object_fit)$feature_id 
  temp_prop = melt(counts[counts$feature_id %in% props, ], 
                id = c("gene_id", "feature_id")) 
  temp_prop = temp_prop[order(temp_prop$gene_id, 
                              temp_prop$variable, 
                              temp_prop$feature_id), ]
  
  # Calculate proportions from counts
  temp_prop = temp_prop %>%
        group_by(gene_id, variable) %>%
        mutate(total = sum(value)) %>%
        group_by(variable, add=TRUE) %>%
        mutate(prop = value/total)

  # Calculate avg proportions from fitted mean of proportions
  prop_fit = proportions(object_fit)

  cat("\nCalculating Delta PSI values...\n")

  calculate_deltapsi = function(x, y) {

    wt_meta = y %>% dplyr::filter(group == condition1) 
    mt_meta = y %>% dplyr::filter(group == condition2) 
    deltapsi_df = mutate(x, 
                          MT = rowMeans(
                            select(x, c(unique(mt_meta$sample_id)))
                            )
                        )
    deltapsi_df = mutate(deltapsi_df, 
                          WT = rowMeans(
                            select(deltapsi_df, c(unique(wt_meta$sample_id)))
                            )
                        )
    
    deltapsi_df = deltapsi_df %>% 
      dplyr::mutate(deltaPSI = MT-WT) %>%
      dplyr::filter(deltaPSI < -delta_psi | deltaPSI > delta_psi)
    deltapsi_df = deltapsi_df[,!names(deltapsi_df) %in% c(y$sample_id)]

    return(deltapsi_df)
  }
  mean_prop_data = calculate_deltapsi(prop_fit, y)

  cat("\nCalculating LFC from the mean of PSI proportions...\n")
  log2fc_data = log2(mean_prop_data[condition2] / mean_prop_data[condition1])
  colnames(log2fc_data) = "log2fc"
  log2fc_data$feature_id = mean_prop_data$feature_id

  padj.df = as.data.frame(stagewise.padj) %>% 
                dplyr::select(feature_id = txID,
                              gene_id = geneID,
                              gene_padj = gene,
                              tx_padj = transcript)

  lfc_prop_data = merge(mean_prop_data, log2fc_data, by="feature_id")
  lfc_prop_temp = lfc_prop_data %>% 
                    dplyr::select(feature_id, gene_id, log2fc, -deltaPSI)
  batch_lfc_complete = merge(lfc_prop_temp, padj.df, 
                              by = c(intersect(names(lfc_prop_temp), names(padj.df)))
                            )
  batch_lfc_complete = annotate_transcripts(batch_lfc_complete, gtf_data, ann_data, 
                                            input_format, release)

  data.table::fwrite(batch_lfc_complete, 
                     paste0(output, "drimseq_batch_lfc_padj_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  cat("\nSaving Delta PSI results...\n") 
  batch_deltapsi = merge(mean_prop_data, padj.df, 
                          by = c(intersect(names(mean_prop_data), names(padj.df)))
                        )
  deltapsi_complete = annotate_transcripts(batch_deltapsi, gtf_data, 
                                           ann_data, input_format, release)
  
  data.table::fwrite(deltapsi_complete, 
                     paste0(output, "drimseq_batch_dPSI_padj_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

}


