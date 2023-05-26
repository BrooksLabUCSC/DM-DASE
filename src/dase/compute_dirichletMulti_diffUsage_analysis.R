
######################################################################################
### DRIMSEQ DU MAIN FUNCTION
######################################################################################

options(warn=-1)
dirichletMulti_diffUsage_analysis = function(x,
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
                                             delta_psi,
                                             speed, 
                                             threads,
                                             colors,
                                             output
                                             ) {
  
  if (input_format=="JB") {
    cat("\nQuantifying JuncBase AS events...\n")
    quantify_JuncBase_as_events(x, colors, output)
  }
  
  cat("\nChecking input data format...\n")
  x = as.data.frame(compute_data_format(x, input_format, select, event_types, output))
  x = filter_funct_as_events(x, gtf_data, input_format, release, output)
  y = as.data.frame(y)
  check_inputs(x, y)
  x = custom_names(x)

  cat("\nComputing global MT/WT analysis...\n")

  cat("\nCreating DrimSeq object...\n")
  dm_object = DRIMSeq::dmDSdata(counts=x, samples=y)
  
  summary_plot = plotData(dm_object) + custom_themes
  png(paste0(output, "transcripts_features_summary_plot.png"),
      width=5, height=5, units="in", res=600)
  print(summary_plot)
  dev.off()
  
  cat("\nFiltering low expressed transcripts...\n")
  dm_object_fil = dmFilter(dm_object, 
                           min_samps_gene_expr=as.numeric(samples_total),  
                           min_samps_feature_expr=as.numeric(samples_alt),
                           min_gene_expr=10,
                           min_feature_expr=5)
  
  #############################################################################
  ### Precision estimate using 10% (0.01) of randomly selected genes
  #############################################################################
  
  cat("\nComputing precision estimates...\n")
  design_full_mtx = stats::model.matrix(~group, 
                                        data=DRIMSeq::samples(dm_object_fil)
                                        )

  set.seed(123)
  dm.object.precision = dmPrecision(dm_object_fil, 
                                    design=design_full_mtx) 
  
  precision_plot = plotPrecision(dm.object.precision) + 
                      geom_point(size=4) + 
                      custom_themes
  png(paste0(output, "precision_estimates_plot.png"),
      width=5, height=5, units="in", res=600)
  print(precision_plot)
  dev.off()
  
  #############################################################################
  ### Proportion estimate
  #############################################################################
  
  cat("\nComputing proportion estimates...\n")
  dm.object.prop = dmFit(dm.object.precision, 
                          design=design_full_mtx, 
                          verbose=1) 

  #############################################################################
  ### Differential transcript usage test
  #############################################################################

  dm_fit_model = dmTest(dm.object.prop, 
                        coef=paste0("group", condition1), 
                        verbose=1)

  cat("\nComputing the null design matrix model...\n")
  design_null = model.matrix(~1, data=DRIMSeq::samples(dm.object.prop))

  cat("\nComputing Diritchlet multinomial test for each gene...\n")
  object_test = dmTest(dm_fit_model, design=design_null)
  object_test_res = DRIMSeq::results(object_test)
  dm_ann.gene = annotate_transcripts(object_test_res, gtf_data, 
                                     ann_data, input_format, release)

  data.table::fwrite(dm_ann.gene, 
                     paste0(output, "drimseq_gene_diffusage_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  cat("\nComputing Diritchlet multinomial test for each transcript...\n")
  dm_object_transcript = as.data.frame(DRIMSeq::results(object_test, level="feature"))
  
  object_test_res = object_test_res[order(object_test_res$pvalue, decreasing=FALSE), ]

  dm_object_transcript = dm_object_transcript[order(dm_object_transcript$pvalue, 
                                              decreasing=FALSE), ]
  dm_ann.transcript = annotate_transcripts(dm_object_transcript, gtf_data, 
                                           ann_data, input_format,release)

  data.table::fwrite(dm_ann.transcript, 
                     paste0(output, "drimseq_transcript_diffusage_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  cat("\nQuantify significant AS events proportions...\n")
  quantify_significant_as_events(dm_ann.gene, colors, output)

  pval_plot = plotPValues(object_test) + 
                custom_themes
  png(paste0(output, "drimseq_dtu_pval_hist.png"),
      width=5, height=5, units="in", res=600)
  print(pval_plot)
  dev.off()
  
  pval_feature_plot = plotPValues(object_test, level="feature") + 
                        custom_themes
  png(paste0(output, "drimseq_dtu_feature_pval_hist.png"),
      width=5, height=5, units="in", res=600)
  print(pval_feature_plot)
  dev.off()

  #############################################################################
  ### Two-stage test
  #############################################################################
  
  cat("\nComputing Two-stage test analysis...\n")
  pval_screen = DRIMSeq::results(object_test)$pvalue
  names(pval_screen) = as.character(DRIMSeq::results(object_test)$`gene_id`)
  
  pval_conf = matrix(DRIMSeq::results(object_test, level="feature")$pvalue, ncol=1)
  rownames(pval_conf) = DRIMSeq::results(object_test, level="feature")$feature_id
  
  transcript2gene = DRIMSeq::results(object_test, 
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
                     paste0(output, "stage_wise_adjpval_diff_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  #############################################################################
  ### Calculate mean proportions, LFC and Delta PSI values
  #############################################################################

  cat("\nCalculate mean proportions...\n")
  counts = as.data.frame(counts(dm_object))
  props = proportions(dm_fit_model)$feature_id 
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
  prop_fit = proportions(dm_fit_model)

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
  lfc_complete = merge(lfc_prop_temp, padj.df, 
                        by = c(intersect(names(lfc_prop_temp), names(padj.df)))
                      )
  lfc_complete = annotate_transcripts(lfc_complete, gtf_data, ann_data, 
                                      input_format, release)

  data.table::fwrite(lfc_complete, 
                     paste0(output, "drimseq_lfc_padj_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  cat("\nSaving Delta PSI results...\n") 
  deltapsi_complete = merge(mean_prop_data, padj.df, 
                              by = c(intersect(names(mean_prop_data), names(padj.df)))
                       	    )
  deltapsi_complete = annotate_transcripts(deltapsi_complete, gtf_data, 
                                           ann_data, input_format, release)
  
  data.table::fwrite(deltapsi_complete, 
                      paste0(output, "drimseq_dPSI_padj_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

}


