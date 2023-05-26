
######################################################################################
### DRIMSEQ DU COMPARISON MAIN FUNCTION
######################################################################################

dirichletMulti_diffUsage_comparison = function(x,
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
                                               output
                                               ) {

  cat("\nChecking input data format...\n")
  x = as.data.frame(compute_data_format(x, input_format, select, event_types, output))
  x = filter_funct_as_events(x, gtf_data, input_format, release, output)
  y = as.data.frame(y)
  check_inputs(x, y)
  x = custom_names(x)

  for (pos in 1:length(condition3)) {

    element = condition3[pos]
    mt_cond = paste0("mut_", element)
    wt_cond = paste0("wt_", element)
    conditions_list = c(mt_cond, wt_cond)

    sub_dir = paste0(mt_cond, "_", wt_cond, "/")
    update_dir = paste0(output, sub_dir)
    if (!dir.exists(update_dir)) { dir.create(update_dir) } 
    
    cat(paste0("\nComputing ", 
               toupper(mt_cond), "/",
               toupper(wt_cond), 
               " comparison analysis...\n")
    )
    cat("\n")
    
    meta_filtered = y %>% dplyr::filter(condition %in% c(conditions_list))
    samples_list = c(unique(meta_filtered$sample_id))
    counts = x %>% dplyr::select(feature_id, gene_id, c(samples_list))

    cat("\nCreating comparison DrimSeq object...\n")
    ds_data = DRIMSeq::dmDSdata(counts=counts, samples=meta_filtered)

    sum_plot = plotData(ds_data) + custom_themes
    png(paste0(output, sub_dir, "comparison_txs_features_summary_plot.png"),
      width=5, height=5, units="in", res=600)
    print(sum_plot)
    dev.off()

    cat("\nFiltering comparison low expressed transcripts...\n")
    nsamples = length(samples_list)
    samples_alt = length(samples_list)/2
    filtered = dmFilter(ds_data,
                        min_samps_gene_expr=as.numeric(nsamples), 
                        min_samps_feature_expr=as.numeric(samples_alt),
                        min_gene_expr=10,
                        min_feature_expr=5)

    #############################################################################
    ### Precision estimate using 10% (0.01) of randomly selected genes
    #############################################################################

    cat("\nComputing comparison precision estimates...\n")
    desing_mtx = stats::model.matrix(~condition, data=DRIMSeq::samples(filtered))
    
    set.seed(123)
    precision = dmPrecision(filtered, design=desing_mtx)
  
    precision_plot = plotPrecision(precision) + 
                        geom_point(size=4) + custom_themes
    png(paste0(output, sub_dir, "comparison_precision_estimates_plot.png"),
        width=5, height=5, units="in", res=600)
    print(precision_plot)
    dev.off()

    #############################################################################
    ### PROPORTION ESTIMATE
    #############################################################################
  
    cat("\nComputing proportion estimates...\n")
    fit_prop = dmFit(precision, design=desing_mtx, verbose=1) 

    #############################################################################
    ### DU TEST
    #############################################################################
  
    cat("\nFitting the Diritchlet multinomial model for comparison analysis...\n")
    
    test_model = dmTest(fit_prop, 
                        coef=paste0("condition", wt_cond), 
                        verbose=1)
  
    cat("\nComputing the null design matrix model...")
    design_null = model.matrix(~1, data=DRIMSeq::samples(fit_prop))

    cat("\nComputing Diritchlet multinomial test for each gene...\n")
    dm_test = dmTest(test_model, design=design_null)

    ds_res = DRIMSeq::results(dm_test)
    ann_gene = annotate_transcripts(ds_res, gtf_data, ann_data, 
                                    input_format, release)

    data.table::fwrite(ann_gene, 
                       paste0(output, 
                        paste0(mt_cond, "_", wt_cond, "_DS_gene_DU_results.tsv")
                        ), 
                      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    check_output(x, y, ds_res)
  
    cat("\nComputing Diritchlet multinomial test for each transcript...\n")
    object_tx = as.data.frame(DRIMSeq::results(dm_test, level="feature"))

    ds_res = ds_res[order(ds_res$pvalue, decreasing=FALSE), ]
    
    object_tx = object_tx[order(object_tx$pvalue, decreasing=FALSE), ]
    ann_tx = annotate_transcripts(object_tx, gtf_data, ann_data, 
                                  input_format, release)

    data.table::fwrite(ann_tx, 
                       paste0(output, 
                        paste0(mt_cond, "_", wt_cond, "_DS_tx_DU_results.tsv")
                        ), 
                      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
    pval_plot = plotPValues(dm_test) + custom_themes
    png(paste0(output, paste0(mt_cond, "_", wt_cond, "_DS_DU_pval_hist.png")),
        width=5, height=5, units="in", res=600)
    print(pval_plot)
    dev.off()

    pval_feature_plot = plotPValues(dm_test, level="feature") + custom_themes
    png(paste0(output, sub_dir, 
          paste0(mt_cond, "_", wt_cond, "_DS_DU_feature_pval_hist.png")),
        width=5, height=5, units="in", res=600)
    print(pval_feature_plot)
    dev.off()

    #############################################################################
    ### TWO-STAGE TEST
    #############################################################################

    cat("\nComputing two-stage test analysis...\n")

    pval_screen = DRIMSeq::results(dm_test)$pvalue
    names(pval_screen) = as.character(DRIMSeq::results(dm_test)$`gene_id`)
    pval_conf = matrix(DRIMSeq::results(dm_test, level="feature")$pvalue, ncol=1)
    rownames(pval_conf) = DRIMSeq::results(dm_test, level="feature")$feature_id
  
    transcript2gene = DRIMSeq::results(dm_test, 
                        level="feature")[, c("feature_id", "gene_id")]

    stagertx_object = stageRTx(pScreen=pval_screen,
                               pConfirmation=pval_conf,
                               pScreenAdjusted=FALSE, 
                               tx2gene=transcript2gene)
    stagertx_adj = stageR::stageWiseAdjustment(stagertx_object,
                                               method="dtu",
                                               alpha=0.05, 
                                               allowNA=TRUE)
    stagewise_padj = stageR::getAdjustedPValues(object=stagertx_adj, 
                                                order=FALSE, 
                                                onlySignificantGenes=FALSE)
    check_output2(counts, meta_filtered, stagewise_padj)
  
    data.table::fwrite(stagewise_padj, 
                        paste0(output, sub_dir, 
                          paste0(mt_cond, "_", wt_cond, "_stage_wise_adjpval_DU_results.tsv")
                        ), 
                        sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    #############################################################################
    ### COMPARISON (MT/WT) ANALYSIS (i.e. Calculates Log2FC from proportions)
    #############################################################################
  
    comp_counts = as.data.frame(counts(ds_data))
    props = proportions(test_model)$feature_id
    comp_prop = melt(comp_counts[comp_counts$feature_id %in% props, ], 
                      id = c("gene_id", "feature_id"))
    comp_prop = comp_prop[order(comp_prop$gene_id, 
                                comp_prop$variable, 
                                comp_prop$feature_id), ]
    comp_prop = comp_prop %>%
          group_by(gene_id, variable) %>%
          mutate(total = sum(value)) %>%
          group_by(variable, add=TRUE) %>%
          mutate(prop = value/total)

    cat("\nCalculating Delta PSI values...\n")

    prop_model = proportions(test_model)
    
    calculate_deltapsi = function(x, y) {

      wt_meta = y %>% dplyr::filter(group == condition1, 
                                    condition %in% c(mt_cond, wt_cond))
      mt_meta = y %>% dplyr::filter(group == condition2, 
                                    condition %in% c(mt_cond, wt_cond))
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
    mean_prop_data = calculate_deltapsi(prop_model, y)

    cat("\nCalculating LFC from the mean of PSI proportions...\n")
    log2fc_data = log2(mean_prop_data[condition2] / mean_prop_data[condition1])
    colnames(log2fc_data) = "log2fc"
    log2fc_data$feature_id = mean_prop_data$feature_id

    padj.df = as.data.frame(stagewise_padj) %>% 
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
                          paste0(output, sub_dir, 
                            paste0(mt_cond, "_", wt_cond, "_DS_LFC_padj_results.tsv")
                          ), 
                      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    cat("\nSaving Delta PSI results...\n") 
    deltapsi_complete = merge(mean_prop_data, padj.df, 
                                by = c(intersect(names(mean_prop_data), names(padj.df)))
                              )
    deltapsi_complete = annotate_transcripts(deltapsi_complete, gtf_data, 
                                              ann_data, input_format, release)

    data.table::fwrite(deltapsi_complete, 
                        paste0(output, sub_dir, 
                          paste0(mt_cond, "_", wt_cond, "_DS_dPSI_padj_results.tsv")
                        ),
                      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  }
}


