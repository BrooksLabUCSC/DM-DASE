
######################################################################################
### DRIMSEQ tuQTL MAIN FUNCTION
######################################################################################

dirichletMulti_diffUsage_tuQTL_analysis = function(x,
                                                   y,
                                                   gene_ranges,
                                                   genotype_data,
                                                   snp_ranges,
                                                   samples_total,
                                                   samples_alt,
                                                   speed,
                                                   threads,
                                                   output
                                                   ) {
  
  cat("\nChecking input data format for tuQTL analysis...")
  x = as.data.frame(x)
  y = as.data.frame(y) 
  check_inputs(x, y)
  
 stqtl_object =  dmSQTLdata(counts=x,
                            gene_ranges=gene_ranges,
                            genotypes=genotype_data,
                            snp_ranges=snp_ranges,
                            samples=y,
                            window=5e3)
  
  png(paste0(output, "drimseq_tuqtl_features_plot.png"),
     width=5, height=5, units="in", res=600)
  print(plotData(stqtl_object, plot_type="features") + custom_themes)
  dev.off()
  
  png(paste0(output, "drimseq_tuqtl_snps_plot.png"),
      width=5, height=5, units="in", res=600)
  print(plotData(stqtl_object, plot_type="snps") + custom_themes)
  de.off()

  png(paste0(output, "drimseq_tuqtl_blocks_plot.png"),
      width=5, height=5, units="in", res=600)
  print(plotData(stqtl_object, plot_type="blocks") + custom_themes)
  dev.off()
  
  stqtl_object_fil = dmFilter(stqtl_object,
                              min_samps_gene_expr=as.numeric(samples_total), 
                              min_samps_feature_expr=as.numeric(samples_alt),
                              minor_allele_freq=5,
                              min_gene_expr=10, 
                              min_feature_expr=10)
  
  set.seed(123)
  stqtl_object_precision = dmPrecision(stqtl_object_fil,
                                       BPPARAM=BiocParallel::MulticoreParam(threads))
  
  png(paste0(output, "drimseq_tuqtl_precision_hist.png"),
      width=5, height=5, units="in", res=600)
  print(plotPrecision(stqtl_object_precision) + custom_themes)
  dev.off()
  
  stqtl_object_prop = dmFit(stqtl_object_precision,
                            BPPARAM=BiocParallel::MulticoreParam(threads))
  
  # Method estimates gene-level null model proportions/coefficients and 
  # likelihoods and performs the likelihood ratio test
  stqtl_test = DRIMSeq::dmTest(stqtl_object_prop, 
                  BPPARAM=BiocParallel::MulticoreParam(threads))
  stqtl_test_res = DRIMSeq::results(stqtl_test)
  
  data.table::fwrite(stqtl_test_res, 
                     paste0(output, "drimseq_tuqtl_diff_results.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}

