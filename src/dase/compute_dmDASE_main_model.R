#!/usr/bin/env Rscript

######################################################################################

# Author: Carlos Arevalo
# Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

# Program computes differential usage (DU) and transcript usage QTL (tuQTL) 
# analyses in RNA-seq data using DrimSeq, computing a Dirichlet Multinomial (DM) distribution 
# across all transcripts of a given gene. The model inputs either a JuncBase output 
# alternative splicing (AS) exclusion-inclusion counts table 
# (JBoutput_AS_exclusion_inclusion_counts_lenNorm.txt) or a DrimSeq formatted input table. 
# If input data is JB format, program will compute DrimSeq format and write an output file.
# Program outputs two DrimSeq DTU results, one for gene differential usage and another for
# transcript differential usage, in the form of txt files. When computing a Two-stage test 
# and/or batch correction, program will output txt files for each analysis, respectively.

### TEST 

# Program was tested on a R virtual environment and the main R libraries were installed using
# Conda prior to running the program:

# Create a conda R-virtual env
# conda create -n "temp_rEnv" r-essentials r-base
# conda activate temp_rEnv
# conda install -c bioconda bioconductor-drimseq
# conda install -c bioconda bioconductor-stager
# conda install -c conda-forge r-reticulate
# conda install -c r r-gridextra
# conda install -c conda-forge argparse

# Package can also be installed directly in R using BiocManager:

# if (!require("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("devtools")
# BiocManager::install("DRIMSeq")
# BiocManager::install("stageR")
# BiocManager::install("reticulate")
# BiocManager::install("ggpubr")
# BiocManager::install("gridExtra")

### PROGRAM USAGE

#Rscript compute_dmDASE_main_model.R \
#  -i jb_table_AS_exclusion_inclusion_counts_lenNorm.txt \
#  -m juncbase_metadata.tsv \
#  -a ensembl_annotation.csv \
#  -g gencode.v38.annotation.gtf \
#  -f "JB" \
#  -s TRUE # (if performing AS event-specific analysis, default == FALSE) \
#  -e "intron_retention" \ # (optional: if s == TRUE)
#  -s1 18 \
#  -s2 9 \
#  -c TRUE \
#  -c1 "MT" \
#  -c2 "WT" \
#  -c3 "lacz" "krasWT" "krasG12V" \ # (i.e. for wt_lacz, mut_lacz, wt_krasWT, mut_krasWT, etc)
#  -t "DU" \
#  -n 8 \
#  -b FALSE \
#  -l "condition" \
#  -o .../output/

######################################################################################
### Required Libraries and Argument Parser
######################################################################################

options(warn=-1)
suppressMessages(if(!require(DRIMSeq)){install.packages("DRIMSeq")})
suppressMessages(if(!require(stageR)){install.packages("stageR")})
suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(xtable)){install.packages("xtable")})
suppressMessages(if(!require(grid)){install.packages("grid")})
suppressMessages(if(!require(gridExtra)){install.packages("gridExtra")})
suppressMessages(if(!require(purrr)){install.packages("purrr")})
suppressMessages(if(!require(stats)){install.packages("stats")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})
suppressMessages(if(!require(reshape2)){install.packages("reshape2")})
source(".../src/dase/compute_dirichletMulti_diffUsage_analysis.R")
source(".../src/dase/compute_batch_dirichletMulti_diffUsage_analysis.R")
source(".../src/dase/compute_dirichletMulti_diffUsage_comparison.R")
source(".../src/dase/compute_dirichletMulti_diffUsage_QTL_analysis.R")
source(".../src/dase/filter_functional_as_events.R")
source(".../src/dase/quantify_JuncBase_as_events.R")
source(".../src/dase/quantify_sig_as_events_proportions.R")
source(".../src/dase/compute_drimseq_format.R")
source(".../src/dase/annotate_ensembl_ids.R")
source(".../src/dase/infer_unannotated_events.R")
source(".../src/dase/barplot_main_program.R")
source(".../src/dase/utils.R")

parser = ArgumentParser(description="compute_dmDASE_main_model.R - Differential usage (DU) and transcript usage QTL (tuQTL) analysis with DRIMSeq")
parser$add_argument("-i", "--input", type="character", required=TRUE, help="JuncBase AS exclusion-inclusion counts or MESA junction counts table in DRIMSeq format") 
parser$add_argument("-i2", "--input2", type="character", help="Gene ranges file for tuQTL analysis in the form of table/data-frame") 
parser$add_argument("-i3", "--input3", type="character", help="SNPs genotypes raw counts/dosages for tuQTL analysis in the form of table/data-frame") 
parser$add_argument("-i4", "--input4", type="character", help="SNPs ranges for tuQTL analysis in the form of table/data-frame") 
parser$add_argument("-a", "--annotation", type="character", required=TRUE, help="Ensembl gene ID annotations data in the form of table/data-frame")
parser$add_argument("-r", "--release", type="character", required=FALSE, help="Reference release date (i.e., ##date: 2019-12-13 for gencode.v33.chr_patch_hapl_scaff.annotation.gtf")
parser$add_argument("-m", "--metadata", type="character", required=TRUE, help="Metadata input in form of table/data-frame")
parser$add_argument("-f", "--format", type="character", required=TRUE, default="DS", help="Specify input data format as either JB, MESA or DS")
parser$add_argument("-s", "--specific", type="logical", required=FALSE, default="FALSE", help="Define AS event-specific analysis")
parser$add_argument("-e", "--events", type="character", required=FALSE, nargs='*', default="intron_retention", help="JuncBase event types (i.e., cassette, mutually_exclusive, etc.)")
parser$add_argument("-s1", "--samples1", type="integer", required=TRUE, help="Total number of samples for filtering")
parser$add_argument("-s2", "--samples2", type="integer", required=TRUE, help="Total number of mutant/knock-down samples for filtering")
parser$add_argument("-c", "--compare", type="logical", required=FALSE, default="TRUE", help="Specify comparisons analyses across sample conditions")
parser$add_argument("-c1", "--condition1", required=TRUE, help="Specify condition one corresponding to samples for DU analysis to compare against c2")
parser$add_argument("-c2", "--condition2", required=FALSE, help="Specify condition two corresponding to samples for DU analysis to compare against c1")
parser$add_argument("-c3", "--condition3", type="character", required=TRUE, nargs='+', help="List of treatments/conditions")
parser$add_argument("-t", "--test", type="character", required=TRUE, default="DU", help="Specify DRIMSeq test to perform (i.e. DU, DUC or tuQTL method)")
parser$add_argument("-n", "--threads", type="integer", required=FALSE, default=8, help="Specify number of threads/multicores to use")
parser$add_argument("-b", "--batch", type="logical", default="FALSE", help="Specify batch correction for DRIMSeq DU analysis")
parser$add_argument("-l", "--layout", type="character", help="If batch correction is apply, specify library_layout column")
parser$add_argument("-g", "--gtf", type="character", required=TRUE, help="Genome reference in the form of .gtf file")
parser$add_argument("-psi", "--psi_thresh", type="character", required=FALSE, default="0.1", help="Delta PSI cutoff filtering (i.e., deltaPSI > 10%, deltaPSI > 20%, etc.)")
parser$add_argument("-z", "--speed", type="character", required=FALSE, default="FALSE", help="Speeds processs by calculating one gene block for all blocks in a gene")
parser$add_argument("-o", "--output_dir", type="character", required=TRUE, help="DRIMSeq DU, DUC or tuQTL output directory")

args = parser$parse_args()
input.file = args$input[1]
generange.file = args$input2[1]
genotype.file = args$input3[1]
snprange.file = args$input4[1]
ann.file = args$annotation[1]
meta.file = args$metadata[1]  
gtf.file = args$gtf[1] 
release = args$release[1]             
format = args$format[1] 
select = args$specific[1]          
event_types = args$events[1]
samples_total = args$samples1[1]
samples_alt = args$samples2[1]
condition1 = args$condition1[1]
condition2 = args$condition2[1]
condition3 = c(args$condition3)
delta_psi = as.numeric(args$psi_thresh)
method = args$test[1]
compare = args$compare[1]
threads = args$threads[1]
batch = args$batch[1]
layout = args$layout[1]
speed = args$speed[1]
output = args$output_dir[1] 

input.data = data.table::fread(input.file)
metadata = data.table::fread(meta.file)
metadata = metadata[order(metadata$group), ]
ann.data = data.table::fread(ann.file)
gtf.data = data.table::fread(gtf.file)
colors = unique(palette_sampling(palette, 15))
if (!dir.exists(output)){ dir.create(output) } 

######################################################################################
### CHECK INPUTS MAIN
######################################################################################

check_inputs = function(x, y) { 
  
  cat("\nChecking input data...\n")
  count_temp = x %>% dplyr::select(-feature_id, -gene_id)
  count_samples = c(unique(colnames(count_temp)))
  meta_samples = c(unique(y$sample_id))
  
  if (is.null(dim(x))) {
    stop("DRIMSeq data has no dimensions!")
  } else if (is.null(dim(y))) {
    stop("Metadata has no dimensions!")
  } else if (is.null(colnames(x))) {
    stop("Column names of the counts table should be given.")
  } else if (is.null(colnames(y))){
    stop("Column names of the metadata table should be given.")
  } else if (all(count_samples %in% meta_samples)==FALSE) { 
    stop("Count samples should be the same as metadata samples.") 
  } else {
    cat("\nConditions met! Computing DRIMSeq analysis...\n")
  }
}

######################################################################################
### CHECK OUTPUTS MAIN 
######################################################################################

### Compares input data with model output results
### Params: Input counts, metadata, DRIMSeq output
check_output = function(x, y, z) { 
  
  counts_fids = c(unique(x$feature_id))
  counts_gids = c(unique(x$gene_id))
  meta_fids = c(unique(x$feature_id))
  meta_gids = c(unique(x$gene_id))
  output_fids = c(unique(x$feature_id))
  output_gids = c(unique(x$gene_id))
  
  if (is.null(dim(x))) {
    stop("DRIMSeq data has no dimensions!")
  } else if (is.null(dim(y))){
    stop("Metadata has no dimensions!")
  } else if (is.null(dim(z))){
    stop("Output data has no dimensions!")
  } else if (all(output_fids %in% counts_fids)==FALSE) { 
    stop("Output transcripts should be the same as input count transcripts!")
  } else if (all(output_gids %in% counts_gids)==FALSE) { 
      stop("Output genes should be the same as input count genes!")
  } else if (all(output_fids %in% meta_fids)==FALSE) { 
    stop("Output transcripts should be the same as metadata transcripts!")
  } else if (all(output_gids %in% meta_gids)==FALSE) { 
    stop("Output genes should be the same as metadata genes!")
  } else {
    cat("\nDRIMSeq analysis computed sucessfully!\n")
  }
}

### Check two-stage test results
check_output2 = function(x, y, z) { 
  
  counts_gids = c(unique(x$gene_id))
  meta_gids = c(unique(x$gene_id))
  output_gids = c(unique(x$geneID))
  
  if (is.null(dim(x))) {
    stop("DRIMSeq data has no dimensions!")
  } else if (is.null(dim(y))){
    stop("Metadata has no dimensions!")
  } else if (is.null(dim(z))){
    stop("Output data has no dimensions!")
  } else if (all(output_gids %in% counts_gids)==FALSE) { 
    stop("Two-stage genes should be the same as input count genes!")
  } else if (all(output_gids %in% meta_gids)==FALSE) { 
    stop("Two-stage genes should be the same as metadata genes!")
  } else {
    cat("\nTwo-stage analysis computed sucessfully!\n")
  }
}

######################################################################################
### DM-DASE MAIN FUNCTION
######################################################################################

if (method=="DU" & compare==FALSE & batch==FALSE) {

  options(warn=-1)
  dirichletMulti_diffUsage_analysis(input.data, metadata, ann.data, gtf.data, release, 
                                    format, select, event_types, samples_total, 
                                    samples_alt, condition1, condition2, condition3, 
                                    delta_psi, speed, threads, colors, output)
  
  cat("\nDRIMSeq DU analysis complete.\n")

} else if (method=="DUC" & compare==TRUE & batch==FALSE) {

  options(warn=-1)
  dirichletMulti_diffUsage_comparison(input.data, metadata, ann.data, gtf.data, release,
                                      format, select, event_types, samples_total, 
                                      samples_alt, condition1, condition2, condition3, 
                                      delta_psi, speed, threads, output)
  
  cat("\nDRIMSeq DU comparison analysis complete.\n")

} else if (method=="DU" & compare==TRUE & batch==FALSE) {

  options(warn=-1)

  dirichletMulti_diffUsage_analysis(input.data, metadata, ann.data, gtf.data, release, 
                                    format, select, event_types, samples_total, 
                                    samples_alt, condition1, condition2, condition3, 
                                    delta_psi, speed, threads, colors, output)
  
  dirichletMulti_diffUsage_comparison(input.data, metadata, ann.data, gtf.data, release, 
                                      format, select, event_types, samples_total, 
                                      samples_alt, condition1, condition2, condition3, 
                                      delta_psi, speed, threads, output)
  
  cat("\nDRIMSeq DU analyses complete.\n")

} else if (method=="DU" & compare==FALSE & batch==TRUE) {

  options(warn=-1)
  dirichletMulti_diffUsage_analysis(input.data, metadata, ann.data, gtf.data, release, 
                                    format, select, event_types, samples_total, 
                                    samples_alt, condition1, condition2, condition3, 
                                    delta_psi, speed, threads, colors, output)
  
  dirichletMulti_diffUsage_batch_analysis(input.data, metadata, ann.data, gtf.data, 
                                          release, format, select, event_types, 
                                          samples_total, samples_alt, condition1, 
                                          condition2, condition3, layout, delta_psi, 
                                          speed, threads, output)
  
  cat("\nDRIMSeq DU and batch corrected analyses complete.\n") 

} else if (method=="DU" & compare==TRUE & batch==TRUE) {

  options(warn=-1)
  dirichletMulti_diffUsage_analysis(input.data, metadata, ann.data, gtf.data, release,
                                    format, select, event_types, samples_total, 
                                    samples_alt, condition1, condition2, condition3, 
                                    delta_psi, speed, threads, colors, output)
  
  dirichletMulti_diffUsage_comparison(input.data, metadata, ann.data, gtf.data, release, 
                                      format, select, event_types, samples_total, 
                                      samples_alt, condition1, condition2, condition3, 
                                      layout, delta_psi, speed, threads, output)
  
  dirichletMulti_diffUsage_batch_analysis(input.data, metadata, ann.data, gtf.data, 
                                          release, format, select, event_types, 
                                          samples_total, samples_alt, condition1, 
                                          condition2, condition3, layout, delta_psi, 
                                          speed, threads, output)

  cat("\nDRIMSeq DU and batch corrected analyses complete.\n")

} else if (method=="tuQTL") { 

  options(warn=-1)
  dirichletMulti_diffUsage_tuQTL_analysis(input.data, metadata, gene_data, 
                                          genotype_data, snp_data, samples_total, 
                                          samples_alt, speed, threads, output) 
  
  cat("\nDRIMSeq tuQTL analysis complete.\n")

} else {

  cat("\nNeeds to specify DRIMSeq method (i.e., DU, DUC or tuQTL)!\n")

}


