#!/usr/bin/env Rscript

######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

### PROGRAM USAGE

### Gene-level
#Rscript qqplot_main_program.R \
#   -i drimseq_gene_diffusage_results.tsv \
#   -o /output/ \
#   -t 10e-8 \
#   -l "gene" \
#   -e "global"

### Transcript-level
#Rscript qqplot_main_program.R \
#   -i drimseq_transcript_diffusage_results.tsv \
#   -o /output/ \
#   -t 10e-8 \
#   -l "txs" \
#   -e "global"

######################################################################################
### Required Libraries and Argument Parser
######################################################################################

options(warn=-1)
suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(tidyverse)){install.packages("tiyverse")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(ggrepel)){install.packages("ggrepel")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})

parser = ArgumentParser(description="Compute QQ plot figure")
parser$add_argument("-i", "--input", type="character", required=TRUE, help="Input data requiring identifier and p-value columns in the form of table/data-frame") 
parser$add_argument("-e", "--event", type="character", required=TRUE, help="Specify a label or event")
parser$add_argument("-l", "--level", type="character", required=TRUE, help="Specify gene or transcript level")
parser$add_argument("-t", "--threshold", required=TRUE, help="P-value threshold for plotting")
parser$add_argument("-o", "--output", type="character", required=TRUE, help="Path to output directory")

args = parser$parse_args()
input.file = args$input[1]
level = args$level[1]
event = args$event[1]
threshold = args$threshold[1]
output.dir = args$output[1]     

######################################################################################
### Clean data
######################################################################################

clean_data = function(x) {

	temp.data = data.table::fread(x) %>%
		dplyr::rename(pvalue=`pvalue`)

	options(warn=-1)
	temp.data = temp.data[!is.na(as.numeric(as.character(temp.data$pvalue))), ]
	temp.data$pvalue = as.numeric(temp.data$pvalue)

	temp.data = temp.data[!(temp.data$pvalue %in% c("pvalue")), ]
	
  return(temp.data)
}

######################################################################################
### Calculate expected p-values
######################################################################################

calculate_expvalues_function = function(x) {

	temp.sorted = x[order(x$pvalue), ]
	raw_p = c(as.numeric(x$pvalue))
	counts = length(!is.na(raw_p))
	exp_count = 1:counts / (counts + 1)
	temp.sorted$pval_exp = exp_count 

	custom.df = temp.sorted[1:20, ]
	custom.ids = c(unique(custom.df$symbol))
	
  temp.sorted = temp.sorted %>%
			dplyr::mutate(custom_label = 
				ifelse(
					symbol %in% c(custom.ids) & as.numeric(temp.sorted$pvalue)<5e-8, 
          TRUE, FALSE
					)
				)

	temp.sorted
}

######################################################################################
### QQ Plot main function
######################################################################################

compute_qqplot_function = function(x, 
                                   output, 
                                   threshold, 
                                   level, 
                                   event
                                   ) {

  cat("\nCleaning data input data...\n")
  x = clean_data(x) 
  
  cat("\nCalculating expected p-values...\n")
  x = calculate_expvalues_function(x)

  options(warn = -1)  
  qq_plot = ggplot() +
    geom_point(aes(x = -log10(as.numeric(x$pval_exp)), 
                   y = -log10(as.numeric(x$pvalue))),
               color = dplyr::case_when(
               						as.numeric(x$pvalue) < 5e-8 ~ "#0D4A70", 
               						as.numeric(x$pvalue) < 5e-6 ~ "#00CD6C",  
               						as.numeric(x$pvalue) >= 5e-6 ~ "lightgrey"  
               						), 
               size=1.5) + #, alpha=0.9
    geom_abline(col="dimgrey", linetype="dashed", lwd=0.8, size=1.2) +
    theme_bw(base_size=20, base_family="", base_line_size=2.0, base_rect_size=2.5) +
    theme(legend.title = element_text(color="black", size=20), 
          legend.text = element_text(size=20),
          axis.text = element_text(color="black", size=20), 
          axis.title.x = element_text(size=20), 
          axis.title.y = element_text(size=20),
          title = element_text(size=20, vjust=0.1),
          plot.title = element_text(hjust=0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "Expected p-value (-log10)", 
         y = "Observed p-value (-log10)") +
    ggrepel::geom_text_repel(aes(-log10(as.numeric(x$pval_exp)), 
                                 -log10(as.numeric(x$pvalue))
                                 ),
                             label = ifelse(
                             			x$custom_label==TRUE & as.numeric(x$pvalue)<as.numeric(threshold), 
                                  x$symbol, 
                                  ""
                                  ),
                             color = ifelse(as.numeric(x$pvalue)<0.01, "#0D4A70", "black"),
                             size = 3,
                             min.segment.length = 0.3,
                             box.padding = unit(0.60, "lines"),
                             point.padding = unit(0.30, "lines"),
                             segment.color = ifelse(as.numeric(x$pvalue)<0.01, 
                                                    "#0D4A70", "lightgrey"
                                                    ),
                             max.overlaps = Inf) + 
    scale_fill_manual(name = "Significance", 
                      values = c("lightgrey", "#9CCEA7", "#0D4A70") ,
                      labels = c("Not Sig.", "P-val < 5e-6", "P-val < 5e-8")) + 
    guides(color = guide_legend(verride.aes=list(size=3))) +
    ylim(0, max(-log10(x$pvalue)) + 0.5) + 
    xlim(0, max(-log10(x$pval_exp)) + 0.5) 
  
  cat("\nWriting figure...\n")
  png(paste0(output, paste0(event, "_", level, "_qqplot.png")),
    width=5, height=5, units="in", res=600)
  return(qq_plot)
  dev.off()
}

cat("\nComputing QQ plot...\n")
compute_qqplot_function(input.file, output.dir, threshold, level, event)


