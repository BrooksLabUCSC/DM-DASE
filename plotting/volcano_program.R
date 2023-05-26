#!/usr/bin/env Rscript

######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

### PROGRAM USAGE

# Rscript volcano_program.R \
#  -i drimseq_lfc_padj_results.tsv \
#  -l "global" \
#  -c "MT/WT" \
#  -t "0.1" \
#  -f "1.5" \
#  -o /output/

######################################################################################

suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(tibble)){install.packages("tibble")})
suppressMessages(if(!require(ggrepel)){install.packages("ggrepel")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})

parser = ArgumentParser(description="volcano_program.R - A program to plot volcano from differential splicing analysis")
parser$add_argument("-i", "--input", type="character", required=TRUE, help="Path to Log2FC of proportions data (containing delta PSI values)")
parser$add_argument("-l", "--layout", type="character", required=TRUE, help="Specify analysis a layout") 
parser$add_argument("-c", "--comparison", type="character", required=TRUE, help="Specify conditions comparison of analysis") 
parser$add_argument("-t", "--threshold", required=TRUE, help="Specify p-adj threshold for gene labels") 
parser$add_argument("-f", "--fc_threshold", required=TRUE, help="Specify p-adj threshold for gene labels") 
parser$add_argument("-o", "--output", type="character", required=TRUE, help="Output directory for all DGE computations")

args = parser$parse_args()
input.file = args$input[1]
layout = args$layout[1] 
comparison = args$comparison[1]  
threshold = args$threshold[1]
fold_threshold = args$fc_threshold[1]  
output = args$output[1] 

input_data = data.table::fread(input.file)

options(warn=-1)
custom_volcano = function(x,  
                          threshold,
                          fold_threshold,
                          comparison,
                          layout,
                          value, 
                          description,
                          output, 
                          ...
                          ) {

    x = x[!is.infinite(x$log2fc), ] 
    volcano = ggplot(x, aes(x = log2fc, y = -log10(gene_padj), 
                        color = dplyr::case_when(gene_padj < 0.05 ~ "grey",
                                                 gene_padj < 0.1 ~ "#FC4E2A",
                                                 gene_padj >= 0.1 ~ "#EEB479"),
                        label = "symbol")) + 
        geom_point(size=2) +
        geom_vline(xintercept = c(-as.numeric(fold_threshold), 
                                  as.numeric(fold_threshold)
                                  ), 
                   col="blue", alpha=0.5, linetype="dashed", lwd=0.8, size=0.5) + 
        theme_bw(base_size=25, base_family="", base_line_size=2.0, base_rect_size=2.5) + 
        labs(x=paste0("Log2FC", " (", comparison, ")"), y="-log10(P-adj)") + 
        theme(title = element_text(size=20, vjust=0.05),
                plot.title = element_text(hjust=0.5),
                legend.title = element_text(size=20, color="black"), 
                legend.text = element_text(size=20, color="black"),
                axis.text = element_text(size=20, color="black"),
                axis.title.x = element_text(size=20, color="black"),
                axis.title.y = element_text(size=20, color="black"),
                strip.text = element_text(size=20, color="black"),
                legend.margin = margin(0, 0, 0, 0),
                legend.box.margin = margin(0.1, 0, 0, 0),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()) + 
        ggrepel::geom_text_repel(
                aes(label = ifelse(gene_padj < as.numeric(threshold) & log2fc <= -as.numeric(fold_threshold) | 
                                   gene_padj < as.numeric(threshold) & log2fc >= as.numeric(fold_threshold),
                                   symbol, "")
                ),
                size = 3,
                color = "black",
                min.segment.length = 0.3,
                box.padding=unit(0.1, "lines"),
                point.padding=unit(0.01, "lines"),
                force=TRUE,
                direction="both",
                max.overlaps=Inf) + 
        scale_color_manual(
                name = "Significance", 
                values = c("grey85", "#EEB479", "#FC4E2A"),
                label = c("Not Sig.", "P-adj < 0.1", "P-adj < 0.05")) + 
        guides(color = guide_legend(override.aes = list(size=3))) + 
        xlim(min(x$log2fc) - 2, max(x$log2fc) + 2) + 
        ylim(0, max(-log10(x$gene_padj)) + 2)

    return(volcano)
}

png(paste0(output, paste0(layout, "_volvano_plot.png")),
	width=7, height=5, units="in", res=600)
print(custom_volcano(input_data, threshold, fold_threshold, 
                     comparison, layout, value, description, 
                     output)
)
dev.off()


