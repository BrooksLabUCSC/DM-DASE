#!/usr/bin/env Rscript

######################################################################################
### GSEA heatmap
######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM USAGE

#Rscript .../gsea/gsea_heatmap_program.R \
#   -i .../global_analysis/mut_KRASG12V_wt_KRASG12V/gsea_krasG12V/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
#   -i2 .../global_analysis/mut_KRASG12V_wt_KRASG12V/gsea_krasG12V/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
#   -i3 .../global_analysis/mut_lacz_wt_lacz/gsea_lacz/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
#   -i4 .../global_analysis/mut_lacz_wt_lacz/gsea_lacz/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
#   -c "U2AF1_S34F_krasG12V" \
#   -c2 "U2AF1_WT_krasG12V" \
#   -c3 "U2AF1_S34F_lacZ" \
#   -c4 "U2AF1_WT_lacZ" \
#   -e global \
#   -s MSigDB_Hallmark_2020 \
#   -t 1 \
#   -pw 9 \
#   -ph 13 \
#   -o .../global_analysis/

######################################################################################

options(warn=-1)
suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(plyr)){install.packages("plyr")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(ggrepel)){install.packages("ggrepel")})
suppressMessages(if(!require(ggpubr)){install.packages("ggpubr")})
suppressMessages(if(!require(Seurat)){install.packages("Seurat")})
suppressMessages(if(!require(gtable)){install.packages("gtable")})
suppressMessages(if(!require(grid)){install.packages("grid")})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer")})

parser = ArgumentParser(description="gsea_heatmap_program.R - Computes heatmap from GSEA results of multiple conditions")
parser$add_argument("-i", "--input", type="character", required=TRUE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i2", "--input2", type="character", required=TRUE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i3", "--input3", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i4", "--input4", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i5", "--input5", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i6", "--input6", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i7", "--input7", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i8", "--input8", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-c", "--condition", type="character", required=TRUE, help="Specify condition for data input") 
parser$add_argument("-c2", "--condition2", type="character", required=TRUE, help="Specify condition for data input 2") 
parser$add_argument("-c3", "--condition3", type="character", required=FALSE, help="Specify condition for data input 3") 
parser$add_argument("-c4", "--condition4", type="character", required=FALSE, help="Specify condition for data input 4") 
parser$add_argument("-c5", "--condition5", type="character", required=FALSE, help="Specify condition for data input 5") 
parser$add_argument("-c6", "--condition6", type="character", required=FALSE, help="Specify condition for data input 6") 
parser$add_argument("-c7", "--condition7", type="character", required=FALSE, help="Specify condition for data input 7") 
parser$add_argument("-c8", "--condition8", type="character", required=FALSE, help="Specify condition for data input 8")
parser$add_argument("-t", "--threshold", required=TRUE, help="P-value threshold for visualization") 
parser$add_argument("-e", "--event", type="character", required=TRUE, help="Specify JuncBASE AS event")
parser$add_argument("-s", "--set", type="character", required=TRUE, help="Specify gene set from prerank analysis")
parser$add_argument("-pw", "--panel_width", required=TRUE, help="Integer or float describing figure panel width")
parser$add_argument("-ph", "--panel_height", required=TRUE, help="Integer or float describing figure panel height")
parser$add_argument("-o", "--output", type="character", required=TRUE,help="Path to output directory")   

args = parser$parse_args()
input.file = args$input[1]
input.file2 = args$input2[1]
input.file3 = args$input3[1]
input.file4 = args$input4[1]
input.file5 = args$input5[1]
input.file6 = args$input6[1]
input.file7 = args$input7[1]
input.file8 = args$input8[1]
condition = args$condition[1]
condition2 = args$condition2[1]
condition3 = args$condition3[1]
condition4 = args$condition4[1]
condition5 = args$condition5[1]
condition6 = args$condition6[1]
condition7 = args$condition7[1]
condition8 = args$condition8[1]
threshold = args$threshold[1]
panel_width = args$panel_width[1]
panel_height = args$panel_height[1]
set = args$set[1]
event = args$event[1]
output = args$output[1]

######################################################################################
### Intersect function 
######################################################################################

intersect_vecs = function(a, 
                          b, 
                          ...
                          ) {

  vectors = c(a, b, ...)
  count_data = length(list(a, b, ...))
  freq_dist = count(vectors)
  inters = freq_dist[which(freq_dist$freq==count_data), "x"]
  
  return(inters)
}

######################################################################################
### Read data and calculate z-scores
######################################################################################

if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
    !is.null(input.file3)==FALSE & !is.null(input.file4)==FALSE & 
    !is.null(input.file5)==FALSE & !is.null(input.file6)==FALSE & 
    !is.null(input.file7)==FALSE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data = rbind(temp.data, input.data2)
  
  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  inters = intersect(terms1, terms2)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==FALSE & 
           !is.null(input.file5)==FALSE & !is.null(input.file6)==FALSE & 
           !is.null(input.file7)==FALSE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data = rbind(temp.data, input.data2, input.data3)
  
  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  inters = intersect_vecs(terms1, terms2, terms3)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==TRUE &
           !is.null(input.file5)==FALSE & !is.null(input.file6)==FALSE & 
           !is.null(input.file7)==FALSE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data4 = data.table::fread(input.file4) %>% 
                  dplyr::mutate(condition=condition4)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4)

  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  terms4 = c(unique(input.data4$Term))
  inters = intersect_vecs(terms1, terms2, terms3, terms4)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==TRUE & 
           !is.null(input.file5)==TRUE & !is.null(input.file6)==FALSE & 
           !is.null(input.file7)==FALSE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data4 = data.table::fread(input.file4) %>% 
                  dplyr::mutate(condition=condition4)
  input.data5 = data.table::fread(input.file5) %>% 
                  dplyr::mutate(condition=condition5)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4, input.data5)

  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  terms4 = c(unique(input.data4$Term))
  terms5 = c(unique(input.data5$Term))
  inters = intersect_vecs(terms1, terms2, terms3, terms4, terms5)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==TRUE & 
           !is.null(input.file5)==TRUE & !is.null(input.file6)==TRUE & 
           !is.null(input.file7)==FALSE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data4 = data.table::fread(input.file4) %>% 
                  dplyr::mutate(condition=condition4)
  input.data5 = data.table::fread(input.file5) %>% 
                  dplyr::mutate(condition=condition5)
  input.data6 = data.table::fread(input.file6) %>% 
                  dplyr::mutate(condition=condition6)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4, input.data5, input.data6)

  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  terms4 = c(unique(input.data4$Term))
  terms5 = c(unique(input.data5$Term))
  terms6 = c(unique(input.data6$Term))
  inters = intersect_vecs(terms1, terms2, terms3, terms4, terms5, terms6)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==TRUE & 
           !is.null(input.file5)==TRUE & !is.null(input.file6)==TRUE & 
           !is.null(input.file7)==TRUE & !is.null(input.file8)==FALSE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data4 = data.table::fread(input.file4) %>% 
                  dplyr::mutate(condition=condition4)
  input.data5 = data.table::fread(input.file5) %>% 
                  dplyr::mutate(condition=condition5)
  input.data6 = data.table::fread(input.file6) %>% 
                  dplyr::mutate(condition=condition6)
  input.data7 = data.table::fread(input.file7) %>% 
                  dplyr::mutate(condition=condition6)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4, 
                     input.data5, input.data6, input.data7)

  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  terms4 = c(unique(input.data4$Term))
  terms5 = c(unique(input.data5$Term))
  terms6 = c(unique(input.data6$Term))
  terms7 = c(unique(input.data7$Term))
  inters = intersect_vecs(terms1, terms2, terms3, terms4, terms5, terms6, terms7)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

} else if (!is.null(input.file)==TRUE & !is.null(input.file2)==TRUE & 
           !is.null(input.file3)==TRUE & !is.null(input.file4)==TRUE & 
           !is.null(input.file5)==TRUE & !is.null(input.file6)==TRUE & 
           !is.null(input.file7)==TRUE & !is.null(input.file8)==TRUE) {

  temp.data = data.table::fread(input.file) %>% 
                dplyr::mutate(condition=condition)
  input.data2 = data.table::fread(input.file2) %>% 
                  dplyr::mutate(condition=condition2)
  input.data3 = data.table::fread(input.file3) %>% 
                  dplyr::mutate(condition=condition3)
  input.data4 = data.table::fread(input.file4) %>% 
                  dplyr::mutate(condition=condition4)
  input.data5 = data.table::fread(input.file5) %>% 
                  dplyr::mutate(condition=condition5)
  input.data6 = data.table::fread(input.file6) %>% 
                  dplyr::mutate(condition=condition6)
  input.data7 = data.table::fread(input.file7) %>% 
                  dplyr::mutate(condition=condition7)
  input.data8 = data.table::fread(input.file8) %>% 
                  dplyr::mutate(condition=condition8)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4, 
                     input.data5, input.data6, input.data7, input.data8)

  terms1 = c(unique(temp.data$Term))
  terms2 = c(unique(input.data2$Term))
  terms3 = c(unique(input.data3$Term))
  terms4 = c(unique(input.data4$Term))
  terms5 = c(unique(input.data5$Term))
  terms6 = c(unique(input.data6$Term))
  terms7 = c(unique(input.data7$Term))
  terms8 = c(unique(input.data8$Term))
  inters = intersect_vecs(terms1, terms2, terms3, terms4, 
                          terms5, terms6, terms7, terms8)

  input.data = input.data %>% 
    dplyr::filter(Term %in% c(inters)) %>%
    dplyr::mutate(z_score = ((NES-mean(NES, na.rm=TRUE)) / sd(NES, na.rm=TRUE))) %>%
    dplyr::mutate(zscore = ifelse(z_score <= -3, -3, z_score)) %>%
    dplyr::mutate(zscore = ifelse(zscore >= 3, 3, zscore))

}

######################################################################################
### Clean data 
######################################################################################

clean_data = function(x) {

	temp.data = as.data.frame(x) %>%
		dplyr::rename(pvalue=`FDR q-val`)

	options(warn=-1)
	temp.data = temp.data[!is.na(as.numeric(as.character(temp.data$pvalue))), ]
	temp.data$pvalue = as.numeric(temp.data$pvalue)
	temp.data = temp.data[!(temp.data$pvalue %in% c("pvalue")), ]
	
  sig.data = temp.data %>% 
              dplyr::filter(pvalue < as.numeric(threshold))
	
  data_sep = sig.data %>% 
              tidyr::separate(Term, into=c("Set", "term"), sep="__")
  
  return(data_sep)
}

######################################################################################
### Main function
######################################################################################

scaled_heatmap = function(x,
                          set,
                          panel_width,
                          panel_height,
                          output,
                          ...
                          ) { 

    x = clean_data(x)
    #term_sif = x %>% dplyr::filter(`NOM p-val`<0.05)
    #term_sif = unique(c(term_sif$term))

    temp = x %>% dplyr::filter(Set==set) #, term %in% c(term_sif))
    temp$condition = factor(temp$condition, 
    						      levels = c("U2AF1_WT_lacZ", 
                                 "U2AF1_S34F_lacZ", 
                                 "U2AF1_WT_krasWT",
    								             "U2AF1_S34F_krasWT", 
                                 "U2AF1_WT_krasG12V", 
                                 "U2AF1_S34F_krasG12V")
                      )
    temp_sig = temp %>% dplyr::filter(`NOM p-val`<0.05)
    scaled_colors = colorRampPalette(brewer.pal(11, "RdBu"))(25) %>% rev()

    plot = ggplot(temp, 
            aes(x=condition, y=reorder(term, -zscore), fill=zscore)
            ) + 
           geom_tile() + 
           geom_tile(data=temp, 
                     color=ifelse(temp$`NOM p-val`<0.05, "black", NA),
                     linewidth=1, alpha=0) +  
    			 scale_fill_gradientn(colours=scaled_colors,
                                limits=c(-3, 3),
                                label=c(">-3", "-2", "-1", "0", "1", "2", "<3"),
                         		 	  guide = guide_colorbar(label=TRUE, 
                                                       draw.ulim=TRUE, 
                                                       draw.llim=TRUE,
                                                       frame.colour="black", 
                                                       ticks=FALSE, 
                                                       nbin=10, 
                                                       frame.linewidth=2, 
                                                       label.position="right",
                                                       title="Z-score",
                                                       direction="vertical")) + 
          theme_minimal() +
          theme(axis.text.x = element_text(color="black", angle=45, 
                                           size=20, hjust=1),
          		  axis.text.y = element_text(color="black", size=20),
          		  legend.text = element_text(color="black", size=16),
          		  legend.title = element_text(color="black", size=20),
          		  axis.title.x = element_blank(), 
          		  axis.title.y = element_blank(),
          		  axis.ticks = element_blank(), 
                plot.margin = margin(t=15, r=15, b=35, l=150),
          		  panel.grid.major = element_blank(), 
          		  panel.grid.minor = element_blank()) +
    	    scale_y_discrete(position="right", limits=rev) 

    png(paste0(output, set, "_prerank_gsea_heatmap.png"), 
      width=as.numeric(panel_width), 
      height=as.numeric(panel_height), 
      units="in", res=600)
    print(plot) 
    dev.off()
}

scaled_heatmap(input.data, set, panel_width, panel_height, output)


