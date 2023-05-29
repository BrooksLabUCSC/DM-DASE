
######################################################################################
### Barplot main program
######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

### PROGRAM USAGE

######################################################################################

options(warn=-1)
suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(tidyverse)){install.packages("tiyverse")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(dplyr)){install.packages("plyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(ggrepel)){install.packages("ggrepel")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})

parser = ArgumentParser(description="barplot_gsea_program.R - computes barplot plot")
parser$add_argument("-i", "--input", type="character", required=TRUE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i2", "--input2", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i3", "--input3", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-i4", "--input4", type="character", required=FALSE, help="Input data from GSEAPY in the form of table/data-frame")
parser$add_argument("-e", "--event", type="character", required=TRUE, help="Specify JuncBASE AS event") 
parser$add_argument("-s", "--set", type="character", required=TRUE, help="Gene set from GSEA") 
parser$add_argument("-p", "--pval", type="character", required=TRUE, help="Specify p-value to use as threshold (i.e. Adjusted P-value, FDR q-val, etc.)") 
parser$add_argument("-t", "--threshold", required=TRUE, help="P-value threshold for visualization")
parser$add_argument("-height", "--height_param", required=TRUE, type="integer", help="Figure panel height")
parser$add_argument("-width", "--width_param", required=TRUE, type="integer", help="Figure panel width")
parser$add_argument("-m", "--method", type="character", required=TRUE, help="GSEA method used for analysis (i.e. prerank or enrichment)")
parser$add_argument("-g", "--group", type="logical", default=FALSE, required=FALSE, help="Specify group analysis if comparing multiple gene sets")
parser$add_argument("-o", "--output", type="character", required=TRUE, help="Path to output directory")   

args = parser$parse_args()
input.file = args$input[1]
input.file2 = args$input2[1]
input.file3 = args$input3[1]
input.file4 = args$input4[1]
event = args$event[1]
pval_name = args$pval[1]
set = args$set[1]
threshold = args$threshold[1]
method = args$method[1]
height = args$height_param[1]
width = args$width_param[1]
group = args$group[1]
output = args$output[1]


if (group==FALSE & 
    !is.null(input.file) & 
    is.null(input.file2) & 
    is.null(input.file3) & 
    is.null(input.file4)) {

  input.data = data.table::fread(input.file)

} else if (group==TRUE & 
           !is.null(input.file) & 
           !is.null(input.file2) & 
           is.null(input.file3) & 
           is.null(input.file4)) {

  temp.data = data.table::fread(input.file)
  input.data2 = data.table::fread(input.file2)
  input.data = rbind(temp.data, input.data2)

} else if (group==TRUE & 
           !is.null(input.file) & 
           !is.null(input.file2) & 
           !is.null(input.file3) & 
           is.null(input.file4)) {

  temp.data = data.table::fread(input.file)
  input.data2 = data.table::fread(input.file2)
  input.data3 = data.table::fread(input.file3)
  input.data = rbind(temp.data, input.data2, input.data3)

} else if (group==TRUE & 
           !is.null(input.file) & 
           !is.null(input.file2) &
           !is.null(input.file3) & 
           !is.null(input.file4)) {

  temp.data = data.table::fread(input.file)
  input.data2 = data.table::fread(input.file2)
  input.data3 = data.table::fread(input.file3)
  input.data4 = data.table::fread(input.file4)
  input.data = rbind(temp.data, input.data2, input.data3, input.data4)

}

######################################################################################
### Clean data 
######################################################################################

clean_data = function(x, 
                      pval_name, 
                      threshold,
                      method
                      ) {

	temp.data = x %>% dplyr::rename(pvalue=.data[[pval_name]]) 

  if (method=="prerank") {

    temp.data = temp.data %>% 
      tidyr::separate(Term, into=c("source", "term"), sep="__")

  } else if (method=="enrichment") {

    temp.data = temp.data %>% 
      tidyr::separate(Term, into=c("term", "feature"), sep=" \\(") %>%
      dplyr::select(-feature)

  }

	options(warn=-1)
	temp.data = temp.data[!is.na(as.numeric(as.character(temp.data$pvalue))), ]
	temp.data$pvalue = as.numeric(temp.data$pvalue)
	temp.data = temp.data[!(temp.data$pvalue %in% c("pvalue")), ]
  sig.data = temp.data %>% dplyr::filter(pvalue < as.numeric(threshold))
  
  if (nrow(sig.data)==0) { 
    stop("No significant terms found.")
  } else {
    return(sig.data)
  }
}

######################################################################################
### Main function
######################################################################################

compute_barplot = function(x,
						               event,
                           set,
                           pval_name,
                           group,
                           method,
                           height,
                           width,  
						               output
                           ) {

    temp_df = clean_data(x, pval_name, threshold, method)
    feature = ifelse(group==TRUE, "Gene_set", "term")
    
    if (group==TRUE) {
      feature_unique = c(unique(x[[feature]]))
      x$Gene_set = factor(x$Gene_set, levels=c(feature_unique))
    }

    plot = ggplot(temp_df, aes(x=-log10(pvalue), 
                               y=reorder(source, -log10(pvalue)), 
                               fill=.data[[feature]])) + 
              geom_bar(stat="identity", position=position_dodge(width=0.9), 
                    width=0.8, linewidth=1.5, color="black", 
                    fill="#ECCBAE") + 
              theme_classic(base_size=20, base_family="", 
                    base_line_size=1.5, base_rect_size=2.0) +
              theme(legend.title = element_text(color="black", size=20),
                    legend.position = ifelse(group==TRUE, "right", "none"),
                    legend.text = element_text(color="black", size=20),
                    axis.text = element_text(color="black", size=20),
                    axis.title.x = element_text(color="black", size=20),
                    axis.title.y = element_blank()) +
              labs(x=paste0(pval_name, " (-log10)")) +
              scale_fill_manual(values = c("#F3A583", "#3C93C2", "#6CBA7D", "#F2AD00")) 

    png(paste0(output, paste0(event, "_", set, "_terms_barplot.png")),
        width=width, height=height, units="in", res=600)
    print(plot)
    dev.off()

}

compute_barplot(input.data, event, set, pval_name, 
                group, method, height, width, output)


