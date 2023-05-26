
######################################################################################
### Quantify AS events significant porportions
######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

### Program filter for significant AS events in DrimSeq results using a the adjusted 
### p-value treshold. Then, it performs a proportion analysis across all AS events 
### inlcuding global, known and novel events, respectively. Program outputs a proportion
### table and a barplot summary of the analysis.
 
### PROGRAM USAGE

#Rscript quantify_sig_as_events_proportions.R \
#  -i drimseq_gene_diffusage_results.tsv \
#  -p 0.1 \
#  -e global \
#  -o /output/

######################################################################################

options(warn=-1)
suppressMessages(if(!require(argparse)){install.packages("argparse")})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})
source(".../plotting/barplot_main_program.R")

parser = ArgumentParser(description="quantify_sig_as_events_proportions.R - quantifies significant AS events porportions from DRIMSeq results.")
parser$add_argument("-i", "--input", type="character", required=TRUE, 
                      help="DRIMSeq output table containing the FWER p-value or p-adj and gene_id columns")
parser$add_argument("-p", "--pval_treshold", required=FALSE, 
                      help="Adjusted p-value threshold for AS events filtering") 
parser$add_argument("-e", "--events", required=FALSE, type="character", nargs='*', 
                      help="JuncBase event types (i.e., cassette, mutually_exclusive, etc.)")
parser$add_argument("-o", "--output", type="character", required=TRUE, help="Output directory")

args = parser$parse_args()
input.file = args$input[1]
event = args$events[1]
pval_treshold = args$pval_treshold[1]
output = args$output[1] 

if (!dir.exists(output)){ dir.create(output) } 
input.data = data.table::fread(input.file)

######################################################################################
### Clean data
######################################################################################

clean_drimseq_data = function(x, 
			      pval_treshold
                              ) {

  temp_sig_df = x %>% 
    dplyr::filter(adj_pvalue < as.numeric(pval_treshold))

  temp_sig_df$gene_id_ann = gsub(pattern="([0-9])_([a-zA-Z])", 
                                  replace="\\1;\\2", 
                                  as.character(temp_sig_df$gene_id)
                                  )
  temp_sig_df$as_event_ann = gsub(pattern=".*;", 
                                  replace="", 
                                  as.character(temp_sig_df$gene_id_ann)
                                  )
  temp_sig_df$as_event_type = gsub(pattern="_*[NK]", 
                                   replace="", 
                                   as.character(temp_sig_df$as_event_ann)
                                   )
  temp_sig_df$annotation = gsub(pattern=".*_", 
                                replace="", 
                                as.character(temp_sig_df$as_event_ann)
                                )
  temp_sig_df$as_event_type = gsub(pattern=" ", 
                                   replace="", 
                                   as.character(temp_sig_df$as_event_type)
                                   )
  temp_sig_df$annotation = gsub(pattern=" ", 
                                replace="", 
                                as.character(temp_sig_df$annotation)
                                )
  
  temp_df = temp_sig_df %>% 
    dplyr::select(gene_id, pvalue, adj_pvalue, 
                  symbol, as_event_type, annotation)
  
  return(temp_df)
}

######################################################################################
### Quantify significant AS events
######################################################################################

quantify_significant_as_events = function(x,
                                          pval_treshold,
                                          output
                                          ) {

  temp.clean_df = clean_drimseq_data(x, pval_treshold)

  options(dplyr.summarise.inform=FALSE)
  temp_global = dplyr::group_by(temp.clean_df, as_event_type) %>% 
    dplyr::summarise(count=n())
  temp_global$percent = round(temp_global$count / sum(temp_global$count)*100, 2)
  temp_global = temp_global %>% 
    dplyr::mutate(annotation="Global") %>%
    dplyr::select(annotation, as_event_type, count, percent)

  options(dplyr.summarise.inform=FALSE)
  temp_ann = dplyr::group_by(temp.clean_df, annotation, as_event_type) %>% 
    dplyr::summarise(count=n())
  temp_split = split(temp_ann, f=temp_ann$annotation)

  temp_ann_events = lapply(temp_split, 
                           function (y) {
                             y$percent = round(y$count / sum(y$count)*100, 2)
                             y
                             }
                           )

  temp_ann = data.table::rbindlist(temp_ann_events, 
                                   use.names=TRUE, 
                                   fill=TRUE, 
                                   idcol="cluster") %>%
             dplyr::select(-cluster)
  temp_ann$annotation[temp_ann$annotation == "K"] = "Known"
  temp_ann$annotation[temp_ann$annotation == "N"] = "Novel"
  complete_data = rbind(temp_global, temp_ann)
  
  data.table::fwrite(complete_data, 
                     paste0(output, "juncbase_sig_as_events_proportions.tsv"), 
                      sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  png(paste0(output, "juncbase_sig_as_events_proportions.png"),
      width=9, height=6.5, units="in", res=600) 
  print(custom_bar_plot(complete_data, colors))
  dev.off()

}
quantify_significant_as_events(input.data, pval_treshold, output)


