
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
source(".../plotting/barplot_main_program.R")

######################################################################################
### Clean data
######################################################################################

clean_sig_drimseq_data = function(x) {

  temp_sig_df = x %>% dplyr::filter(adj_pvalue < 0.1)

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
    dplyr::select(gene_id, pvalue, adj_pvalue, symbol, 
                  as_event_type, annotation)
  
  return(temp_df)
}

######################################################################################
### Quantify significant AS events proportions
######################################################################################

quantify_significant_as_events = function(x,
                                          colors,
                                          output
                                          ) {

  temp.clean_df = clean_sig_drimseq_data(x)

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
                             y$percent = round(y$count/sum(y$count)*100, 2)
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
#quantify_significant_as_events(input.data, pval_treshold, output)

