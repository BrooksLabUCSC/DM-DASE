
######################################################################################
### QUANTIFY JUNCBASE AS EVENTS
######################################################################################

suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(readr)){install.packages("readr")})
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
source(".../plotting/barplot_main_program.R")


quantify_JuncBase_as_events = function(
                                      x,
                                      colors,
                                      output
                                      ) {
  

  options(dplyr.summarise.inform=FALSE)
  temp_global = dplyr::group_by(x, as_event_type) %>% 
    dplyr::summarise(count = n())
  temp_global$percent = round(temp_global$count / sum(temp_global$count)*100, 2)
  temp_global = temp_global %>% 
    dplyr::mutate(annotation = "Global") %>%
    dplyr::select(annotation, as_event_type, count, percent)

  options(dplyr.summarise.inform=FALSE)
  temp_ann = dplyr::group_by(x, 
      annotation = `#Contains_Novel_or_Only_Known(Annotated)_Junctions`,
      as_event_type) %>% 
    dplyr::summarise(count = n())

  temp_split = split(
    temp_ann, 
    f=temp_ann$annotation
    )

  temp_ann_events = lapply(temp_split, 
                           function (y) {
                             y$percent = round(y$count / sum(y$count)*100, 2)
                             y
                            }
                          )

  temp_ann = data.table::rbindlist(
                          temp_ann_events, 
                          use.names=T,
                          fill=T, 
                          idcol="cluster"
                          ) %>% 
  dplyr::select(-cluster)
  temp_ann$annotation[temp_ann$annotation == "K"] = "Known"
  temp_ann$annotation[temp_ann$annotation == "N"] = "Novel"
  complete_data = rbind(temp_global, temp_ann)
  
  data.table::fwrite(complete_data, 
                     paste0(output, "juncbase_as_events_proportions.tsv"), 
                     sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

  png(paste0(output, "juncbase_as_events_proportions.png"),
      width=9, height=6.5, units="in", res=600)
  print(custom_bar_plot(complete_data, colors))
  dev.off()

}
            

