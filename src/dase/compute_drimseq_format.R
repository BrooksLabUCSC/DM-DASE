
######################################################################################
### FORMAT DATA FUNCTION
######################################################################################

### Check input data format. If data is in JuncBase table format, program will 
### converted to DrimSeq format. If the data format it is not specify, program 
### will stop and call an error

juncbase2drimseq_format = function(x) {
      
      temp_collapse = apply(x, 1, 

        function(x) {
        
          collapse_rows = c(x["gene_name"], 
                            x["exclusion_junctions"],
                            x["inclusion_junctions"], 
                            x["as_event_type"], 
                            x["junc_ann"])
          collapsed_data = paste(collapse_rows, collapse=";")
        
          return(collapsed_data)
        }
      )

      temp_df2 = x %>% 
        dplyr::mutate(gene_id = c(temp_collapse), 
                      feature = exclusion_junctions) %>%
        dplyr::select(-junc_ann, -as_event_type, -gene_name, -exclusion_junctions, 
                      -inclusion_junctions, -chr, -strand, -exclusion_exons, 
                      -inclusion_exons, -neighboring_constitutive_exons, 
                      -`intron-exon_junctions`)

      temp_df3 = temp_df2 %>% dplyr::select(-feature, -gene_id)
      sample_list = c(colnames(temp_df3))
      
      compute_as_events = function(x, 
                                   sample_list
                                   ) {
        
        new.list = list()
        for (pos in 1:length(sample_list)) {
          
          sample = sample_list[pos]
          event1 = paste0(sample, "_inclusions")
          event2 = paste0(sample, "_exclusions")
          
          events_df = x %>% 
            dplyr::select(feature, gene_id, sample) %>%
            tidyr::separate(sample, into=c(event1, event2), sep=";") %>%
            tidyr::gather(key="condition", value="counts", 3:4) %>%
            dplyr::mutate(feature_id = case_when(
              condition==event1 ~ paste0(feature, "_inclusion"),
              condition==event2 ~ paste0(feature, "_exclusion"))) %>%
            dplyr::select(feature_id, gene_id, counts) 
          
          events_df$counts = as.numeric(events_df$counts)  
          events_df = events_df %>% dplyr::rename(!!sample:=counts)
          new.list[[length(new.list) + 1]] = events_df
        }

        return(new.list)
      }
      
      complete_df = compute_as_events(temp_df2, sample_list) %>%
                      purrr::reduce(left_join, by=c("feature_id", "gene_id"))
      
      complete_df = complete_df[order(complete_df$feature_id), ]
      rownames(complete_df) = NULL

      return(complete_df)
    }

######################################################################################
### MAIN
######################################################################################

compute_data_format = function(x, 
                               input_format, 
                               select, 
                               event_list, 
                               output
                               ) {  

  if (input_format=="DS") {
    
    cat("\nInput data in DRIMSeq format...\n")
    names(x)[names(x) == "gene"] = "gene_id"
    return(x)
    
  } else if (input_format=="JB" & select==FALSE){

    cat("\nFormatting JuncBase data to DRIMSeq format...\n")
    temp_df = as.data.frame(x) %>% 
      dplyr::rename(
        junc_ann=`#Contains_Novel_or_Only_Known(Annotated)_Junctions`
      ) 

    drimseq_data = juncbase2drimseq_format(temp_df)
    drimseq_data= drimseq_data[!duplicated(drimseq_data$feature_id), ]

    data.table::fwrite(drimseq_data, 
                       paste0(output, 
                        "jb_table_AS_excl_incl_counts_lenNorm.DrimSeq_format.tsv"), 
                       sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    return(drimseq_data)

  } else if (input_format=="JB" & select==TRUE) {

    cat("\nFormatting JuncBase data to DRIMSeq format...\n")
    temp_df = as.data.frame(x) %>% 
      dplyr::rename(
        junc_ann=`#Contains_Novel_or_Only_Known(Annotated)_Junctions`
        ) %>%
      dplyr::filter(as_event_type %in% c(event_list))

    drimseq_data = juncbase2drimseq_format(temp_df)
    drimseq_data= drimseq_data[!duplicated(drimseq_data$feature_id), ]

    data.table::fwrite(drimseq_data, 
                       paste0(output, 
                        "jb_table_AS_excl_incl_counts_lenNorm.DrimSeq_format.tsv"), 
                       sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    cat(paste0("\nComputing ", event_list, " specific analysis...\n"))
    return(drimseq_data)

  } else {
    stop("Input data format must be specified!") 
  }
}


