
######################################################################################
### Annotate unannotated and novel junctions using the genome reference file (.gtf)
######################################################################################

clean_gtf_feature_ids = function(x) {

  for (pos in 1:ncol(x)) {

    if (pos == 5) {
      x$ensembl = gsub("gene_id", "", as.character(x$ensembl))
      x$ensembl = gsub("transcript_id", "", as.character(x$ensembl))
      x$ensembl = gsub("gene_type", "", as.character(x$ensembl))
      x$ensembl = gsub("gene_name", "", as.character(x$ensembl))
      x$ensembl = gsub(" ", "", as.character(x$ensembl))
      x$ensembl = gsub('"', "", as.character(x$ensembl))
    } else if (pos == 6) {
      x$gene_type = gsub("gene_id", "", as.character(x$gene_type))
      x$gene_type = gsub("transcript_id", "", as.character(x$gene_type))
      x$gene_type = gsub("gene_type", "", as.character(x$gene_type))
      x$gene_type = gsub("gene_name", "", as.character(x$gene_type))
      x$gene_type = gsub('"', "", as.character(x$gene_type))
    } else if (pos == 7) {
      x$gene_name = gsub("gene_id", "", as.character(x$gene_name))
      x$gene_name = gsub("transcript_id", "", as.character(x$gene_name))
      x$gene_name = gsub("gene_type", "", as.character(x$gene_name))
      x$gene_name = gsub("gene_name", "", as.character(x$gene_name))
      x$gene_name = gsub('"', "", as.character(x$gene_name))
    }
  }

  x$ensembl = as.character(x$ensembl)
  x$gene_type = as.factor(x$gene_type)
  x$gene_name = as.factor(x$gene_name)
  
  return(x)
}

######################################################################################
### Clean GTF data
######################################################################################

clean_gtf_data = function(x, release) {

	release_date = paste0("##date: ", release) 
	temp_df = x %>% 
    dplyr::select(chr=release_date, source=V2, feature=V3, 
                  start=V4, end=V5, strand=V7, attribute=V9) %>% 
  	dplyr::filter(source=="ENSEMBL")%>%
    dplyr::select(-source) %>% 
  	tidyr::separate(attribute, 
                    into=c("ensembl", "transcript_id", 
                           "gene_type", "gene_name"), 
                    sep=";"
                    )
  
  temp.clean = as.data.frame(clean_gtf_feature_ids(temp_df)) %>%
    tidyr::separate(ensembl, 
                    into=c("ensembl_id", "empty"), 
                    remove=FALSE) %>%
    tidyr::separate(gene_name, 
                    into=c("name", "symbol2"), 
                    remove=FALSE) %>%
    dplyr::select(-empty, -name)
	
  return(temp.clean)
}


######################################################################################
### Clean DrimSeq data
######################################################################################


clean_drimseq_data = function(x) {

	temp_df = x %>%
		tidyr::separate(gene_id, 
                    into=c("ensembl", "chr", "start", "stop"), 
                    sep="_", remove=FALSE) %>%
    tidyr::separate(stop, 
                    into=c("stop", "empty"), 
                    sep=",", remove=TRUE) %>%
    dplyr::select(-empty)
  
  temp_df = temp_df %>% 
    dplyr::mutate(across(c("ensembl"), 
                  ~ifelse(.=="", NA, as.character(.))))
  
  return(temp_df)
}

clean_mesa_data2 = function(x) {

  x$feature_id = gsub(":", "_", as.character(x$gene_id))
  x$feature_id = gsub("-", "_", as.character(x$gene_id))

  temp_sep = x %>% 
    tidyr::separate(feature_id, 
      into=c("chr", "start", "stop"), sep="_", remove=FALSE)
  
  return(temp_sep)
}

######################################################################################
### Infer gene names
######################################################################################

annotate_transcripts = function(x,
                                y,
                                annotations,
                                format,
                                release,
                                ...
                                ) {

  y = clean_gtf_data(y, release)
  annotations = annotations %>% dplyr::rename(ensembl_id=gene_id)
  y$chr = as.character(y$chr)
  y$start = as.numeric(y$start)
  y$end = as.numeric(y$end)

  if (format=="DS") {

    x = clean_mesa_data2(x)

    cat("\nMapping AS events to reference...\n")
    
    data.split = split(x, f=x$gene_id) 
    ann.split = lapply(data.split,

                          function(x) {

                            chromosome = as.character(x$chr)
                            start_pos = as.numeric(x$start)
                            stop_pos = as.numeric(x$stop)
                            event_region = paste0(chromosome, ":", start_pos, "-", stop_pos)
                            temp_ann = y %>% 
                              dplyr::filter(chr==chromosome, start<=start_pos, end>=stop_pos) 
                       
                            if (!nrow(temp_ann)==0) {

                              temp_df = as.data.frame(head(temp_ann, 1))
                              x$ensembl = temp_df$ensembl
                              x$ensembl_id = temp_df$ensembl_id
                              x$symbol = temp_df$symbol2
                          
                            } else if (nrow(temp_ann)==0) {
                          
                              x$ensembl = event_region
                              x$ensembl_id = event_region
                              x$symbol = event_region
                          
                            }
                            return(x)
                          }
                      )
    ann.temp = rbindlist(ann.split, use.names=TRUE, fill=TRUE, idcol="id") %>%
      dplyr::select(-id, -chr, -start, -stop)

  } else if (format=="JB") {

    x = clean_drimseq_data(x)

    cat("\nMapping unannotated events to reference...\n")

    data.split = split(x, f=x$gene_id) 
    ann.split = lapply(data.split,

                     	    function(x) {

                       	    chromosome = as.character(x$chr)
                       	    start_pos = as.numeric(x$start)
                            stop_pos = as.numeric(x$stop)
                            event_region = paste0(chromosome, ":", start_pos, "-", stop_pos)
                            temp_ann = y %>% 
                              dplyr::filter(chr==chromosome, start<=start_pos, end>=stop_pos) 
                       
                       	    if (!nrow(temp_ann)==0) {
                              
                              temp_df = as.data.frame(head(temp_ann, 1))
                              x$ensembl = temp_df$ensembl
                              x$ensembl_id = temp_df$ensembl_id
                              x$symbol = temp_df$symbol2
                          
                            } else if (nrow(temp_ann)==0) {
                          
                              x$ensembl = event_region
                              x$ensembl_id = event_region
                              x$symbol = event_region
                          
                            }
                       	    return(x)
                       	  }
                      )
	  ann.temp = rbindlist(ann.split, use.names=TRUE, fill=TRUE, idcol="id") %>%
      dplyr::select(-id, -chr, -start, -stop)

  }

  ann.data = dplyr::left_join(ann.temp, annotations, by="ensembl_id")

  cat("\nAS events have been annotated.\n")

  return(ann.data)
}


