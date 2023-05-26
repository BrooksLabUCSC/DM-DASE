
######################################################################################
### Filter CDS AS events only
######################################################################################

######################################################################################
### Clean data
######################################################################################

clean_feature_ids = function(x) {

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

clean_gtf_data = function(x, 
                          release
                          ) {
	
  release_date = paste0("##date: ", release)
	temp_df = x %>% 
    dplyr::select(chr=.data[[release_date]], source=V2, feature=V3, 
                  start=V4, end=V5, strand=V7, attribute=V9) %>% 
  	dplyr::filter(source=="ENSEMBL")%>%
    dplyr::select(-source) %>% 
  	tidyr::separate(attribute, 
                    into=c("ensembl", "transcript_id", "gene_type", "gene_name"), 
                    sep=";")
  
  temp.clean = as.data.frame(clean_feature_ids(temp_df)) %>%
    tidyr::separate(ensembl, into=c("ensembl_id", "empty"), remove=FALSE) %>%
    tidyr::separate(gene_name, into=c("name", "symbol2"), remove=FALSE) %>%
    dplyr::select(-empty, -name)
	
  return(temp.clean)
}

######################################################################################
### Clean DrimSeq data
######################################################################################

clean_ds_data = function(x) {

	temp_df = x %>%
    tidyr::separate(gene_id, 
                    into=c("ensembl", "regions"), 
                    sep=";", remove=FALSE) %>%
    tidyr::separate(regions, 
                    into=c("event", "empty1"), 
                    sep=",", remove=TRUE) %>%
    tidyr::separate(event, 
                    into=c("chr", "positions"), 
                    sep=":", remove=FALSE) %>%
    tidyr::separate(positions, 
                    into=c("start", "stop"), 
                    sep="-", remove=FALSE) %>%
		tidyr::separate(ensembl, 
                    into=c("ensembl", "ensembl2"), 
                    sep=",", remove=TRUE) %>%
		tidyr::separate(ensembl, 
                    into=c("ensembl_id", "empty2"), 
                    remove=FALSE) %>%
		dplyr::select(-empty2, -ensembl2, -empty1, -positions, -event) 
  
  temp_df = temp_df %>%
    	dplyr::mutate(across(c("ensembl", "ensembl_id"), 
                    ~ifelse(.=="", NA, as.character(.))))
  	
  return(temp_df)
  	
}

clean_mesa_data = function(x) {

  temp_df = x %>% 
    tidyr::separate(gene_id, 
                    into=c("empty1", "empty2", "event"), 
                    sep="_", remove=FALSE) %>%
    tidyr::separate(event, 
                    into=c("chr", "region", "strand"), 
                    sep=":", remove=TRUE) %>%
    tidyr::separate(region, 
                    into=c("start", "stop"), 
                    sep="-", remove=TRUE) %>%
    dplyr::select(-empty1, -empty2)
  
  return(temp_df)

}

######################################################################################
### Infer gene names
######################################################################################

filter_funct_as_events = function(x,
                                  y,
                                  format,
                                  release,
                                  output
                                  ) {

    y = clean_gtf_data(y, release)
    y$chr = as.character(y$chr)
    y$start = as.numeric(y$start)
    y$end = as.numeric(y$end)

    if (format=="DS") {

      x = clean_mesa_data(x)
      data.split = split(x, f=x$feature_id) 
      
      cat("\nMapping functional AS events to annotation reference...\n")
      ann.split = lapply(data.split,

                        function(x) {
                          
                          ensembl = x$ensembl
                            chromosome = as.character(x$chr)
                            start_pos = as.numeric(x$start)
                            stop_pos = as.numeric(x$stop)
                            event_region = paste0(chromosome, ":", start_pos, "-", stop_pos)
                            temp_ann = y %>% dplyr::filter(chr == chromosome, 
                                                           start <= start_pos, 
                                                           end >= stop_pos, 
                                                           strand == strand)

                            if (!nrow(temp_ann)==0) {
                              
                              temp_df = as.data.frame(head(temp_ann, 1))
                              gene_region = paste0(as.character(temp_df$chr), ":", 
                                                   as.character(temp_df$start), "-", 
                                                   as.character(temp_df$end))
                              x$gene_id = gene_region
                              return(x)
                            } 
                        }
                    )
      functional_data = as.data.frame(rbindlist(ann.split, 
                                      use.names=TRUE, fill=TRUE,
                                       idcol="id")) %>%
        dplyr::select(-chr, -start, -stop, -strand, -id)

    } else if (format=="JB") {

      x = clean_ds_data(x)
      data.split = split(x, f=x$gene_id)

      cat("\nMapping functional AS events to annotation reference...\n")
      ann.split = lapply(data.split,

                        function(x) {

                          ensembl = x$ensembl
                            chromosome = as.character(x$chr)
                            start_pos = as.numeric(x$start)
                            stop_pos = as.numeric(x$stop)
                            strand = as.character(x$strand)
                            event_region = paste0(chromosome, ":", start_pos, "-", stop_pos)
                            temp_ann = y %>% dplyr::filter(chr == chromosome, 
                                                           start <= start_pos, 
                                                           end >= stop_pos) 
                            
                            if (!nrow(temp_ann)==0) {
                              return(x)
                            } 
                        }
                    )
      functional_data = as.data.frame(rbindlist(ann.split, 
                                      use.names=TRUE, 
                                      fill=TRUE, 
                                      idcol="id")) %>%
              dplyr::select(-ensembl, -ensembl_id, -chr, -start, -stop, -id)
    }

  	data.table::fwrite(functional_data, 
  		paste0(output, "functional_as_events_drimseq.tsv"), 
  		sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    cat("\nFunctional AS events filtered.\n")
  	
    return(functional_data)
  	
}
#annotate_transcripts(input.data, gtf.data, format, release, output) 



