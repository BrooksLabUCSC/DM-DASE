
## Motivation

Identifying differentially spliced genes and transcripts in cancer requires prior identification and annotation of AS events in tumors. To identify and annotate novel AS events (at either the cancer or pan-cancer or single-cell level), it is required to quantify and classify AS events from sequencing reads. JuncBASE (Junction-Based Analysis of Splicing Events) is a tool to quantify and classify AS events using splice junctions reads from RNA-Seq alignments and annotated exon coordinates. Another approach includes detecting and quantifying AS events by mutually exclusive junctions. MESA (Mutually Exclusive Splicing Analysis) computes percent splice in (PSI) values applying a mutually exclusive junction approach. To identify transcripts or novel AS events with differential usage or splicing, it is necessary to perform the correct differential analysis based on sample size and counts distribution. DRIMSeq performs a Dirichlet multinomial distribution approach to identify significant genes and transcripts with differential usage. In addition, it could be implemented to perform differential splicing analysis (DS). The DM-DASE model workflow requires the quantification and classification of AS events using either JuncBASE or MESA, and takes a PSI counts file from either JuncBASE or MESA to performe DS analysis at the global or event-type level.

</h1>
<img src= "https://github.com/caeareva/DM-DASE/blob/616157d31dc6ea292eab4abaaa7167a2e2781f8c/figures/as_event_types.png"
</h1>

## JuncBASE

JuncBASE takes in aligned and trimmed reads in the form of BAM files.

##### STEP 0: Remove non-standard chromosomes on each sample
```bash
for index in sample1 sample2 sample3 sample4 sample5 sample6 sample7
	do
	echo "Processing "${index}".sortedByCoord.out.bam ..."
	samtools view -b "${index}".sortedByCoord.out.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > /Path_to_Dir/juncbase_trimmed_bams/"${index}".sortedByCoord.trimmed.bam
	echo "Processing "${index}".sortedByCoord.out.bam completed..."
done
```

#### Build de novo transcript DB 

```bash
python runCufflinks.py \
	-i juncbase_samp2bam.tsv \
	--txt_ref UCSC.GRCh38.knownGene_w_gene_symbol_EnsemblChr.gtf \
	--out_dir /path/to/cufflinks_denovo/ \
	--num_processes 2

```

```bash
python runCuffmerge.py \
	-f cufflinks_denovo_IlluminaHBM2.0.txt \
	--no_single_exons \
	-g UCSC.GRCh38.knownGene_w_gene_symbol_EnsemblChr.gtf \
	-o /cuffmerge_out \
	-s Homo_sapiens_assembly_GRCh38.fasta \
	--tmp_dir /tmp_dir
```

Initialize
```bash
python build_DB_FromGTF.py \
	-g gencode.GRCh38.annotation.chr.basic.gtf  \
	-d sqlitedb \
	--sqlite_db_dir /output_dir \
	--initialize
```

Build de novo DB
```bash
python build_DB_FromGTF.py \
	-g gencode.GRCh38.annotation.chr.basic.gtf  \
	-d sqlitedb \
	--sqlite_db_dir /output_dir
```

##### Plot entropy scores
```bash
python plotEntropyScores.py \
	-s juncbase_trimmed_bams/sample_test.bam \
	-o /juncbase_data/sample_test \
	--known_junctions /sqlite_db_dir/genecode.GRCh38.as.an_intron_coordinate.txt

```


##### STEP 1: Process SAM/BAMs files
```bash
python run_preProcess_by_chr_step1.py \
	-i juncbase_samp2bam.tsv \
	-o /juncbase_trimmed_bams/ \
	--preProcess_options "–unique -j .../genecode.GRCh38.as.an_intron_coordinate.txt" \
	-p 20
```


##### STEP 1B: Desambiguate splice juntion orientations
```bash
python disambiguate_junctions.py \
	-i /juncbase_trimmed_bams/ \
	-g GRCh38.v1.2020_04_01.fa \
	--by_chr --majority_rules
```


##### STEP 2: Identify all junctions
```bash
python preProcess_getASEventReadCounts_by_chr_step2.py \
	-i /juncbase_trimmed_bams/ \
	--by_chr
```


##### STEP 3: Create exon-intron junction count files
```bash
python run_preProcess_step3_by_chr.py \
	--input_dir /juncbase_trimmed_bams/ \
	--num_processes 20 --force
```


##### STEP 4: Create a pseudo/"all junction" sample
```bash
python createPseudoSample.py \
	-i /juncbase_trimmed_bams \
	-s sample1 --by_chr
```

##### STEP 5: Identify AS events and quantify events from each sample
```bash
python run_getASEventReadCounts_multiSample.py \
	-s sample1,sample2,sample3,sample4,sample5,sample6,sample7 \
	-i /juncbase_trimmed_bams \
	-o /getASEventReadCounts \
	--sqlite_db_dir /output_dir/sqlite_db_dir \
	--txt_db1 genecode.v33_pr.as.an_db1 \
	--txt_db2 genecode.v33_basic.ann_db1 \
	--txt_db3 genecode.v33_pr.as.an_db1 \
	--jcn_seq_len 240 -p 20 --by_chr
```


##### STEP 6 Part 1: Create tables of raw and length-normalized read counts of exclusion and inclusion isoforms
```bash
python run_createAS_CountTables.py \
	-d /getASEventReadCounts/ \
	-i /juncbase_trimmed_bams/ \
	--jcn_seq_len 240 \
	-s sample1,sample2,sample3,sample4,sample5,sample6,sample7 \
	--num_processes 30
```


##### STEP 6 Part 2: Combine tables of raw and length-normalized read counts
```bash
python combine_createAS_CountTables_by_chr.py \
	-d /getASEventReadCounts/ \
	-o /juncbase_tables/
```

## MESA

Create a virtual env and install MESA
```bash
### virtual env
mkdir mesaEnv
python3 -m pip install --user virtualenv 
python3 -m virtualenv mesaEnv
source mesaEnv/bin/activate

### installation
git clone https://github.com/BrooksLabUCSC/mesa.git
cd mesa/
pip install .
mesa --version
mesa -h
deactivate
```

Convert BAMs to BEDs with junctions
```bash
mesa bam_to_junc_bed \
    -m manifest_bam.tsv \
	-o /MESA_HG38/ \
	-n 10 -s inferCombine \
	-a gencode.v38.annotation.gtf \
	-g GRCh38.u2af1_fix.v1.2020_04_01.fa
```

Quantify MESA events
```bash
mesa quant \
    -m manifest_bed.tsv \
	-o /MESA_U2AF1_HG38/ \
	--drim --maxLength 50000 --minLength 50 --minOverhang 5 \
	--minUnique 5 --lowCoverageNan --minEntropy 1
```

## DASE Model

Compute global DS/DU analysis from JuncBase counts
```bash
conda create -n "temp_rEnv" r-essentials r-base
conda activate temp_rEnv

Rscript .../src/dase/compute_dmDASE_main_model.R  \
   -i JB_AS_exclusion_inclusion_counts_lenNorm.txt \
   -m metadata.tsv \
   -g gencode.v38.annotation.gtf \
   -a ensembl_annotation.csv \
   -r "2021-03-12" \
   -f "JB" \
   -s FALSE \
   -s1 12 \
   -s2 6 \
   -c TRUE \
   -c1 "MT" \
   -c2 "WT" \
   -c3 "condition1" "condition2" "condition3" \
   -t "DU" \
   -n 12 \
   -b FALSE \
   -l "condition" \
   -o /PATH_TO_OUTPUT_DIR/global_analysis/
```

AS event-specific DS/DU analysis from JuncBase counts. Some event types include `cassette`, `intron_retention`, `alternative_acceptor`, `alternative_donor`, `jcn_only_AD`, `jcn_only_AA`, `coord_cassette`, `mutually_exclusive`, `alternative_last_exon`, and `alternative_first_exon`.
```bash
Rscript .../src/dase/compute_dmDASE_main_model.R  \
   -i JB_AS_exclusion_inclusion_counts_lenNorm.txt \
   -m metadata.tsv \
   -g gencode.v38.annotation.gtf \
   -a ensembl_annotation.csv \
   -r "2021-03-12" \
   -f "JB" \
   -s TRUE \
   -e "cassette" \
   -s1 18 \
   -s2 9 \
   -c TRUE \
   -c1 "MT" \
   -c2 "WT" \
   -c3 "condition1" "condition2" "condition3" \
   -t "DU" \
   -n 8 \
   -b FALSE \
   -l "condition" \
   -o /PATH_TO_OUTPUT_DIR/casette_events/
```

Compute DS/DU analysis from MESA
```bash
Rscript .../src/dase/compute_dmDASE_main_model.R \
   -i MESA_drimTable.tsv \
   -m metadata.tsv \
   -g gencode.v38.annotation.gtf \
   -a ensembl_annotation.csv \
   -r "2021-03-12" \
   -f "DS" \
   -s FALSE \
   -s1 18 \
   -s2 9 \
   -c TRUE \
   -c1 "MT" \
   -c2 "WT" \
   -c3 "condition1" "condition2" "condition3" \
   -t "DU" \
   -n 8 \
   -b FALSE \
   -l "condition" \
   -o /PATH_TO_OUTPUT_DIR/MESA_HG38/
```

Quantify significant AS events

```bash
Rscript quantify_sig_as_events_proportions.R \
  -i drimseq_gene_diffusage_results_annotated.tsv \
  -p 0.1 \
  -e global \
  -o /juncbase_global/
```

## Visualizations 

Plot results for significant genes and transcripts, and global DU/DS analysis. `qqplot_main_program.R` takes as input gene or transcript (txs) level DU results with the screened `p-adj` values. When runing results from JuncBASE, program provides the option to specify which AS event(s) are analyzing (i.e. "global", "intron_retention", "casette", etc). `volcano_program.R` takes as input DM-DASE results with `p-adj` and `LFC` values.

Gene-level QQ plot
```bash
# Global analysis:
Rscript qqplot_main_program.R \
   -i drimseq_gene_diffusage_results.tsv \
   -o /output/ \
   -t 10e-8 \
   -l "gene" \
   -e "global"

# Event-level analysis:
Rscript qqplot_main_program.R \
   -i drimseq_gene_diffusage_results.tsv \
   -o /output/ \
   -t 10e-8 \
   -l "gene" \
   -e "global
```

Transcript-level QQ plot
```bash
Rscript qqplot_main_program.R \
   -i drimseq_transcript_diffusage_results_annotated.tsv \
   -o output/ \
   -t 10e-8 \
   -l "txs" \
   -e "global"
```

Volcano plots: MT/WT volcano plot
```bash
# Global analysis:
Rscript volcano_program.R \
  -i drimseq_lfc_padj_results.tsv \
  -l "global" \
  -c "MT/WT" \
  -t "0.01" \
  -f "1.5" \
  -o /output/

# Event-level analysis:
Rscript volcano_program.R \
  -i /casette_events/drimseq_lfc_padj_results.tsv \
  -l "casette" \
  -c "MT/WT" \
  -t "0.01" \
  -f "1.5" \
  -o /casette_events/
  
# For comparisons:
Rscript volcano_program.R \
  -i /condition1/condition1_MT_WT_DS_LFC_padj_results.tsv \
  -l "mut_wt_condition1" \
  -c "Condition1 MT/WT" \
  -t "0.01" \
  -f "1.5" \
  -o /condition1/
```


## GSEA 

Perform GSEA and/or prerank-GSEA on DM-DASE output tables. If running the prerank analysis, the program requires a table with `p-adj` and `LFC` values to rank genes based on significance and expression. Some of the latest gene set libraries include `MSigDB_Hallmark_2020`, `Reactome_2022`, `GO_Biological_Process_2021`, `GO_Cellular_Component_2021`, and `KEGG_2021_Human`

```bash
python compute_gsea_program.py \
	-i /global_analysis/drimseq_lfc_padj_results.tsv \
	-s "Human" \
	-c "global" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 1.5 \
	-e 0.5 \
	-n 10 \
	-o /gsea_global/
	
# For comparisons:
python compute_gsea_program.py \
	-i /condition1/condition1_mut_wt_DS_LFC_padj_results.tsv \
	-s "Human" \
	-c "condition1_MT_WT" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 0 \
	-e 0 \
	-p 1 \
	-n 10 \
	-o /gsea_condition1/
```

Prerank analysis barplot 
```bash
Rscript gsea_barplot_program.R \
   -i /gsea_global/gseapy.gene_set.prerank.report.csv \
   -p "FDR q-val" \
   -e "global" \
   -s "all" \
   -t 0.1 \
   -m "prerank" \
   -o /gsea_global/prerank_report/
```

Enrichment analysis barplot 
```bash
Rscript gsea_barplot_program.R \
   -i /gsea_global//GO_Biological_Process_2021.Human.enrichr.reports.txt \
   -p "Adjusted P-value" \
   -e "global" \
   -s "all" \
   -t 0.1 \
   -m "enrichment" \
   -o /gsea_global/prerank_report/
```

Multiple prerank libraries barplot 
```bash
Rscript gsea_barplot_program.R \
   -i /gsea_output/casette/global_enrichr_analysis/GO_Biological_Process_2021.Human.enrichr.reports.txt \
   -i2 /gsea_output/casette/global_enrichr_analysis/GO_Cellular_Component_2021.Human.enrichr.reports.txt \
   -i3 /gsea_output/casette/global_enrichr_analysis/GO_Molecular_Function_2021.Human.enrichr.reports.txt \
   -i4 /gsea_output/casette/global_enrichr_analysis/MSigDB_Hallmark_2020.Human.enrichr.reports.txt \
   -p "Adjusted P-value" \
   -e "global" \
   -s "all" \
   -t 0.05 \
   -m "enrichment" \
   -g TRUE \
   -o /gsea_global/prerank_report/
```

Heatmap summary
```bash
Rscript gsea_heatmap_program.R \
   -i .../gsea_global/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../gsea_condition1/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../gsea_condition2/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../gsea_condition3/prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "global" \
   -c2 "condition1" \
   -c3 "condition2" \
   -c4 "condition3" \
   -e global \
   -s MSigDB_Hallmark_2020 \
   -t 1 \
   -pw 8.5\
   -ph 12 \
   -o .../gsea_global/
```

Genotypes summary heatmap
```bash
Rscript gsea_heatmap_program.R \
   -i .../gsea_condition1/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../gsea_condition1/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../gsea_condition2/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../gsea_condition2/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i5 .../gsea_condition3/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i6 .../gsea_condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "Condition1_pos" \
   -c2 "Condition1_neg" \
   -c3 "Condition2_pos" \
   -c4 "Condition2_neg" \
   -c5 "Condition3_pos" \
   -c6 "Condition4_neg" \
   -e "global" \
   -s MSigDB_Hallmark_2020 \
   -t 1 \
   -pw 11 \
   -ph 14 \
   -o .../gsea_global/
```
