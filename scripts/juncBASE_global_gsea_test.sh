
#################################################################################################
### JuncBASE Global GSEA
#################################################################################################

# Global analysis
python compute_gsea_program.py \
	-i drimseq_lfc_padj_results.tsv \
	-s "Human" \
	-c "global" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 0 \
	-e 0 \
	-p 1 \
	-n 10 \
	-o /PATH_TO_OUTPUT_DIR/

# Condition1 MT/WT
python compute_gsea_program.py \
	-i condition1_mut_wt_DS_LFC_padj_results.tsv \
	-s "Human" \
	-c "condition1_MT_WT" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 0 \
	-e 0 \
	-p 1 \
	-n 10 \
	-o /PATH_TO_OUTPUT_DIR/

# Condition2 MT/WT
python compute_gsea_program.py \
	-i condition2_mut_wt_DS_LFC_padj_results.tsv \
	-s "Human" \
	-c "condition2_MT_WT" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 0 \
	-e 0 \
	-p 1 \
	-n 10 \
	-o /PATH_TO_OUTPUT_DIR/

# Condition3 MT/WT
python compute_gsea_program.py \
	-i condition3_mut_wt_DS_LFC_padj_results.tsv \
	-s "Human" \
	-c "condition3_MT_WT" \
	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
	-t 0 \
	-e 0 \
	-p 1 \
	-n 10 \
	-o /PATH_TO_OUTPUT_DIR/


#################################################################################################
### Heatmap
#################################################################################################

# Global analysis
Rscript gsea_heatmap_program.R \
   -i .../jb_global/gsea_mt_wt/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../jb_global/condition1_mut_wt/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../jb_global/condition2_mut_wt/prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../jb_global/condition2_mut_wt/prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "global" \
   -c2 "Condition1" \
   -c3 "Condition2" \
   -c4 "Condition3" \
   -e global \
   -s MSigDB_Hallmark_2020 \
   -t 1 \
   -pw 8.5\
   -ph 12 \
   -o /PATH_TO_OUTPUT_DIR/


Rscript gsea_heatmap_program.R \
   -i .../Condition1/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../Condition1/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../Condition2/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i5 .../Condition3/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i6 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "genotype1_pos" \
   -c2 "genotype1_neg" \
   -c3 "genotype2_pos" \
   -c4 "genotype2_neg" \
   -c5 "genotype3_pos" \
   -c6 "genotype4_neg" \
   -e global \
   -s MSigDB_Hallmark_2020 \
   -t 1 \
   -pw 11 \
   -ph 14 \
   -o /PATH_TO_OUTPUT_DIR/

Rscript gsea_heatmap_program.R \
   -i .../Condition1/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../Condition1/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../Condition2/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i5 .../Condition3/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i6 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "genotype1_pos" \
   -c2 "genotype1_neg" \
   -c3 "genotype2_pos" \
   -c4 "genotype2_neg" \
   -c5 "genotype3_pos" \
   -c6 "genotype4_neg" \
   -e global \
   -s Reactome_2022 \
   -t 1 \
   -pw 11 \
   -ph 14 \
   -o /PATH_TO_OUTPUT_DIR/

Rscript gsea_heatmap_program.R \
   -i .../Condition1/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../Condition1/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../Condition2/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i5 .../Condition3/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i6 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "genotype1_pos" \
   -c2 "genotype1_neg" \
   -c3 "genotype2_pos" \
   -c4 "genotype2_neg" \
   -c5 "genotype3_pos" \
   -c6 "genotype4_neg" \
   -e global \
   -s GO_Cellular_Component_2021 \
   -t 1 \
   -pw 11 \
   -ph 14 \
   -o /PATH_TO_OUTPUT_DIR/

Rscript gsea_heatmap_program.R \
   -i .../Condition1/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i2 .../Condition1/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i3 .../Condition2/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i4 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i5 .../Condition3/positive_prerank_report/gseapy.gene_set.prerank.report.csv \
   -i6 .../Condition3/negative_prerank_report/gseapy.gene_set.prerank.report.csv \
   -c "genotype1_pos" \
   -c2 "genotype1_neg" \
   -c3 "genotype2_pos" \
   -c4 "genotype2_neg" \
   -c5 "genotype3_pos" \
   -c6 "genotype4_neg" \
   -e global \
   -s KEGG_2021_Human \
   -t 1 \
   -pw 11 \
   -ph 14 \
   -o /PATH_TO_OUTPUT_DIR/





