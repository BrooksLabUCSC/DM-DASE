
######################################################################################
### Gene set enrichment analysis with GSEAPY
######################################################################################

### Author: Carlos Arevalo
### Email: carevalo0170@gmail.com

### PROGRAM DESCRIPTION

### Program takes as input a DTU table from DrimSeq or any other differential expression/usage analysis 
### containing "symbol", "gene_padj" and "log2fcp" and perform gene-set enrichment analysis (GSEA) and 
### pre-rank GSEA (pGSEA) using the GSEAPY package. Program outputs a main table for GSEA and pGSEA, 
### respectively. Analysis can be perform on both individual or a list of libraries. Some Human test 
### libabries include 'GO_Biological_Process_2021', 'GO_Molecular_Function_2021', 'GO_Cellular_Component_2021', 
### 'KEGG_2021_Human', 'Reactome_2022', and 'MSigDB_Hallmark_2020'

### INSTALLATION

### pip install gseapy 
### conda install -c bioconda gseapy

### PROGRAM USAGE

#python .../gsea/compute_gsea.py \
#	-i drimseq_deltaPSI_padj_results.txt \
#	-s "Human" \
#	-c "global" \
#	-l GO_Biological_Process_2021 GO_Molecular_Function_2021 GO_Cellular_Component_2021 KEGG_2021_Human Reactome_2022 MSigDB_Hallmark_2020 \
#	-t 1.5 \
#	-e 0.5 \
#	-n 10 \
#	-o /gsea_output/

######################################################################################

import gseapy as gp
from gseapy.plot import barplot, dotplot
from gseapy import enrichment_map
from gseapy.plot import gseaplot
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np
import sys
import os
import argparse
import warnings

class CommandLine():

	def __init__(self, inOpts=None):
		
		self.parser = argparse.ArgumentParser(
			description = "compute_gsea.py - computes GSEA using GSEAPY package",
			epilog = "Program epilog - please provide corrections and implementations for the program",
			add_help = True,
			prefix_chars = "-",
			usage = "%(prog)s [options] -option1[default] <input>output")

		self.parser.add_argument("-i", "--input", type=str, required=True, help='Input data in the form of table/data-frame')
		self.parser.add_argument("-s", "--organism", type=str, required=True, help='Organism for analysis')
		self.parser.add_argument("-c", "--condition", type=str, required=True, help='Condition label')
		self.parser.add_argument("-l", "--library", type=str, nargs='+', required=True, help='Library list')
		self.parser.add_argument("-t", "--threshold", type=float, required=True, help='Log2 FC threshold for significant genes selection')
		self.parser.add_argument("-e", "--enrich_threshold", type=float, required=True, help='Threshold for enrichment analysis')
		self.parser.add_argument("-p", "--pval_threshold", type=float, required=True, help='Log2 FC threshold for significant genes selection')
		self.parser.add_argument("-n", "--terms", type=int, required=True, help='Number of top terms to plot')
		self.parser.add_argument("-o", "--output", type=str, required=True, help='Output directory for results')
		
		if inOpts is None:
			self.args = self.parser.parse_args()
		else:
			self.args = self.parser.parse_args(inOpts)

class computeGSEA():

	def readData(inFile, threshold, pval):
		'''
		Read input data
		'''
		warnings.filterwarnings("ignore")
		df = pd.read_csv(inFile, header=None, sep="\t")
		df.columns = df.iloc[0]
		df = df[1:]
		df.reset_index(drop=True, inplace=True)
		
		temp_df = df[["symbol", "gene_padj", "log2fc", "feature_id"]]
		temp_df['gene_padj'] = temp_df['gene_padj'].astype(float)
		temp_df['log2fc'] = temp_df['log2fc'].astype(float)
		
		temp_sig = temp_df.loc[temp_df.gene_padj < pval] #0.05] 
		down_df = temp_sig[(temp_sig.log2fc < -threshold)]
		up_df = temp_sig[(temp_sig.log2fc > threshold)]
		data = pd.concat([up_df, down_df])
		
		return data

	def readPrerank(data, threshold, pval):
		'''
		Read input data
		'''
		warnings.filterwarnings("ignore")
		df = pd.read_csv(data, header=None, index_col=0, sep="\t")
		df.columns = df.iloc[0]
		df = df[1:]
		df.reset_index(drop=True, inplace=True)
		
		temp_df = df[["symbol", "gene_padj", "log2fc"]]
		temp_df['gene_padj'] = temp_df['gene_padj'].astype(float)
		temp_df['log2fc'] = temp_df['log2fc'].astype(float)
		
		temp_sig = temp_df.loc[temp_df.gene_padj < pval] #0.05]
		down_df = temp_sig[(temp_sig.log2fc < -threshold)]
		up_df = temp_sig[(temp_sig.log2fc > threshold)]
		data = pd.concat([up_df, down_df])

		data['rank'] = data['gene_padj']*data['log2fc']
		temp_sorted = data.sort_values('rank', ascending=False)
		uniq_df = temp_sorted.drop_duplicates(subset=['symbol'])

		data = pd.DataFrame()
		data[1] = list(uniq_df["rank"])
		data.index = list(uniq_df["symbol"])
		
		return data

	def readPositives(data, threshold, pval):
		'''
		Get genes with positive LFC values
		'''
		warnings.filterwarnings("ignore")
		df = pd.read_csv(data, header=None, index_col=0, sep="\t")
		df.columns = df.iloc[0]
		df = df[1:]
		df.reset_index(drop=True, inplace=True)
		temp_df = df[["symbol", "gene_padj", "log2fc"]]
		temp_df['gene_padj'] = temp_df['gene_padj'].astype(float)
		temp_df['log2fc'] = temp_df['log2fc'].astype(float)
		temp_sig = temp_df.loc[temp_df.gene_padj < pval] #0.05]
		pos_df = temp_sig[(temp_sig.log2fc > threshold)]
		pos_df['rank'] = pos_df['gene_padj']*pos_df['log2fc']
		temp_sorted_pos = pos_df.sort_values('rank', ascending=False)
		uniq_pos = temp_sorted_pos.drop_duplicates(subset=['symbol'])
		pos_data = pd.DataFrame()
		pos_data[1] = list(uniq_pos["rank"])
		pos_data.index = list(uniq_pos["symbol"])
		
		return pos_data

	def readNegatives(data, threshold, pval):
		'''
		Get genes with negative LFC values
		'''
		warnings.filterwarnings("ignore")
		df = pd.read_csv(data, header=None, index_col=0, sep="\t")
		df.columns = df.iloc[0]
		df = df[1:]
		df.reset_index(drop=True, inplace=True)
		temp_df = df[["symbol", "gene_padj", "log2fc"]]
		temp_df['gene_padj'] = temp_df['gene_padj'].astype(float)
		temp_df['log2fc'] = temp_df['log2fc'].astype(float)
		temp_sig = temp_df.loc[temp_df.gene_padj < pval] #0.05]
		down_df = temp_sig[(temp_sig.log2fc < -threshold)]
		down_df['rank'] = (-1)*down_df['gene_padj']*down_df['log2fc']
		temp_sorted = down_df.sort_values('rank', ascending=False)
		uniq_down = temp_sorted.drop_duplicates(subset=['symbol'])
		down_data = pd.DataFrame()
		down_data[1] = list(uniq_down["rank"])
		down_data.index = list(uniq_down["symbol"])
		
		return down_data

	def df2List(df):
		"""
		Convert dataframe or series to list
		"""
		warnings.filterwarnings("ignore")
		temp_df = pd.DataFrame()
		gene_list = pd.unique(list(df["symbol"])).tolist()
		gene_list = [x for x in gene_list if str(x) != 'nan']
		temp_df[0] = gene_list
		return(temp_df)

	def enrichR(gene_list, gene_set, organism, threshold, layout, output): 
		"""
		Peforms enrichr analysis
		"""
		enr = gp.enrichr(gene_list=gene_list, 
                 		 gene_sets=gene_set, 
                 		 organism=organism, 
                 		 no_plot=True,
                 		 cutoff=threshold, 
                 		 outdir='{output_dir}{label}_enrichr_analysis'.format(output_dir=output, label=layout)
                		)
		enr_results = enr.results
		enr_results = enr_results[enr_results["Adjusted P-value"] < 0.05]
		return(enr_results)

	def barPlot(df, title, top_term, layout, output):
		"""
		Plots pathways barplot
		"""
		plot = barplot(df, column="Adjusted P-value", size=10, 
					   top_term=top_term, title=title)
		plot.figure.savefig('{output_dir}{label}_enrichment_barplot.png'.format(
								output_dir=output, label=layout),
							bbox_inches="tight",
							dpi=600)

	def dotPlot(df, title, top_term, layout, output):
		"""
		Plots pathways dotplot
		"""
		plot2 = dotplot(df, size=10, top_term=top_term, title=title, 
						marker='o', show_ring=False, cmap="seismic",)
		plot2.figure.savefig(
			'{output_dir}{label}_enrichment_dotplot.png'.format(
				output_dir=output, label=layout),
			bbox_inches="tight",
			dpi=600)

	def enrichmentPlot(df, top_term, output):
		"""
		Plot enrichment analysis for a defined number of terms
		"""
		results = df.sort_index().head()
		terms = df.Term
		for i in range(1, top_term):
			term = terms[i]
			plot = gseaplot(rank_metric=df.ranking, 
							term = df.Term[i],
							**df[terms[i]])
			plot.figure.savefig(
				'{output_dir}term_{label}_gsea_plot.png'.format(output_dir=output, label=term),
				bbox_inches="tight",
				dpi=600)

	def prerankGSEA(rank_df, gset, top_term, layout, output):
		"""
		Enrichr libraries are supported by prerank module. Just provide the name
		use 4 process to acceralate the permutation speed
		"""
		prerank = gp.prerank(rnk=rank_df, 
							 gene_sets=gset,
							 threads=4,
							 min_size=10,
                     		 max_size=1000,
                     		 processes=4,
                     		 permutation_num=100, 
                     		 ascending=False,
                     		 outdir='{output_dir}{label}_prerank_report'.format(output_dir=output, label=layout),
                     		 format='png', 
                     		 seed=6,
                     		 verbose=True)
		return prerank 

	def enrichmentMap(df, layout, output):
		"""
		Performs enrichment mapping
		"""
		nodes, edges = enrichment_map(df)
		graph = nx.from_pandas_edgelist(edges,
                            source='src_idx',
                            target='targ_idx',
                            edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
		fig, ax = plt.subplots(figsize=(7, 7))
		pos = nx.layout.spiral_layout(graph)
		nx.draw_networkx_nodes(graph,
                       	pos=pos,
                       	cmap=plt.cm.RdYlBu,
                       	node_color=list(nodes.NES),
                       	node_size=list(nodes.Hits_ratio*1000))
		nx.draw_networkx_labels(graph,
                        	pos=pos,
                        	labels=nodes.Term.to_dict())
		edge_weight = nx.get_edge_attributes(graph, 'jaccard_coef').values()
		nx.draw_networkx_edges(graph,
                       	pos=pos,
                       	width=list(map(lambda x: x*10, edge_weight)),
                       	edge_color='#CDDBD4')
		plt.savefig(
			'{output_dir}{label}_pca_projection.png'.format(output_dir=output, label=layout),
			bbox_inches="tight",
			dpi=600)

def main(incl=None):
	
	if incl is None:
		command_line = CommandLine()
		
		if command_line.args.input:
			
			output = command_line.args.output
			if not os.path.isdir(output):
				os.mkdir(output)

			inFile = command_line.args.input
			organism = command_line.args.organism
			library = command_line.args.library
			condition = command_line.args.condition
			lfc_threshold = command_line.args.threshold
			enrich_threshold = command_line.args.enrich_threshold
			pval_threshold = command_line.args.pval_threshold
			top_terms = command_line.args.terms
			sets = '_'.join(library)
			
			print("Computing global enrichment analysis...")
			input_data = computeGSEA.readData(inFile, threshold=lfc_threshold, pval=pval_threshold)
			gene_list = computeGSEA.df2List(input_data)
			enrich = computeGSEA.enrichR(gene_list=gene_list, 
										 gene_set=library, 
										 organism=organism, 
										 threshold=enrich_threshold,
										 layout=condition,
										 output=output)

			computeGSEA.barPlot(df=enrich, title="test", top_term=top_terms, 
								layout=condition, output=output)
			computeGSEA.dotPlot(df=enrich, title="test", top_term=top_terms, 
								layout=condition, output=output)
			prerank_df = computeGSEA.readPrerank(data=inFile, threshold=lfc_threshold, pval=pval_threshold)
			prerank = computeGSEA.prerankGSEA(rank_df=prerank_df, gset=library, 
											  top_term=top_terms, layout="all", output=output)
			prerank_res = prerank.res2d
			
			print("Computing positives and negatives only enrichment analysis...")
			pos_prerank_df = computeGSEA.readPositives(data=inFile, threshold=lfc_threshold, pval=pval_threshold)
			neg_prerank_df = computeGSEA.readNegatives(data=inFile, threshold=lfc_threshold, pval=pval_threshold)
			pos_prerank = computeGSEA.prerankGSEA(rank_df=pos_prerank_df, gset=library,
					 							  top_term=top_terms, layout="positive", output=output)
			neg_prerank = computeGSEA.prerankGSEA(rank_df=neg_prerank_df, gset=library,
					 							  top_term=top_terms, layout="negative", output=output)
			pos_prerank_res = pos_prerank.res2d
			neg_prerank_res = neg_prerank.res2d

			prerank_sig = prerank_res.loc[prerank_res["FDR q-val"]<0.05]
			pos_prerank_sig = pos_prerank_res.loc[pos_prerank_res["FDR q-val"]<0.05]
			neg_prerank_sig = neg_prerank_res.loc[neg_prerank_res["FDR q-val"]<0.05]

			if not prerank_sig.empty:
				computeGSEA.enrichmentMap(df=prerank, layout=condition, output=output)
				computeGSEA.enrichmentPlot(df=prerank, top_term=10, output=output)
			if prerank_sig.empty:
				print("\033[1mNo significant pre-rank terms found.\n")

			if not pos_prerank_sig.empty:
				computeGSEA.enrichmentMap(df=pos_prerank, layout=condition, output=output)
				computeGSEA.enrichmentPlot(df=pos_prerank, top_term=10, output=output)
			if pos_prerank_sig.empty:
				print("\033[1mNo significant positive pre-rank terms found.\n")
			
			if not neg_prerank_sig.empty:
				computeGSEA.enrichmentMap(df=neg_prerank, layout=condition, output=output)
				computeGSEA.enrichmentPlot(df=neg_prerank, top_term=10, output=output)
			if neg_prerank_sig.empty:
				print("\033[1mNo significant negative pre-rank terms found.\n")

			print("\033[1mGSEA and prerank GSEA analysis complete.\n")

if __name__ == '__main__':
    main()

