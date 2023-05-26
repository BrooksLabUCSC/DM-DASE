#!/usr/bin/env Rscript

################################################################################

# Author: Carlos Arevalo
# Email: carevalo0170@gmail.com

### DESCRIPTION

# Heatmap for single cell DGE analysis. Input a text/tsv file with at least the 
# following columns: gene, avg_log2FC, p_val, p_val_adj, cluster, and condition.

### PROGRAM USAGE

#Rscript custom_heatmap_program.R \
#  gene_expression_test_data.tsv \
#  /output/

################################################################################

options(warn=-1)
suppressMessages(if(!require(ggplot2)){install.packages("ggplot2")})
suppressMessages(if(!require(dplyr)){install.packages("dplyr")})
suppressMessages(if(!require(data.table)){install.packages("data.table")})
suppressMessages(if(!require(tidyverse)){install.packages("tidyverse")})
suppressMessages(if(!require(tidyr)){install.packages("tidyr")})
suppressMessages(if(!require(ggrepel)){install.packages("ggrepel")})
suppressMessages(if(!require(ggpubr)){install.packages("ggpubr")})
suppressMessages(if(!require(RColorBrewer)){install.packages("RColorBrewer")})
suppressMessages(if(!require(gtable)){install.packages("gtable")})
suppressMessages(if(!require(grid)){install.packages("grid")})

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
output = args[2]

main_df = fread(filename, header=TRUE)

### PRE-PROCESS DATA
main_df = main_df %>% 
  filter(p_val_adj < 0.1) %>%
  dplyr::mutate(
    z_score = ((avg_log2FC-mean(avg_log2FC, na.rm=TRUE)) / sd(avg_log2FC, na.rm=TRUE))
    ) %>% 
  #dplyr::mutate(z_score = ifelse(zscore < -2, -2, zscore)) %>% 
  select(gene, condition, cluster, z_score) 
 
### Join conditions and clusters as identifier: 
main_df$Ident = paste(main_df$condition, main_df$cluster, sep=" ")

### heatmap main function
custom_heatmap = function(x, 
						              ...
						              ) {
  
  ggplot(x, aes(x=reorder(gene, -z_score), y=cluster, fill=z_score,)) + 
  	geom_tile() + 
  	theme_bw(base_size=25, base_family="",
  			 base_line_size=1.5, base_rect_size=1.5) +
  	theme(panel.grid.major = element_blank(),
    	  panel.grid.minor = element_blank(),
    	  strip.switch.pad.grid = unit(0, "cm"),
    	  axis.ticks.x = element_blank(), 
    	  axis.ticks.y = element_blank(),
    	  axis.title.x = element_blank(),
    	  axis.title.y = element_blank(),
    	  axis.text.x = element_blank(),
    	  axis.text.y = element_text(hjust=1, vjust=0.5, size=20, color="black"),
    	  legend.text = element_text(color="black", size=20),
    	  legend.justification = "top",
    	  panel.spacing = unit(0.05, "lines"),
    	  strip.background = element_blank(),
    	  strip.text = element_blank(),
    	  strip.text.x = element_text(margin = margin(0.001, 0, 0.001, 0, "cm")),
    	  panel.background = element_rect(fill="white")) + 
  	scale_fill_gradientn(name="z-score",
  		colors=rev(brewer.pal(11, "RdBu")),
  		limits=c(-2.5, 2.5),
  		guide = guide_colorbar(label=TRUE, draw.ulim=TRUE, draw.llim=TRUE,
      		frame.colour="black", ticks=FALSE, nbin=10, frame.linewidth=1.5, 
      		label.position="right", title="z-score", direction="vertical")) + 
  	facet_grid(vars(cluster), vars(condition), drop=TRUE, scales="free", space="free")
}

### top annotations
annotations = function(x, 
					   ...
					   ) {
  
  ggplot(x,  aes(x=gene, y=cluster)) + 
  	facet_grid(.~condition, drop = TRUE, scales = "free", space = "free") + 
  	geom_rect(aes(fill = condition), xmin = -Inf, xmax = Inf, 
          ymin = -Inf, ymax = Inf) + theme_classic() + 
  	theme(strip.text = element_blank(),
          panel.spacing = unit(0.05, "lines"),
          strip.background = element_blank(),
          legend.title = element_text(color="black", size=20),
          legend.text = element_text(color="black", size=20), 
          strip.text.x = element_blank()) + 
  	scale_fill_manual(name = "DGE comparison", 
  		values = c("#FC4E2A", "#CCCCCC", "#F3A583", "#1770AB", "#39B185"),  
        label = c("Condition (MT/WT)",
          		  "Vehicle (MT/WT)", 
          		  "Condition_WT/Vehicle_WT",
          		  "Condition_MT/Vehicle_MT",
          		  "Condition_MT/Vehicle_WT")) + 
  	guides(color=guide_legend(verride.aes=list(size=5)))
}
  
### plot heatmap and column annotations
heatmap = custom_heatmap(main_df)
main_ann = annotations(main_df)
g1 = ggplotGrob(heatmap)
g2 = ggplotGrob(main_ann)

gtable_select = function (x,
						  ...
						  ) {
  matches = c(...)
  x$layout = x$layout[matches, , drop=FALSE]
  x$grobs = x$grobs[matches]
  x
}

panels = grepl(pattern="panel", g2$layout$name)
strips = grepl(pattern="strip-t", g2$layout$name)
g2$layout$t[panels] = g2$layout$t[panels]-1
g2$layout$b[panels] = g2$layout$b[panels]-1
new_strips = gtable_select(g2, panels | strips)


### combine heatmap and top annotations
gtable_stack = function(g1, g2){
	g1$grobs = c(g1$grobs, g2$grobs)
	g1$layout = transform(g1$layout, z=z-max(z), name="g2")
	g1$layout = rbind(g1$layout, g2$layout)
	g1
}
temp_heatmap = gtable_stack(g1, new_strips) 


png(paste0(output, "heatmap_test.png"), width=10, 
    height=6, units="in", res=600)
grid.draw(temp_heatmap)
dev.off()

png(paste0(output, "labels_heatmap_test.png"), 
    width=10, height=6, units="in", res=600)
print(main_ann)
dev.off()


