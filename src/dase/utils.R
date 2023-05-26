
######################################################################################
### Custom col and elements names for raw counts/PSI values data
######################################################################################

custom_names = function(x) {

  for (pos in 1:ncol(x)) {
    
    colnames(x)[pos] = sub(":", "_", as.character(colnames(x)[pos]))
    colnames(x)[pos] = sub("-", "_", as.character(colnames(x)[pos]))
    colnames(x)[pos] = sub(";", "_", as.character(colnames(x)[pos]))

    if (pos == 1) {
      x$feature_id = gsub(":", "_", as.character(x$feature_id))
      x$feature_id = gsub("-", "_", as.character(x$feature_id))
      x$feature_id = gsub(";;", "_", as.character(x$feature_id))
      x$feature_id = gsub(";", "_", as.character(x$feature_id))
    } else if (pos == 2) {
      x$gene_id = gsub(":", "_", as.character(x$gene_id))
      x$gene_id = gsub("-", "_", as.character(x$gene_id))
      x$gene_id = gsub(";;", "_", as.character(x$gene_id))
      x$gene_id = gsub(";", "_", as.character(x$gene_id))
    }
  }
  return(x)
}

#####################################################################################
### Custom visualization themes
######################################################################################

custom_themes = theme_bw(base_size=20, base_family="",
                            base_line_size=1, base_rect_size=1.5) +
  theme(legend.title = element_text(color="black", size=14), 
        legend.text = element_text(size=14),
        axis.text = element_text(color="black", size=14), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        title = element_text(size=14, vjust=0.1),
        plot.title = element_text(hjust=0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

#####################################################################################
### Color sampling
######################################################################################

palette = c("#003147", "#D69C4E", "#9CCB86", "#F2AD00",
            "#E32977", "#46AEA0", "#ECCBAE", "#F08F6E",
            "#046C9A", "#FF0000", "#CCCCCC", "#FC4E2A",
            "#053061", "#B7E6A5", "#5BBCD6", "#ABDDDE", 
            "#E9E29C", "#00A08A", "#3C93C2", "#F98400", 
            "#7C1D6F", "#FABF7B", "#000000", "#1770AB", 
            "#F7F7F7", "#D1E5F0", "#F9D8E6", "#CDE5D2",
            "#EEB479", "#F3A583", "#FDDBC7", "#E1F2E3",
            "#053061", "#D12959", "#6CBA7D", "#E4F1F7",
            "#C40F5B")
             
palette_sampling = function(palette, ncolors) {

  color_list = list()
  
  for (pos in seq(1, ncolors, by=1)) {
    color = sample(palette, 1)
    color_list[[length(color_list) + 1]] = as.character(color)
  }

  colors = unlist(color_list)
  return(colors)
}


