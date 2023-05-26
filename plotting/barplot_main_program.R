
######################################################################################
### CUSTOM BAR PLOT MAIN
######################################################################################

custom_bar_plot = function(x, colors) { 
  
  bar_plot = ggplot(x, aes(
                        x=annotation,
                        y=percent, 
                        fill=as_event_type)
                        ) +
    geom_bar(stat="identity", color="black", linewidth=1) + 
    theme_classic(base_size=20, base_family="",
                  base_line_size=1.5, base_rect_size=1.5) +
    theme(title = element_text(size=20),
          legend.title = element_text(size=20),
          legend.text = element_text(size=20),
          axis.text.x = element_text(color="black", size=20),
          axis.text.y = element_text(color="black", size=20),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_text(size=20), 
          strip.background = element_blank()) + 
  scale_fill_manual(name="Event type", values=c(colors)) + 
  scale_y_continuous(labels = function(x) paste0(x, '%'))
  
  return(bar_plot)
}


