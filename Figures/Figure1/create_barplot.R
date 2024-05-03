create_barplot <- function(data, y_var, ylab_text) {
  ggplot(data = data, aes(x = dummy, y = {{y_var}}, fill = Vegetation)) +
    geom_bar(stat = "identity", position = "fill") +
    coord_flip() +
    theme_classic(base_size = 10) +
    scale_y_continuous(expand = c(0,0), labels = function(x) paste0(x*100)) +
    scale_x_discrete(expand = c(0,0)) +
    theme(legend.position = ifelse(ylab_text == "Studies Percentage (%)", "none", "bottom"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.x = element_text(face = "bold")) +
    ylab(ylab_text) +
    scale_fill_manual(values = vegetationColors)
}