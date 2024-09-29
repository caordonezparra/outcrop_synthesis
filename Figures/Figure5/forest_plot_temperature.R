forest_plot_temperature <- function(data, temperature) {
  data %>%
    dplyr::filter(Temperature == {{temperature}}) %>%
    ggplot(aes(y = Mods,
               x = post_mean,
               xmin = Lower,
               xmax = Upper,
               col = Significant)) +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.8) +
    geom_point(size = 1.8, shape = 15) +
    geom_errorbarh(height=.1) +
    theme_bw(base_size = 8) +
    scale_colour_manual(values = scale_color_temperature) +
    xlim(-12, 9) +
    theme(axis.text.y = if (temperature %in% c(15, 30, 40)) element_blank() else element_text(),
          axis.ticks.y = if (temperature %in% c(15, 30, 40)) element_blank() else element_line(),
          axis.text.x = if (temperature %in% c(10, 15, 20, 30)) element_blank() else element_text(),
          axis.ticks.x = if (temperature %in% c(10, 15, 20, 30)) element_blank() else element_line(),
          axis.title.x = if (temperature %in% c(10, 15, 20, 30)) element_blank() else element_text(face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_blank(), 
          legend.position = "none", 
          axis.title.y = element_blank(),
          axis.text = element_text(color = "black")) +
    labs(title = paste({{temperature}}, "°C"), x = "Effect size") +
    facet_wrap(~Variable)
}

