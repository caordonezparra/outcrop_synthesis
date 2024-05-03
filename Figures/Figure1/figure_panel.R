figure_panel <- function(Fig1tar, Fig2tar) { 
  Figure <- plot_grid(Fig1tar, Fig2tar, # The plots we want to put in the grid.
                       ncol = 1, nrow = 2, # The number of columns and rows we want in our grid, respectively.
                       labels = c("E", NULL), # The labels of our plots. 
                       rel_heights = c(1,1.75)) # The relative height of each plot. Since Figure 1F has the legend, we need to adjust its height so that bars have the same height.
return(Figure)
  
  }


