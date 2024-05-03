figure_panel_alternate <- function(Fig1tar, Fig2tar) { 
  Figure <- plot_grid(Fig1tar, Fig2tar, # The plots we want to put in the grid.
                      ncol = 2, # The number of columns we want in our grid.
                      rel_widths = c(1, 0.9)) 
  
  return(Figure)
  
}
