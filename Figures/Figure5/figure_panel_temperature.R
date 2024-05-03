figure_panel_temperature <- function(Fig1tar, Fig2tar, Fig3tar, Fig4tar, Fig5tar, Fig6tar) { 
  Figure <- plot_grid(Fig1tar, Fig2tar, Fig3tar, Fig4tar, Fig5tar, Fig6tar, # The plots we want to put in the grid.
                      ncol = 2, nrow = 3, # The number of columns and rows we want in our grid, respectively.
                      rel_heights = c(1,1,1.105),
                      rel_widths = c(1, 0.9)) 

  return(Figure)
  
}
