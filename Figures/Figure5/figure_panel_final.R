figure_panel_final <- function(Fig1tar, Fig2tar) { 
  Figure <- plot_grid(Fig1tar, Fig2tar, # The plots we want to put in the grid.
                      nrow = 2, # The number of rows we want in our grid.
                      rel_heights = c(4,1),
                      labels = "AUTO") 
  
  return(Figure)
  
}
