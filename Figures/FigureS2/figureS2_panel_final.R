figureS2_panel_final <- function(FigureS2Atar, FigureS2Btar) { 
  Figure <- plot_grid(FigureS2Atar, FigureS2Btar, # The plots we want to put in the grid.
                      ncol = 1, nrow = 2, # The number of columns and rows we want in our grid, respectively
                      rel_heights = c(1.2,1), labels = "AUTO")    
  return(Figure)
  
}