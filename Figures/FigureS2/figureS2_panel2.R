figureS2_panel2 <- function(FigureS2B_1tar, FigureS2B_2tar, FigureS2B_3tar) { 
  Figure <- plot_grid(FigureS2B_1tar, FigureS2B_2tar, FigureS2B_3tar, # The plots we want to put in the grid.
                      ncol = 1, nrow = 3,
                      rel_heights = c(0.8, 0.8, 1.5)) # The number of columns and rows we want in our grid, respectively.
  
  return(Figure)
  
}