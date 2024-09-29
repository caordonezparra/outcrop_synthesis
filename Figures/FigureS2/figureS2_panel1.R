figureS2_panel1 <- function(FigureS2A_1tar, FigureS2A_2tar, FigureS2A_3tar, FigureS2A_4tar) { 
  Figure <- plot_grid(FigureS2A_1tar, FigureS2A_2tar, FigureS2A_3tar, FigureS2A_4tar, 
                      ncol = 1, nrow = 4) 
  
  return(Figure)
  
}