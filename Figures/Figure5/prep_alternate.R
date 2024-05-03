prep_alternate <- function(fdat) {
  constant_temperature <- fdat
  
  constant_temperature$Mods <- factor(constant_temperature$Mods, 
                                      levels = c("Seed mass",
                                                 "Shrub",
                                                 "Herb",
                                                 "Growth form",
                                                 "All"))
  
  constant_temperature$Variable <- factor(constant_temperature$Variable, 
                                          levels = c("Percentage",
                                                     "Time"))
  
  return(constant_temperature)
  
}
