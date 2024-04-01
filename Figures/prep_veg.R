prep_veg <- function(fdat){
  vegetation <- fdat %>% # This function will take the dataset (fdat)...
    mutate(Vegetation = as.factor(Vegetation), # ..., ensure R recognizes the 'Vegetation' column as a factor, ...
           dummy = "dummy") # ... and create a column called 'dummy' where every observation has the value 'dummy'.
  
  vegetation$Vegetation <- ordered(vegetation$Vegetation, # The following lines ensure that vegetation types appear in   the order we want them to appear.
                                   levels = c("Inselberg",
                                              "Campo de altitude",
                                              "Canga",
                                              "Campo rupestre"))
  return(vegetation)
}
