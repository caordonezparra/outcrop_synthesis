figure_season <- function(traits_all) {
  # Process the traits_all dataset
  season_database <- traits_all %>%
    select(Family2, Species_acceptedLCVP, Dispersal_period) %>% # Select these two columns.
    filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
           !is.na(Dispersal_period)) %>% # Filter out species without information on dispersal period.
    distinct() %>% # Keep only distinct observations.
    mutate(Record = 1) %>% # Create a new column called 'Record' where every observation has 1 as value.
    pivot_wider(names_from = Dispersal_period, values_from = Record, values_fill = 0) %>% # Reshape the data set to     take the 'Dispersal_period' column and create four new columns (based on 'Dispersal_period'). It will fill each      observation with the value of 'Rercord' or zero if no 'Record' value is available.
    mutate(sum = ED + ER + LD + LR, # Create a new column with the sum of the values of the four columns we have just   created.
           Dispersal_period2 = case_when(sum > 1 ~ "More than 1", # Create a new column, called 'Dispersal_period'              based on the values on the 'sum' column. This will allow us to identify whether a species disperse their             seeds in more than one period.
                                         sum == 1 & LD == 1 ~ "LD",
                                         sum == 1 & ED == 1 ~ "ED",
                                         sum == 1 & ER == 1 ~ "ER",
                                         sum == 1 & LR == 1 ~ "LR"),
           Dispersal_period2 = as.factor(Dispersal_period2)) %>% # Ensure that R understands this new column as a factor. 
    group_by(Family2, Dispersal_period2) %>%
    summarise(counts = n(), .groups = "drop") # Count occurrences
  
  season_database$Dispersal_period2 <- factor(season_database$Dispersal_period2,
                                              levels = c("ED", "LD", "ER", "LR", "More than 1"))
  
  # Create the plot
  Figure <- season_database %>%
    ggplot(aes(x = Family2, y = counts, fill = Dispersal_period2)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_season) +
    scale_y_continuous(expand = c(0, 0), labels = percent_format()) +
    theme_linedraw(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_text(face = "bold")
    ) +
    labs(x = "Family", y = "Percentage of species")
  
  # Return the plot
  return(Figure)
}