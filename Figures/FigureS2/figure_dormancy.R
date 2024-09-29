figure_dormancy <- function(traits_all) {
  # Process the traits_all dataset
  dormancy_database <- traits_all %>%
    select(Species_acceptedLCVP, Family2, Dormancy) %>% # Select specific columns
    filter(!is.na(Species_acceptedLCVP), # Remove rows with missing species
           !is.na(Dormancy), # Remove rows with missing dormancy
           Dormancy != "NC") %>% # Remove non-conclusive dormancy records
    distinct() %>% # Remove duplicate rows
    mutate(Record = 1) %>% # Add a column with 1 as value for each row
    pivot_wider(
      names_from = Dormancy, 
      values_from = Record, 
      values_fill = 0 # Fill missing values with 0
    ) %>%
    mutate(Dormancy_PA = case_when(
      ND == 1 & D == 0 ~ "ND", # ND without D
      ND == 0 & D == 1 ~ "D", # D without ND
      ND == 1 & D == 1 ~ "Both" # Both ND and D
    )) %>%
    group_by(Family2, Dormancy_PA) %>%
    summarise(counts = n(), .groups = "drop") # Count occurrences
  
  dormancy_database$Dormancy_PA <- factor(dormancy_database$Dormancy_PA,
                                          levels = c("ND", "D", "Both"))

  # Create the plot
  Figure <- dormancy_database %>%
    ggplot(aes(x = Family2, y = counts, fill = Dormancy_PA)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_dormancy) +
    scale_y_continuous(expand = c(0, 0), labels = percent_format()) +
    theme_linedraw(base_size = 8) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
    ) +
    labs(x = "Family", y = "Percentage of species")
  
  # Return the plot
  return(Figure)
}

