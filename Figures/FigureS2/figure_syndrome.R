figure_syndrome <- function(traits_all) {
  # Process the traits_all dataset
  syndrome_database <- traits_all %>%
    select(Species_acceptedLCVP, Family2, Dispersal_syndrome) %>% # Select these two columns.
    filter(!is.na(Species_acceptedLCVP), # Filter out species with no species record.
           !is.na(Dispersal_syndrome)) %>% # Filter out species without information on dispersal syndrome.
    distinct() %>%
    group_by(Family2, Dispersal_syndrome) %>%
    summarise(counts = n(), .groups = "drop") # Count occurrences
  
  # Create the plot
  Figure <- syndrome_database %>%
    ggplot(aes(x = Family2, y = counts, fill = Dispersal_syndrome)) +
    geom_bar(position = "fill", stat = "identity") +
    scale_fill_manual(values = color_syndrome) +
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