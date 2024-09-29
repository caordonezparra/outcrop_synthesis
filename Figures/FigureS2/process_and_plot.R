process_and_plot <- function(traits_all, trait_name) {
  traits_all_transformed <- traits_all %>%
    mutate(
      Transformed_trait = if (trait_name == "Dry_mass") {
        log(!!sym(trait_name))
      } else {
        !!sym(trait_name)
      }
    )
  
  summary_traits <- traits_all_transformed %>%
    group_by(Family2) %>%
    summarise(
      Median = median(Transformed_trait, na.rm = TRUE),
      SD = sd(Transformed_trait, na.rm = TRUE),
      n = n()
    )
  
  largest_families <- summary_traits %>%
    filter(Family2 != "All") %>%
    arrange(desc(Median)) %>%
    slice_head(n = 3)
  
  smallest_families <- summary_traits %>%
    filter(Family2 != "All") %>%
    arrange(Median) %>%
    slice_head(n = 3) %>%
    arrange(desc(Median))
  
  filtered_data <- summary_traits %>%
    filter(Family2 == "All") %>%
    bind_rows(largest_families) %>%
    bind_rows(smallest_families)
  
  filtered_data$Family2 <- factor(
    filtered_data$Family2, 
    levels = c("All", largest_families$Family2, smallest_families$Family2)
  )
  
  traits_all_transformed$Family2 <- factor(
    traits_all_transformed$Family2,
    levels = levels(filtered_data$Family2)
  )
  
  y_label <- switch(trait_name,
                    "Dry_mass" = "Log(Dry mass, mg)",
                    "Water_content" = "Water content (%)",
                    "Embryoless" = "Embryoless seeds (%)",
                    "Viable" = "Viable seeds (%)",
                    paste(trait_name, "(unit)"))
  
  Figure <- traits_all_transformed %>%
    filter(Family2 %in% filtered_data$Family2) %>%
    ggplot(aes(x = Family2, y = Transformed_trait, color = Family2)) +
    geom_jitter(alpha = 0.15, width = 0.2) +
    geom_point(data = filtered_data, aes(x = Family2, y = Median), size = 3) +
    geom_linerange(data = filtered_data, aes(x = Family2, y = Median, ymin = Median - SD, ymax = Median + SD)) +
    theme_linedraw(base_size = 8) +
    scale_color_manual(values = scale_color_families) +
    theme(legend.position = "none", 
          axis.title.x = if (trait_name == "Viable") element_text(face = "bold") else element_blank()) +
    labs(
      x = "Family",
      y = y_label
    )
  
  return(Figure)
}