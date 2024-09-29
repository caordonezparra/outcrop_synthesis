joined_datasets <- function(traits_joined) {
  selected_families_traits <- traits_joined %>%
    filter(Family %in% selected_families) %>%
    mutate(Family2 = Family)

  traits_all <- traits_joined %>%
    mutate(Family2 = "All") %>%
    bind_rows(selected_families_traits)
  
  return(traits_all)
}