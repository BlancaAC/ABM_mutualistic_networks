library(tidyverse)

files <- list.files(path = "Model/Data/Temporal_slices_individuals/", pattern = "network_plot_.*\\.csv", full.names = TRUE)

# Extract plot and slice info from filenames
file.info <- tibble(file = files) %>%
  mutate(
    plot = str_extract(file, "plot_[A-F]"),
    slice = str_extract(file, "slice_[1-3]")
  )

species.id <- read.csv("Model/Data/coords_plot_month.csv", sep=";")

aggregate.to.species <- function(file, plant_lookup) {
  
  df <- read.csv(file, row.names = 1, sep=";")
  
  df_long <- df %>%
    rownames_to_column("Plant_id") %>%
    pivot_longer(-Plant_id, names_to = "Pollinator", values_to = "visits") %>%
    filter(visits > 0)
  
  # Join species info
  df_joined <- df_long %>%
    left_join(plant_lookup, by = "Plant_id")
  
  # Aggregate to species level
  df_species <- df_joined %>%
    group_by(Plant_sp, Pollinator) %>%
    summarise(visits = sum(visits), .groups = "drop")
  
  return(df_species)
}


results <- file.info %>%
  mutate(data = map(file, ~aggregate.to.species(.x, species.id)))


results <- results %>%
  mutate(data = map(data, ~ .x %>%
                      mutate(interaction = paste(Plant_sp, Pollinator, sep = "__"))))

# test differences in interaciton compoisiton across slices within plots

plot_dfs <- results %>%
  unnest(data) %>%
  group_by(plot, slice, interaction) %>%
  summarise(visits = sum(visits), .groups = "drop") %>%
  pivot_wider(names_from = interaction, values_from = visits, values_fill = 0)

plot_list <- plot_dfs %>%
  group_split(plot)

names(plot_list) <- unique(plot_dfs$plot)




library(vegan)

meta <- plot_dfs %>%
  mutate(
    slice = factor(slice),
    plot = factor(plot),
    sample = paste(plot, slice, sep = "_")
  ) %>%
  dplyr::select(sample, plot, slice)

comm <- plot_dfs %>%
  dplyr::select(-plot, -slice) %>%
  as.data.frame()
rownames(comm) <- meta$sample


dist_mat <- vegdist(comm, method = "bray")

nmds <- metaMDS(comm, distance = "bray", k = 2)

plot(nmds$points, col = meta$slice, pch = 19)
legend("topright", legend = levels(meta$slice), col = 1:3, pch = 19)

adonis2(comm ~ slice,
        data = meta,
        method = "bray",
        permutations = how(blocks = meta$plot))


# test differences in pollinator composition across slices within plots

pollinator_dfs <- results %>%
  unnest(data) %>%
  group_by(plot, slice, Pollinator) %>%
  summarise(visits = sum(visits), .groups = "drop") %>%
  pivot_wider(names_from = Pollinator,
              values_from = visits,
              values_fill = 0)

meta <- pollinator_dfs %>%
  mutate(
    slice = factor(slice),
    plot = factor(plot),
    sample = paste(plot, slice, sep = "_")
  ) %>%
  dplyr::select(sample, plot, slice)

comm_poll <- pollinator_dfs %>%
  dplyr::select(-plot, -slice) %>%
  as.data.frame()

rownames(comm_poll) <- meta$sample

adonis2(comm_poll ~ slice,
        data = meta,
        method = "bray",
        permutations = how(blocks = meta$plot))


# test differences in flowering composition across slices within plots

flowering.data <- read.csv("Model/Data/synchrony_data_plant_sp.csv", sep=",")
flowering.data %<>% 
  dplyr::select(Plot, Week, Plant_sp, N_flowers) %>%
  mutate(
    slice = case_when(
      Week >= 1  & Week <= 6  ~ "slice_1",
      Week >= 7  & Week <= 12 ~ "slice_2",
      Week >= 13 & Week <= 19 ~ "slice_3"
    )
  ) %>%
  group_by(Plot, slice, Plant_sp) %>%
  summarise(N_flowers = sum(N_flowers), .groups = "drop")

flowering_comm <- flowering.data %>%
  pivot_wider(names_from = Plant_sp,
              values_from = N_flowers,
              values_fill = 0)
  
meta <- flowering_comm %>%
  mutate(
    slice = factor(slice),
    plot = factor(Plot),
    sample = paste(Plot, slice, sep = "_")
  ) %>%
  dplyr::select(sample, plot, slice)
  
comm <- flowering_comm %>%
  dplyr::select(-Plot, -slice) %>%
  as.data.frame()

rownames(comm) <- meta$sample

adonis2(comm ~ slice,
        data = meta,
        method = "bray",
        permutations = how(blocks = meta$plot))
