library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)


all_data <- readRDS("processed_data/allData.rds") |>
  mutate(
    genus = stringr::word(species, 1)
  ) 

plot_data <- all_data |>
  group_by(source) |>
  summarise(
    n_unique_pairs = n_distinct(paste(AA_pos_Ecoli, AA_mutation, sep = "_")),
    n_species = n_distinct(species),
    .groups = "drop"
  ) 

source_palette <- c(
  CARD    = "#dfcce3",
  CRyPTIC = "#ffdfba",
  DRAGdb  = "#fbb4ae",
  RifReport = "#bae1ff",
  Yang    = "#cdecc6",
  PointFinder     = "#C8AAAA",  
  "Reference Gene Catalog" = "#DEE791"
)

# A. Number of unique mutations per source
plot_pairs <- ggplot(
  plot_data,
  aes(
    x = reorder(source, n_unique_pairs),
    y = n_unique_pairs,
    fill = source
  )
) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = source_palette) +
  labs(
    x = NULL,
    y = "Number of mutations"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(face = "bold")
  )

# B. Number of species per source
plot_n_species <- ggplot(plot_data, 
                       aes(x = reorder(source, n_species), 
                           y = n_species, 
                           fill = source)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = source_palette) +
  labs(
    x = NULL,
    y = "Number of species",
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(face = "bold")
  )

# C. Number of mutations per species, stacked by source — top 25 species by
# TRUE unique mutation count (deduplicated across sources, matching the
# per-species counts in reports.R). 
mut_species_source <- all_data |>
  mutate(mut_id = paste(AA_pos_Ecoli, AA_mutation, sep = "_")) |>
  distinct(species, mut_id, source)

mut_species_totals <- mut_species_source |>
  group_by(species, mut_id) |>
  summarise(n_sources = n_distinct(source), primary_source = dplyr::first(source), .groups = "drop")

species_data <- mut_species_totals |>
  mutate(source = if_else(n_sources > 1, "Multiple", primary_source)) |>
  count(species, source, name = "n_unique_pairs")

# Identify top 25 species by true unique mutation total
top25_species <- mut_species_totals |>
  count(species, name = "total") |>
  slice_max(total, n = 25) |>
  arrange(desc(total)) |>
  pull(species)

source_palette_c <- c(source_palette, "Multiple" = "#404040")
source_palette_c <- source_palette_c[names(source_palette_c) != "None"]

# Filter and set factor levels (ordered by total, descending)
species_data_top25 <- species_data |>
  filter(species %in% top25_species) |>
  mutate(
    species = factor(species, levels = top25_species),
    source  = factor(source, levels = names(source_palette_c))
  )

plot_species <- ggplot(
  species_data_top25,
  aes(x = species, y = n_unique_pairs, fill = source)
) +
  geom_col(width = 0.7) +
  labs(
    x    = "Species",
    y    = "Number of mutations",
    fill = "Source"
  ) +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = source_palette_c, drop = FALSE) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, face = "italic"),
    plot.margin     = margin(t = 5.5, r = 5.5, b = 5.5, l = 40),
    legend.position = "bottom",
    axis.title      = element_text(face = "bold")
  )


# Combine panels
# Arrange A and B in one row
top_row <- ggarrange(
  plot_pairs,
  plot_n_species,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  hjust = 0,
  vjust = 1
)

# Combine top row with C
overview_barplot <- ggarrange(
  top_row,
  plot_species,
  nrow = 2,
  heights = c(1, 1.4),
  labels = c("", "C"),
  common.legend = TRUE,
  legend = "bottom",
  hjust = 0,
  vjust = 1
)

ggsave(
  filename = "plots/overview_data.pdf",
  plot = overview_barplot,
  width = 12,
  height = 10
)

