library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

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
  ) |>
  arrange(source) |>
  mutate(
    source = factor(source, levels = source)
  )

# A. Number of unique mutations per source
plot_pairs <- ggplot(plot_data, aes(x = source, y = n_unique_pairs, fill = source)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of mutations",
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.title = element_text(face = "bold")
  )

# B. Number of species per source
plot_species <- ggplot(plot_data, aes(x = source, y = n_species, fill = source)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
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

# C. Number of mutations per genus, stacked by source
genus_data <- all_data |>
  group_by(genus, source) |>
  summarise(
    n_unique_pairs = n_distinct(paste(AA_pos_Ecoli, AA_mutation, sep = "_")),
    .groups = "drop"
  ) |>
  arrange(genus, source) |>
  mutate(
    genus = factor(genus, levels = sort(unique(genus))),
    source = factor(source, levels = sort(unique(source)))
  )

plot_genus <- ggplot(genus_data, aes(x = genus, y = n_unique_pairs, fill = source)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  labs(
    x = "Genus",
    y = "Number of mutations",
    fill = "Source"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 40),
    legend.position = "bottom",
    axis.title = element_text(face = "bold")
  )


# Combine panels
# Arrange A and B in one row
top_row <- ggarrange(
  plot_pairs,
  plot_species,
  nrow = 1,
  ncol = 2,
  labels = c("A", "B"),
  hjust = 0,
  vjust = 1
)

# Combine top row with C
overview_barplot <- ggarrange(
  top_row,
  plot_genus,
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

