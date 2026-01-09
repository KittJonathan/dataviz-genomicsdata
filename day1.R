# Day 1

# Load packages ----

library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(datasauRus)
library(gt)
library(MetBrewer)
library(colorspace)
library(patchwork)
library(plotly)

theme_set(theme_minimal())

# Why do we need data visualisations? ----

# Anscombe's quartet

anscombe_long <- anscombe |>
  rowid_to_column(var = "observation") |>
  pivot_longer(
    cols = -observation,
    names_to = c(".value", "set"),
    names_pattern = "(.)(.)"
  )

anscombe_long |>
  ggplot(aes(x, y)) +
  geom_point() +
  facet_wrap(~ set, labeller = labeller(set = \(x) paste("Dataset", x)))

anscombe_long |>
  ggplot(aes(x, y)) +
  geom_point(color = "#28A87D", size = 3.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, color = "#1E8B68", linewidth = 1.2) +
  facet_wrap(~ set, labeller = labeller(set = \(x) paste("Dataset", x)))

anscombe_long |>
  summarise(mean_x = mean(x),
            mean_y = mean(y),
            var_x = var(x),
            var_y = var(y),
            cor = cor(x, y),
            .by = set)

# DatasauRus dozen

datasaurus_dozen |>
  ggplot(aes(x, y, color = dataset)) +
  geom_point(size = 1.5, alpha = 0.7, show.legend = FALSE) +
  scale_color_manual(values = rep(c("#28A87D", "#FFD166", "#7BC5A8"), length.out = 13)) +
  facet_wrap(~dataset)

# Calculate summary statistics for each dataset
summary_stats <- datasaurus_dozen |>
  group_by(dataset) |>
  summarise(
    Mean_X = mean(x),
    Mean_Y = mean(y),
    SD_X = sd(x),
    SD_Y = sd(y),
    Correlation = cor(x, y),
    .groups = "drop"
  ) |>
  arrange(dataset)

# Create beautiful gt table
summary_stats |>
  gt() |>
  tab_header(
    title = "Summary Statistics Across All 13 Datasets",
    subtitle = "Notice how similar the values are despite vastly different shapes"
  ) |>
  cols_label(
    dataset = "Dataset",
    Mean_X = "Mean X",
    Mean_Y = "Mean Y",
    SD_X = "SD X",
    SD_Y = "SD Y",
    Correlation = "Correlation"
  ) |>
  fmt_number(
    columns = c(Mean_X, Mean_Y, SD_X, SD_Y, Correlation),
    decimals = 2
  ) |>
  data_color(
    columns = c(Correlation),
    colors = scales::col_numeric(
      palette = c("#E8F5F1", "#28A87D"),
      domain = NULL
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_text(size = px(12)),
    locations = cells_body()
  ) |>
  tab_options(
    table.font.size = px(14),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(14),
    table.width = pct(90),
    table.align = "center"
  )

# Big data in biology ----

# Color palettes in R

# RColorBrewer
# viridis
# ggsci

# MetBrewer
# wesanderson
# paletteer
# colorspace

# genefilter
# ComplexHeatmap

# https://projects.susielu.com/viz-palette
# ColorBrewer: "colorblind safe" option
# colorspace::deutan() simulator

# Color palette fundamentals

# Define a custom palette
primary_colors <- c("#28A87D", "#A83D28", "#2853a8")
colorspace::demoplot(primary_colors, type = "pie")
colorspace::demoplot(primary_colors, type = "bar")

# Create gradient from color to neutral
seq_colors <- c("#28A87D", "grey98")
sequential_pal <- colorRampPalette(seq_colors)(35)
colorspace::demoplot(sequential_pal, type = "map")
colorspace::demoplot(sequential_pal, type = "heatmap")

# Create diverging palette through neutral center
div_colors <- c("#7754BF", "grey98", "#208462")
diverging_pal <- colorRampPalette(div_colors)(35)
colorspace::demoplot(diverging_pal, type = "map")
colorspace::demoplot(diverging_pal, type = "heatmap")

# Darken primary color for emphasis
dark_colors <- c(colorspace::darken("#28A87D", .35), "grey98")
adjusted_pal <- colorRampPalette(dark_colors)(35)
colorspace::demoplot(adjusted_pal, type = "map")
colorspace::demoplot(adjusted_pal, type = "heatmap")

# Compare color spaces
hcl_adjusted <- colorspace::darken("#28A87D", .35, space = "HCL")
hls_adjusted <- colorspace::darken("#28A87D", .35, space = "HLS")
colorspace::demoplot(hcl_adjusted, type = "map")
colorspace::demoplot(hls_adjusted, type = "map")

# Introduction to ggplot2 ----

set.seed(123)
demo_data <- data.frame(
  x = rnorm(100, 50, 10),
  y = rnorm(100, 50, 10),
  group = sample(c("A", "B"), 100, replace = TRUE)
)

cross_colors <- met.brewer("Cross", 5)

ggplot(demo_data, aes(x = x, y = y, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = c(cross_colors[1], cross_colors[3])) +
  labs(title = "Built with Layers") +
  theme(legend.position = "bottom")

# Building a plot

ggplot(data = iris,
       aes(x = Sepal.Length,
           y = Sepal.Width))

ggplot(data = iris,
       aes(x = Sepal.Length,
           y = Sepal.Width)) +
  geom_point()

ggplot(data = iris,
       aes(x = Sepal.Length,
           y = Sepal.Width,
           colour = Species)) +
  geom_point(size = 3, alpha = 0.7)

ggplot(data = iris,
       aes(x = Sepal.Length,
           y = Sepal.Width,
           colour = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_colour_manual(values = cross_colors) +
  labs(
    title = "Iris Sepal Measurements",
    x = "Sepal Length (mm)",
    y = "Sepal Width (mm)"
  )

# Create custom diverging palette
library(colorspace)
div_colors <- c("#7754BF", "grey98", "#208462")
my_palette <- colorRampPalette(div_colors)(3)

ggplot(iris,
       aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = my_palette) +  # Custom colors!
  labs(title = "Iris Sepal Measurements",
       x = "Sepal Length (cm)",
       y = "Sepal Width (cm)")

div_colors <- c("#7754BF", "grey98", "#208462")
my_palette <- colorRampPalette(div_colors)(3)

ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = my_palette) +
  labs(
    title = "Step 5: Custom Diverging Palette",
    subtitle = "Created with colorRampPalette() - purple to grey to green",
    x = "Sepal Length (cm)",
    y = "Sepal Width (cm)"
  ) +
  theme(legend.position = "right")

cross_colors <- met.brewer("Cross", 3)

# Steps 2-5: Build the plot layer by layer
ggplot(iris,
       aes(x = Sepal.Length,        # Step 1: X variable
           y = Sepal.Width,          # Step 1: Y variable
           color = Species)) +       # Step 3: Color by group
  geom_point(size = 3, alpha = 0.7) +  # Step 2: Add points
  scale_color_manual(values = cross_colors) +  # Step 4: Custom colors
  labs(
    title = "Iris Sepal Measurements",
    subtitle = "Comparing three species",
    x = "Sepal Length (cm)",
    y = "Sepal Width (cm)"
  ) +
  theme(legend.position = "bottom")

# Common geometries

set.seed(456)
demo <- data.frame(
  category = rep(c("A", "B", "C"), each = 30),
  value = c(rnorm(30, 10, 2), rnorm(30, 15, 3), rnorm(30, 12, 2)),
  x = rep(1:30, 3),
  y = c(cumsum(rnorm(30, 0.1, 1)), cumsum(rnorm(30, 0.15, 1)), cumsum(rnorm(30, 0.12, 1)))
)

cross_colors <- met.brewer("Cross", 3)

p1 <- ggplot(demo, aes(x = category, y = value, fill = category)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = cross_colors) +
  labs(title = "geom_boxplot()", subtitle = "Distributions") +
  theme(legend.position = "none")

p2 <- ggplot(demo, aes(x = x, y = y, color = category)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = cross_colors) +
  labs(title = "geom_line()", subtitle = "Trends over time") +
  theme(legend.position = "none")

p3 <- ggplot(demo |> group_by(category) |> summarise(mean_val = mean(value)),
             aes(x = category, y = mean_val, fill = category)) +
  geom_col(alpha = 0.7) +
  scale_fill_manual(values = cross_colors) +
  labs(title = "geom_col()", subtitle = "Comparisons") +
  theme(legend.position = "none")

p4 <- ggplot(demo, aes(x = value, y = x, color = category)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = cross_colors) +
  labs(title = "geom_point()", subtitle = "Relationships") +
  theme(legend.position = "none")

(p1 | p2) / (p3 | p4)

# Our custom theme

custom_theme <- function(base_size = 12,
                         title_size_rel = 1.5,
                         subtitle_size_rel = 1.3,
                         axis_title_size_rel = 1.25) {

  # Start with theme_bw and Asap font
  theme_bw(base_size = base_size, base_family = "Asap") +
    theme(
      # Center and size titles
      plot.title = element_text(hjust = 0.5,
                                size = rel(title_size_rel)),
      plot.subtitle = element_text(hjust = 0.5,
                                   size = rel(subtitle_size_rel)),

      # Axis customization
      axis.title.x = element_text(size = rel(axis_title_size_rel),
                                  margin = margin(t = 12.5)),
      axis.title.y = element_text(size = rel(axis_title_size_rel),
                                  margin = margin(r = 12.5)),

      # Grid and margins
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(linetype = "dashed",
                                        linewidth = 0.3,
                                        colour = "grey80"),
      plot.margin = unit(c(2, 2, 1, 1), "cm")
    )
}

# Set as default theme for all plots
theme_set(custom_theme())

ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, color = Species)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = met.brewer("Cross", 3)) +
  labs(
    title = "Using Custom Theme",
    subtitle = "Consistent, professional styling",
    x = "Sepal Length (cm)",
    y = "Sepal Width (cm)"
  ) +
  theme(legend.position = "bottom")

# Quick exercise: your first ggplot

data(mtcars)

mtcars |>
  ggplot() +
  geom_point(aes(x = wt, y = mpg,
                 colour = as.factor(cyl), size = hp)) +
  geom_smooth(method = "lm",aes(x = wt, y = mpg),
              se = FALSE, show.legend = FALSE,
              colour = "lightgrey") +
  scale_colour_manual(values = cross_colors) +
  labs(title = "mtcars dataset",
       x = "Weight (lbs)",
       y = "Miles / Gallon",
       colour = "Number of cylinders",
       size = "Horse power")

# Common types of data visualization in biological sciences ----

## Scatter plots

# Generate PCA-like data
set.seed(42)
pca_data <- data.frame(
  PC1 = rnorm(200, mean = 0, sd = 3),
  PC2 = rnorm(200, mean = 0, sd = 2),
  group = sample(c("Control", "Treatment"),
                 200, replace = TRUE)
)

cross_colors <- met.brewer("Cross", 5)

ggplot(pca_data, aes(x = PC1, y = PC2,
                     color = group)) +
  geom_point(size = 3.5, alpha = 0.65) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "gray50") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "gray50") +
  stat_ellipse(level = 0.95, linewidth = 1.2) +
  scale_color_manual(
    values = c(cross_colors[1], cross_colors[3]),
    name = ""
  ) +
  labs(title = "PCA Plot - Genomic Features",
       x = "PC1 (44.3%)", y = "PC2 (28.1%)") +
  theme(legend.position = "bottom")

# Exercise: Create your own scatter plot

# Simulate gene expression vs copy number data
set.seed(123)
gene_data <- data.frame(
  gene_id = paste0("Gene_", 1:100),
  expression = rnorm(100, mean = 50, sd = 15),
  copy_number = rnorm(100, mean = 2, sd = 0.5),
  mutation_status = sample(c("WT", "Mutant", "Deletion"),
                           100, replace = TRUE,
                           prob = c(0.6, 0.3, 0.1))
)

# Add correlation between expression and copy number
gene_data$expression <- gene_data$expression +
  (gene_data$copy_number - 2) * 20 + rnorm(100, 0, 5)

head(gene_data)

gene_data |>
  ggplot(aes(x = expression, y = copy_number)) +
  geom_point(aes(colour = mutation_status)) +
  scale_colour_manual(values = cross_colors) +
  geom_smooth(method = "lm", se = TRUE, colour = "red",
              level = 0.95, linewidth = 0.5) +
  labs(title = "Gene expression vs copy number",
       x = "Gene expression",
       y = "Copy number",
       colour = "Mutation status")

## Bar chart

# Step 1: Generate data
set.seed(123)
dose_data <- data.frame(
  dose = rep(c("0.5", "1", "2"), each = 6),
  supplement = rep(c("OJ", "VC"), times = 9),
  length = c(
    rnorm(3, 13, 2), rnorm(3, 8, 1.5),
    rnorm(3, 23, 2), rnorm(3, 17, 2),
    rnorm(3, 26, 2), rnorm(3, 27, 2)
  )
)

# Step 2: Calculate summary statistics
dose_summary <- dose_data |>
  group_by(dose, supplement) |>
  summarise(mean_len = mean(length),
            sd_len = sd(length), .groups = "drop")

# Step 3: Create plot
cross_colors <- met.brewer("Cross", 5)

ggplot(dose_summary,
       aes(x = dose, y = mean_len, fill = supplement)) +
  geom_col(position = position_dodge(0.9),
           width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_len - sd_len,
        ymax = mean_len + sd_len),
    position = position_dodge(0.9), width = 0.25
  ) +
  scale_fill_manual(
    values = c("OJ" = cross_colors[2],
               "VC" = cross_colors[4]),
    labels = c("Orange Juice", "Vitamin C")
  ) +
  labs(x = "Dose (mg/day)", y = "Length", fill = "") +
  theme(legend.position = "top")

# Exercise: RNA-seq bar chart

# Simulate RNA-seq data
set.seed(456)
rnaseq_data <- expand.grid(
  tissue = c("Brain", "Liver", "Heart", "Kidney"),
  gene = c("Gene_A", "Gene_B", "Gene_C"),
  replicate = 1:4
) |>
  mutate(
    expression = case_when(
      gene == "Gene_A" & tissue == "Brain" ~ rnorm(n(), 150, 20),
      gene == "Gene_A" ~ rnorm(n(), 50, 10),
      gene == "Gene_B" & tissue == "Liver" ~ rnorm(n(), 200, 25),
      gene == "Gene_B" ~ rnorm(n(), 60, 15),
      gene == "Gene_C" & tissue == "Heart" ~ rnorm(n(), 180, 22),
      TRUE ~ rnorm(n(), 55, 12)
    )
  )

# Step 2: Calculate summary statistics
rnaseq_summary <- rnaseq_data |>
  group_by(tissue, gene) |>
  summarise(mean_expr = mean(expression),
            sd_expr = sd(expression),
            .groups = "drop")

# Step 3: Create plot
cross_colors <- met.brewer("Cross", 5)

ggplot(rnaseq_summary,
       aes(x = tissue, y = mean_expr, fill = gene)) +
  geom_col(position = position_dodge(0.9),
           width = 0.8) +
  geom_errorbar(
    aes(ymin = mean_expr - sd_expr,
        ymax = mean_expr + sd_expr),
    position = position_dodge(0.9), width = 0.25
  ) +
  scale_fill_manual(
    values = c("Gene_A" = cross_colors[2],
               "Gene_B" = cross_colors[3],
               "Gene_C" = cross_colors[4]),
    labels = c("Gene A", "Gene B", "Gene C")
  ) +
  labs(x = "Tissue", y = "Expression", fill = "",
       title = "Gene expression") +
  theme(legend.position = "top")

## Box plots and violin plots

# Step 1: Generate distributions
set.seed(42)
dist_data <- data.frame(
  type = rep(c("Normal", "Log-transformed"), each = 200),
  values = c(rnorm(200, 80, 15),
             exp(rnorm(200, 4, 0.3)))
)

# Step 2: Create combined plot
cross_colors <- met.brewer("Cross", 5)

ggplot(dist_data, aes(x = type, y = values, fill = type)) +
  geom_violin(alpha = 0.6,
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width = 0.25, alpha = 0.8,
               outlier.shape = 16, outlier.size = 2) +
  scale_fill_manual(values = c(cross_colors[2],
                               cross_colors[4])) +
  labs(title = "Distribution Comparison",
       x = "", y = "Values") +
  theme(legend.position = "none")

# Exercise: Gene Expression distributions

# Simulate gene expression data
set.seed(789)
cell_types <- c("T_cells", "B_cells", "NK_cells",
                "Monocytes", "Dendritic")

expression_data <- data.frame(
  cell_type = rep(cell_types, each = 50),
  expression = c(
    rnorm(50, 100, 15),  # T cells
    rnorm(50, 80, 20),   # B cells - higher variance
    rgamma(50, 5, 0.1),  # NK cells - skewed
    rnorm(50, 120, 10),  # Monocytes
    rlnorm(50, 4, 0.5)   # Dendritic - log-normal
  ),
  condition = rep(c("Healthy", "Disease"), length.out = 250)
)

cross_colors <- met.brewer("Cross", 5)

# Calculate overall mean expression
overall_mean <- mean(expression_data$expression)

# Perform statistical tests (t-tests comparing each cell type between conditions)
stat_results <- expression_data %>%
  group_by(cell_type) %>%
  summarise(
    p_value = t.test(expression[condition == "Healthy"],
                     expression[condition == "Disease"])$p.value,
    max_expr = max(expression)
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = max_expr * 1.1  # Position stars above the data
  )

ggplot(expression_data,
       aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_violin(position = position_dodge(width = 0.9),
              alpha = 0.6) +
  geom_boxplot(position = position_dodge(width = 0.9),
               width = 0.2, alpha = 0.8, show.legend = FALSE) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1, show.legend = FALSE) +
  geom_hline(yintercept = overall_mean, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  geom_text(data = stat_results,
            aes(x = cell_type, y = y_pos, label = significance),
            inherit.aes = FALSE, size = 6, fontface = "bold") +
  scale_fill_manual(values = cross_colors) +
  labs(title = "Distribution comparison",
       x = "Condition", y = "Expression", fill = "Cell type") +
  facet_wrap(~condition) +
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

## Stacked bar plots

# Generate simulated microbiome data
set.seed(42)
columns <- c("Sample", paste0("OTU", 1:11))
samples <- paste0("Day", 1:12)
data <- matrix(runif(12 * 11, min = 0, max = 100), nrow = 12, ncol = 11)
simulated_df <- data.frame(Sample = samples, data)
names(simulated_df) <- columns

# Transform to long format and calculate percentages
simulated_long_df <- simulated_df %>%
  pivot_longer(cols = -Sample, names_to = "OTU", values_to = "Value") %>%
  group_by(Sample) %>%
  mutate(Percentage = Value / sum(Value) * 100) %>%
  ungroup() %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

# Get Cross palette for OTUs
color_palette <- met.brewer("Cross", 11)

# Create stacked bar plot
ggplot(simulated_long_df, aes(x = Sample, y = Percentage, fill = OTU)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Taxonomic Composition Across Time Points",
    subtitle = "Relative abundance of 11 OTUs over 12 sampling days",
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "OTU"
  ) +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Same data, different geom
ggplot(simulated_long_df, aes(x = Sample, y = Percentage, fill = OTU, group = OTU)) +
  geom_area(position = "stack", alpha = 0.8) +
  labs(
    title = "Community Dynamics Over Time",
    subtitle = "Stacked area plot emphasizes continuity and trends",
    x = "Sample",
    y = "Relative Abundance (%)",
    fill = "OTU"
  ) +
  scale_fill_manual(values = color_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Exercise: Microbiome Stacked Plots

# Simulate microbiome data across body sites
set.seed(202)
body_sites <- c("Gut", "Skin", "Oral")
timepoints <- paste0("Week_", 1:8)
taxa <- c("Firmicutes", "Bacteroidetes", "Proteobacteria",
          "Actinobacteria", "Verrucomicrobia", "Other")

microbiome_data <- expand.grid(
  site = body_sites,
  timepoint = timepoints,
  taxon = taxa
) |>
  mutate(
    abundance = case_when(
      taxon == "Firmicutes" & site == "Gut" ~
        runif(n(), 40, 60),
      taxon == "Bacteroidetes" & site == "Gut" ~
        runif(n(), 20, 35),
      taxon == "Proteobacteria" & site == "Skin" ~
        runif(n(), 30, 50),
      TRUE ~ runif(n(), 5, 20)
    )
  ) |>
  group_by(site, timepoint) |>
  mutate(percentage = abundance / sum(abundance) * 100) |>
  ungroup()

# 4. Reorder taxa by average abundance
# We calculate mean abundance first to reorder the factor levels
ordered_taxa <- microbiome_data %>%
  group_by(taxon) %>%
  summarize(mean_abund = mean(percentage)) %>%
  arrange(desc(mean_abund)) %>%
  pull(taxon)

microbiome_data$taxon <- factor(microbiome_data$taxon, levels = ordered_taxa)

# Prepare labels for dominant taxa (>20%)
microbiome_data <- microbiome_data |>
  mutate(label = ifelse(percentage > 20, paste0(round(percentage, 1), "%"), ""))

# Get Cross palette for OTUs
color_palette <- met.brewer("Cross", 6)

# Create stacked bar plot
p_bar <- ggplot(microbiome_data, aes(x = as.numeric(as.factor(timepoint)), y = percentage, fill = taxon)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  labs(
    title = "Microbiome Composition Across Time Points",
    subtitle = "Relative abundance of 6 taxa over 8 weeks",
    x = "Week",
    y = "Relative Abundance (%)",
    fill = "Taxa"
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~site) +
  theme(axis.text.x = element_text(hjust = 0.5),
        legend.position = "top")

# Create stacked area plot
p_area <- ggplot(microbiome_data, aes(x = as.numeric(as.factor(timepoint)), y = percentage, fill = taxon)) +
  geom_area(alpha = 0.8, size = 0.5, colour = "white") +
  labs(
    title = "Microbiome Composition Across Time Points",
    subtitle = "Relative abundance of 6 taxa over 8 weeks",
    x = "Week",
    y = "Relative Abundance (%)",
    fill = "Taxa"
  ) +
  scale_fill_manual(values = color_palette) +
  facet_wrap(~site) +
  theme(axis.text.x = element_text(hjust = 0.5),
        legend.position = "top")

ggplotly(p_bar)
ggplotly(p_area)

# Save all simulated datasets for exercises
save(
  gene_data,           # Scatter plot exercise
  rnaseq_data,         # Bar chart exercise
  expression_data,     # Box/violin exercise
  # expr_matrix,         # Heatmap exercise
  # sample_info,         # Heatmap annotations
  microbiome_data,     # Stacked bar/area exercise
  file = "genomics_viz_exercise_data.RData"
)

# To use in exercises:
load("genomics_viz_exercise_data.RData")
