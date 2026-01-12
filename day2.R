# Day 2

# Load packages ----

library(tidyverse)
library(MetBrewer)
library(umap)
library(Rtsne)
library(ggrepel)
library(EnhancedVolcano)
library(ggvolc)

# Set up ----

# Define a custom theme function for ggplot2 to enhance plot aesthetics
custom_theme <- function(base_size = 12, title_size_rel = 1.5,
                         subtitle_size_rel = 1.25,
                         axis_title_size_rel = 1.25,
                         axis_text_size_rel = 1.25,
                         strip_text_size_rel = 1.1) {
  # Use the theme_bw base theme with a specified font family ("Asap")
  theme_bw(base_size = base_size, base_family = "Asap") +
    theme(
      # Customize the title: center-align and set relative size
      plot.title = element_text(hjust = 0.5, size = rel(title_size_rel)),
      # Customize the subtitle: center-align and set relative size
      plot.subtitle = element_text(hjust = 0.5, size = rel(subtitle_size_rel)),
      # Customize the x-axis title: set size and add a top margin
      axis.title.x = element_text(size = rel(axis_title_size_rel),
                                  margin = margin(t = 12.5)),
      # Customize the y-axis title: set size and add a right margin
      axis.title.y = element_text(size = rel(axis_title_size_rel),
                                  margin = margin(r = 12.5)),
      # Customize the x-axis text: set size
      axis.text.x = element_text(size = rel(axis_text_size_rel)),
      # Customize the y-axis text: set size
      axis.text.y = element_text(size = rel(axis_text_size_rel)),
      # Customize facet strip text size
      strip.text = element_text(size = rel(strip_text_size_rel)),
      # Add margins around the plot t,r,b,l
      plot.margin = unit(c(2, 2, 1, 1), "cm"),
      # Remove major grid lines on the x-axis
      panel.grid.major.x = element_blank(),
      panel.grid.major.y  = element_blank(),
      # Remove minor grid lines
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_line(
        linetype = "dashed",
        linewidth = 0.3,
        colour = "grey80"
      ),      # Remove background from facet strips
      strip.background = element_blank()
    )
}

theme_set(custom_theme())

# PCA ----

# Create simulated data with 3 groups
set.seed(123)
pca_demo <- data.frame(
  PC1 = c(rnorm(100, -3, 1), rnorm(100, 0, 1), rnorm(100, 3, 1)),
  PC2 = c(rnorm(100, 2, 1), rnorm(100, 0, 1), rnorm(100, -2, 1)),
  Group = rep(c("Group A", "Group B", "Group C"), each = 100)
)

cross_colors <- met.brewer("Cross", 5)

# Visualise with ellipses
ggplot(pca_demo, aes(x = PC1, y = PC2, colour = Group)) +
  geom_point(size = 2.5, alpha = 0.6) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  scale_colour_manual(values = cross_colors[c(1, 3, 5)]) +
  labs(title = "PCA: Linear Separation",
       x = "PC1 (45.2%)",
       y = "PC2 (23.8%)") +
  theme(legend.position = "bottom")

# UMAP ----

# Simulate UMAP-like structure
set.seed(456)
theta <- seq(0, 2*pi, length.out = 100)
umap_demo <- data.frame(
  UMAP1 = c(2*cos(theta) + rnorm(100, 0, 0.3),
            5 + 1.5*cos(theta) + rnorm(100, 0, 0.3),
            2.5 + 1.5*sin(theta*1.5) + rnorm(100, 0, 0.3)),
  UMAP2 = c(2*sin(theta) + rnorm(100, 0, 0.3),
            1.5*sin(theta) + rnorm(100, 0, 0.3),
            5 + 1.5*cos(theta*1.5) + rnorm(100, 0, 0.3)),
  Group = rep(c("Group A", "Group B", "Group C"), each = 100)
)

# Visualize
ggplot(umap_demo, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = cross_colors[c(1, 3, 5)]) +
  labs(title = "UMAP: Preserves Local & Global Structure")

# Simulate t-SNE-like structure
set.seed(789)
tsne_demo <- data.frame(
  tSNE1 = c(rnorm(100, -5, 0.5), rnorm(100, 0, 0.5), rnorm(100, 5, 0.5)),
  tSNE2 = c(rnorm(100, 3, 0.5), rnorm(100, -3, 0.5), rnorm(100, 0, 0.5)),
  Group = rep(c("Group A", "Group B", "Group C"), each = 100)
)

# Visualize
ggplot(tsne_demo, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.6) +
  scale_color_manual(values = cross_colors[c(1, 3, 5)]) +
  labs(title = "t-SNE: Emphasizes Local Clusters",
       x = "t-SNE1", y = "t-SNE2")

# PCA with prcomp() ----

# Simulate gene expression matrix (samples × genes)
set.seed(42)
n_samples <- 60
n_genes <- 500

# Create expression with tissue-specific patterns
expr_data <- matrix(
  rnorm(n_samples * n_genes, mean = 50, sd = 20),
  nrow = n_samples,
  ncol = n_genes
)

# Add tissue-specific variation
expr_data[1:20, 1:200] <- expr_data[1:20, 1:200] + 30
expr_data[21:40, 201:400] <- expr_data[21:40, 201:400] + 30
expr_data[41:60, 401:500] <- expr_data[41:60, 401:500] + 30

sample_info <- data.frame(
  sample_id = paste0("Sample_", 1:n_samples),
  tissue = rep(c("Brain", "Liver", "Heart"), each = 20)
)

pca_result <- prcomp(expr_data, scale. = TRUE)

pca_scores <- pca_result$x |>
  as.data.frame() |>
  bind_cols(sample_info)

var_exp <- summary(pca_result)$importance[2,]

ggplot(pca_scores,
       aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse() +
  scale_color_manual(values = cross_colors[c(1, 3, 5)]) +
  labs(
    title = "PCA of Gene Expression",
    x = paste0("PC1 (", round(var_exp[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_exp[2]*100, 1), "%)")
  ) +
  theme(legend.position = "bottom")

# Exercise 1: PCA on Gene Expression Data ----

# 1. Create gene expression matrix (samples × genes)
set.seed(123)
n_samples <- 40
n_genes <- 1000

gene_expr <- matrix(
  rpois(n_samples * n_genes, lambda = 50),
  nrow = n_samples,
  ncol = n_genes
)

# Add column names (gene IDs)
colnames(gene_expr) <- paste0("Gene_", 1:n_genes)

# Create sample metadata
sample_meta <- data.frame(
  sample_id = paste0("Sample_", 1:n_samples),
  condition = rep(c("Control", "Treatment"), each = 20),
  batch = rep(1:4, times = 10)
)

# 2. Perform PCA with scaling
pca_result <- prcomp(gene_expr, scale. = TRUE)
summary(pca_result)

# 3. Extract scores and combine with metadata
pca_df <- pca_result$x |>
  as.data.frame() |>
  bind_cols(sample_meta)

# 4. Create PCA plot colored by condition
ggplot(pca_df, aes(x = PC1, y = PC2,
                   color = condition,
                   shape = factor(batch))) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(color = condition)) +
  scale_color_manual(values = c("Control" = "steelblue",
                                "Treatment" = "coral")) +
  labs(
    title = "PCA of RNA-seq Data",
    subtitle = "Check for batch effects and treatment separation"
  )

# Exercise 2 - Scree plot and gene loadings ----

# 1. Create scree plot
variance_explained <- summary(pca_result)$importance[2,]

scree_data <- data.frame(
  PC = paste0("PC", 1:length(variance_explained)),
  Variance = variance_explained * 100
)

ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_col(fill = "steelblue") +
  geom_line(aes(group = 1), color = "#d1422f", linewidth = 1) +
  geom_point(color = "#d1422f", size = 3) +
  labs(title = "Scree Plot: Variance Explained")

# 2. Find top contributing genes to PC1
loadings_pc1 <- pca_result$rotation[, 1]
top_genes <- sort(abs(loadings_pc1), decreasing = TRUE)[1:10]
print("Top 10 genes driving PC1:")
print(top_genes)

# 3. Create loading plot for PC1 vs PC2
loadings_df <- pca_result$rotation[, 1:2] |>
  as.data.frame() |>
  mutate(gene = rownames(pca_result$rotation),
         total_loading = sqrt(PC1^2 + PC2^2)) |>
  arrange(desc(total_loading)) |>
  slice_head(n = 20)

ggplot(loadings_df, aes(x = PC1, y = PC2, label = gene)) +
  geom_segment(aes(xend = 0, yend = 0),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_text(size = 3, hjust = -0.1) +
  labs(title = "Gene Loadings: PC1 vs PC2")

# Exercise 3 - Compare PCA, UMAP and t-SNE ----

# 1. Scale the gene expression data
scaled_expr <- scale(gene_expr)

# 2. PCA
pca_res <- prcomp(scaled_expr)
pca_df <- data.frame(
  Dim1 = pca_res$x[,1],
  Dim2 = pca_res$x[,2],
  Method = "PCA",
  condition = sample_meta$condition,
  batch = factor(sample_meta$batch)
)

# 3. UMAP
set.seed(123)
umap_res <- umap(scaled_expr, n_neighbors = 15)
umap_df <- data.frame(
  Dim1 = umap_res$layout[,1],
  Dim2 = umap_res$layout[,2],
  Method = "UMAP",
  condition = sample_meta$condition,
  batch = factor(sample_meta$batch)
)

# 4. t-SNE
set.seed(123)
tsne_res <- Rtsne(scaled_expr, perplexity = 10, check_duplicates = FALSE)
tsne_df <- data.frame(
  Dim1 = tsne_res$Y[,1],
  Dim2 = tsne_res$Y[,2],
  Method = "t-SNE",
  condition = sample_meta$condition,
  batch = factor(sample_meta$batch)
)

# 5. Combine and plot
all_methods <- bind_rows(pca_df, umap_df, tsne_df)

ggplot(all_methods, aes(x = Dim1, y = Dim2,
                        color = condition, shape = batch)) +
  geom_point(size = 3, alpha = 0.7) +
  facet_wrap(~Method, scales = "free") +
  scale_color_manual(values = c("Control" = "steelblue",
                                "Treatment" = "coral")) +
  labs(title = "Comparison: PCA vs UMAP vs t-SNE") +
  theme(legend.position = "bottom")

# Exercise 4 - Single-cell RNA-seq simulation ----

# 1. Simulate single-cell RNA-seq data
set.seed(456)
n_cells <- 500
n_genes <- 2000

scrna_expr <- matrix(
  rnorm(n_cells * n_genes, mean = 0, sd = 1),
  nrow = n_cells,
  ncol = n_genes
)

# Add cell type-specific expression patterns
scrna_expr[1:150, 1:500] <- scrna_expr[1:150, 1:500] + 3      # T cells
scrna_expr[151:300, 501:1000] <- scrna_expr[151:300, 501:1000] + 3  # B cells
scrna_expr[301:450, 1001:1500] <- scrna_expr[301:450, 1001:1500] + 3  # Monocytes
scrna_expr[451:500, 1501:2000] <- scrna_expr[451:500, 1501:2000] + 3  # NK cells

# Cell metadata
cell_meta <- data.frame(
  cell_id = paste0("Cell_", 1:n_cells),
  cell_type = c(rep("T cell", 150),
                rep("B cell", 150),
                rep("Monocyte", 150),
                rep("NK cell", 50)),
  donor = rep(paste0("Donor_", 1:5), each = 100)
)

# 2. Perform PCA
pca_sc <- prcomp(scrna_expr, scale. = TRUE)

# 3. Extract PC scores
pca_sc_df <- pca_sc$x |>
  as.data.frame() |>
  bind_cols(cell_meta)

var_explained <- summary(pca_sc)$importance[2,]

# 4. Plot PCA colored by cell type
ggplot(pca_sc_df, aes(x = PC1, y = PC2, color = cell_type)) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_manual(values = cross_colors[c(1, 2, 3, 5)]) +
  labs(
    title = "PCA of Single-Cell RNA-seq",
    x = paste0("PC1 (", round(var_explained[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2]*100, 1), "%)"),
    color = "Cell Type"
  )

# 5. BONUS: Check for donor effects
ggplot(pca_sc_df, aes(x = PC1, y = PC2, color = donor)) +
  geom_point(size = 2, alpha = 0.6) +
  labs(title = "Check for Donor Batch Effects")

# Volcano plots ----

# Simulate differential expression data
set.seed(123)
n_genes <- 5000

volcano_demo <- data.frame(
  gene = paste0("Gene_", 1:n_genes),
  log2FC = rnorm(n_genes, mean = 0, sd = 2),
  pvalue = rbeta(n_genes, 0.1, 1)
) |>
  mutate(log10p = -log10(pvalue))

head(volcano_demo)

# Create basic plot (no colors yet)
ggplot(volcano_demo, aes(x = log2FC, y = log10p)) +
  geom_point(alpha = 0.5, size = 1.5, color = "gray50") +
  labs(title = "Basic Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)")

# Add significance classification
volcano_demo <- volcano_demo |>
  mutate(
    significance = case_when(
      abs(log2FC) > 1 & pvalue < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

# Color by significance
ggplot(volcano_demo, aes(x = log2FC, y = log10p, color = significance)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Not significant" = "gray70",
                                "Significant" = "#d1422f")) +
  labs(title = "Volcano Plot with Colors",
       x = "log2 Fold Change", y = "-log10(p-value)")

# Classify as Up, Down, or Not significant
volcano_demo <- volcano_demo |>
  mutate(
    status = case_when(
      log2FC > 1 & pvalue < 0.05 ~ "Upregulated",
      log2FC < -1 & pvalue < 0.05 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Color by direction
ggplot(volcano_demo, aes(x = log2FC, y = log10p, color = status)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#d1422f",
                                "Downregulated" = "#1a5b5b",
                                "Not significant" = "gray70")) +
  labs(title = "Up vs Down Regulation", x = "log2 Fold Change", y = "-log10(p-value)")

# Add threshold lines
ggplot(volcano_demo, aes(x = log2FC, y = log10p, color = status)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Upregulated" = "#d1422f",
                                "Downregulated" = "#1a5b5b",
                                "Not significant" = "gray70")) +
  labs(title = "Volcano Plot with Thresholds",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = NULL) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "bottom")

# Select top 10 genes to label
top_genes <- volcano_demo |>
  filter(status != "Not significant") |>
  arrange(pvalue) |>
  slice_head(n = 10)

# Add labels with ggrepel
ggplot(volcano_demo, aes(x = log2FC, y = log10p, color = status)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_point(data = top_genes, size = 3, alpha = 1) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  scale_color_manual(values = c("Upregulated" = "#d1422f",
                                "Downregulated" = "#1a5b5b",
                                "Not significant" = "gray70")) +
  labs(title = "Volcano Plot with Gene Labels", x = "log2 Fold Change", y = "-log10(p-value)")

# Method 2 - EnhancedVolcano ----

# Create enhanced volcano plot
EnhancedVolcano(volcano_demo,
                lab = volcano_demo$gene,
                x = 'log2FC',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'EnhancedVolcano Example',
                pointSize = 2,
                labSize = 3)

# Method 3 - ggvolc ----

data(all_genes)
head(all_genes, 3)

data(attention_genes)
head(attention_genes, 3)

# Basic plot without segments
ggvolc(all_genes,
       add_seg = FALSE)

# Add significance segments
ggvolc(all_genes,
       add_seg = TRUE) +
  labs(title = "ggvolc with Significance Segments")

# Highlight specific genes
ggvolc(all_genes,
       attention_genes,  # Genes to highlight
       add_seg = TRUE)

# Customize with ggplot2 layers
ggvolc(all_genes, attention_genes, add_seg = TRUE) +
  labs(title = "Highlighted Genes of Interest",
       subtitle = "Using ggvolc with attention_genes") +
  theme_minimal(base_size = 14)

# Exercise 1 - Your first volcano plot ----

# Create data
set.seed(789)
my_deg <- data.frame(
  gene = paste0("Gene_", 1:1000),
  log2FC = rnorm(1000, mean = 0, sd = 2),
  pvalue = rbeta(1000, 0.1, 1)
) |>
  mutate(
    padj = p.adjust(pvalue, method = "BH"),
    log10p = -log10(padj),
    status = case_when(
      log2FC > 1.5 & padj < 0.01 ~ "Upregulated",
      log2FC < -1.5 & padj < 0.01 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Create plot
ggplot(my_deg, aes(x = log2FC, y = log10p, color = status)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "#d1422f",
                                "Downregulated" = "#1a5b5b",
                                "Not significant" = "gray70")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed") +
  labs(title = "My First Volcano Plot",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)")

# Count genes
table(my_deg$status)

# Exercise 2 - Multi-level classification ----

# Multi-level classification
my_deg <- my_deg |>
  mutate(
    detailed = case_when(
      log2FC > 2 & padj < 0.001 ~ "Highly up",
      log2FC > 1 & padj < 0.05 ~ "Moderately up",
      log2FC < -2 & padj < 0.001 ~ "Highly down",
      log2FC < -1 & padj < 0.05 ~ "Moderately down",
      TRUE ~ "Not significant"
    )
  )

# Custom colors
my_colors <- c("Highly up" = "#d73027",
               "Moderately up" = "#fc8d59",
               "Not significant" = "#cccccc",
               "Moderately down" = "#91bfdb",
               "Highly down" = "#4575b4")

# Plot
ggplot(my_deg, aes(x = log2FC, y = log10p, color = detailed)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = my_colors) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.001)),
             linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(-2, -1, 1, 2), linetype = "dashed", alpha = 0.5) +
  labs(title = "Multi-Level Classification",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)")
