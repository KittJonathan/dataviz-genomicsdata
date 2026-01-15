# Day 2

# Load packages ----

library(tidyverse)
library(MetBrewer)
library(umap)
library(Rtsne)
library(ggrepel)
library(EnhancedVolcano)
library(ggvolc)
library(karyoploteR)
library(GenomicRanges)
library(Gviz)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(VennDiagram)
library(UpSetR)
library(RColorBrewer)
library(ggtree)
library(ape)
library(ggtreeExtra)
library(plotly)

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

# Exercise 3 - Faceted comparisons ----

# Create 3 treatment comparisons
set.seed(999)
create_comp <- function(name) {
  data.frame(
    gene = paste0("Gene_", 1:1000),
    log2FC = rnorm(1000, 0, 2),
    pvalue = rbeta(1000, 0.1, 1),
    comparison = name
  ) |>
    mutate(padj = p.adjust(pvalue, "BH"),
           log10p = -log10(padj),
           sig = ifelse(abs(log2FC) > 1 & padj < 0.05, "Sig", "NS"))
}

comp_data <- bind_rows(
  create_comp("Treatment A"),
  create_comp("Treatment B"),
  create_comp("Treatment C")
)

# Faceted plot
ggplot(comp_data, aes(x = log2FC, y = log10p, color = sig)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Sig" = "#d1422f", "NS" = "gray80")) +
  facet_wrap(~comparison, ncol = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  labs(title = "Multiple Treatment Comparisons",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)") +
  theme_minimal()

# Peak visualisation ----

# Simulate genomic signal data
set.seed(123)
genomic_pos <- 1:1000
signal <- abs(rnorm(1000, mean = 2, sd = 3))
signal[450:550] <- signal[450:550] + rnorm(101, mean = 10, sd = 2)

# Create data frame
track_data <- data.frame(
  position = genomic_pos,
  signal = signal
)
head(track_data)

# Plot peak track
ggplot(track_data, aes(x = position, y = signal)) +
  geom_area(fill = "steelblue", alpha = 0.6) +
  geom_line(color = "darkblue", linewidth = 0.5) +
  annotate("rect", xmin = 450, xmax = 550, ymin = 0, ymax = Inf,
           fill = "#d1422f", alpha = 0.1) +
  annotate("text", x = 500, y = 18, label = "Peak", color = "#d1422f", size = 5) +
  labs(title = "Example Peak Track",
       x = "Genomic Position (kb)", y = "Signal Intensity")+
  custom_theme()

# Method 1 - karyoploteR ----

custom.genome <- data.frame(
  chr = "chr1",
  start = 1,
  end = 4600000
)

plotKaryotype(genome = custom.genome)

# karyoploteR: Basic Karyoplot

# Create example peak data
set.seed(123)
peaks <- GRanges(
  seqnames = paste0("chr", sample(1:22, 100, replace = TRUE)),
  ranges = IRanges(start = sample(1:200000000, 100),
                   width = sample(200:2000, 100)),
  score = runif(100, 0, 10)
)

# Create karyoplot and add peaks
kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:22))
kpPlotRegions(kp, data = peaks, col = "steelblue")

# karyoploteR: multi-track

# Create karyoplot
kp <- plotKaryotype(genome = "hg38", chromosomes = paste0("chr", 1:5))

# Track 1: Sample A peaks (top)
kpPlotRegions(kp, data = peaks,
              r0 = 0, r1 = 0.3,
              col = "#d1422f")

# Track 2: Sample B peaks (middle)
kpPlotRegions(kp, data = peaks,
              r0 = 0.35, r1 = 0.65,
              col = "#1a5b5b")

# Track 3: Density (bottom)
kpPlotDensity(kp, data = peaks,
              r0 = 0.7, r1 = 1,
              col = "#f4ab5c")

# r0/r1: vertical position (0=bottom, 1=top)

# karyoploteR: density plot

set.seed(456)
peaks_chr1 <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(
    start = sample(1:200000000, 200),
    width = sample(200:2000, 200)
  )
)

# Create karyoplot
kp <- plotKaryotype(genome = "hg38", chromosomes = "chr1")

# Add peak density with window
kpPlotDensity(kp,
              data = peaks,
              window.size = 1e7,  # 10Mb windows
              col = "steelblue")

# karyoploteR: coverage plot

# Create coverage data
set.seed(789)
coverage <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = seq(1e6, 5e6, by = 1000),
                   width = 1000),
  score = abs(rnorm(4001, mean = 50, sd = 20))
)

# Plot coverage
kp <- plotKaryotype(genome = "hg38", chromosomes = "chr1")
kpPlotCoverage(kp,
               data = coverage,
               col = "forestgreen",
               r0 = 0, r1 = 0.5)

# karyoploteR: chromosome zoom

# Focus on specific region
kp <- plotKaryotype(
  genome = "hg38",
  chromosomes = "chr1",
  plot.type = 4
)

# Add multiple data types
kpPlotRegions(kp, data = peaks, col = "#d1422f", r0 = 0, r1 = 0.45)
kpPlotDensity(kp, data = peaks, col = "#1a5b5b", r0 = 0.55, r1 = 1)

# Add horizontal line at specific position
kpAbline(kp, h = 0.5, col = "gray", lty = 2)

# Method 2: Gviz ----

# Basic tracks

# Create example ChIP-seq peaks
chipseq_peaks <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = c(1000000, 1005000, 1010000),
                   end = c(1002000, 1007000, 1012000)),
  score = c(25, 45, 35)
)

# Create tracks
genome_axis <- GenomeAxisTrack()
ideogram <- IdeogramTrack(genome = "hg38", chromosome = "chr1")
peaks_track <- AnnotationTrack(
  chipseq_peaks, name = "Peaks", fill = "steelblue"
)

# Plot
plotTracks(
  list(ideogram, genome_axis, peaks_track),
  from = 1000000,
  to = 1015000,
  chromosome = "chr1"
)

# Data tracks

# Create coverage data for sample
set.seed(123)
coverage_data <- GRanges(
  seqnames = "chr1",
  ranges = IRanges(start = seq(1e6, 1.5e6, by = 100),
                   width = 100),
  score = abs(rnorm(5001, mean = 50, sd = 20))
)

# Create DataTrack
coverage_track <- DataTrack(
  coverage_data,
  name = "ChIP Signal",
  type = "histogram",
  fill = "darkgreen",
  col = "darkgreen"
)

# Plot with axis
plotTracks(
  list(ideogram, genome_axis, coverage_track),
  from = 1e6,
  to = 1.5e6,
  chromosome = "chr1"
)

# Multi-sample comparison

# Sample 1: Control
sample1_track <- DataTrack(
  coverage_data,
  name = "Control",
  type = "histogram",
  fill = "gray60"
)

# Sample 2: Treatment (higher signal)
sample2_track <- DataTrack(
  coverage_data,
  name = "Treatment",
  type = "histogram",
  fill = "firebrick"
)

# Plot both samples
plotTracks(
  list(ideogram, genome_axis,
       sample1_track,
       sample2_track),
  from = 1e6,
  to = 1.5e6,
  chromosome = "chr1"
)

# Different plot types

# Type: polygon (filled area)
track1 <- DataTrack(coverage_data, name = "Polygon",
                    type = "polygon", fill = "lightblue")

# Type: line
track2 <- DataTrack(coverage_data, name = "Line",
                    type = "l", col = "darkblue", lwd = 2)

# Type: points
track3 <- DataTrack(coverage_data, name = "Points",
                    type = "p", col = "#d1422f")

# Type: smooth (smoothed line)
track4 <- DataTrack(coverage_data, name = "Smooth",
                    type = "smooth", col = "purple")

# Plot all
plotTracks(list(ideogram, genome_axis, track1, track2, track3, track4),
           from = 1e6, to = 1.5e6, chromosome = "chr1")

# Overlay tracks

# Create overlapping data
set.seed(999)
data1 <- GRanges("chr1", IRanges(seq(1e6, 1.5e6, 100), width = 100),
                 score = abs(rnorm(5001, 50, 20)))
data2 <- GRanges("chr1", IRanges(seq(1e6, 1.5e6, 100), width = 100),
                 score = abs(rnorm(5001, 40, 15)))

# Create overlay track
overlay_track <- OverlayTrack(
  trackList = list(
    DataTrack(data1, type = "l", col = "#d1422f", name = "Sample1"),
    DataTrack(data2, type = "l", col = "#1a5b5b", name = "Sample2")
  )
)

# Plot
plotTracks(
  list(ideogram, genome_axis, overlay_track),
  from = 1e6,
  to = 1.5e6,
  chromosome = "chr1"
)

# Advanced visualizations ----

# Basic heatmap

# Simulate gene expression matrix
set.seed(123)
n_genes <- 50
n_samples <- 10

expr_matrix <- matrix(
  rnorm(n_genes * n_samples, mean = 50, sd = 20),
  nrow = n_genes,
  dimnames = list(paste0("Gene_", 1:n_genes),
                  paste0("Sample_", 1:n_samples))
)

# Add tissue-specific patterns
expr_matrix[1:20, 1:5] <- expr_matrix[1:20, 1:5] + 40  # Brain
expr_matrix[21:40, 6:10] <- expr_matrix[21:40, 6:10] + 40  # Liver

# Create heatmap
pheatmap(expr_matrix,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Gene Expression Heatmap")

# Annotated heatmap

# Create sample annotations
sample_anno <- data.frame(
  Tissue = rep(c("Brain", "Liver"), each = 5),
  Treatment = rep(c("Control", "Treated"), times = 5),
  row.names = colnames(expr_matrix)
)

# Define annotation colors
anno_colors <- list(
  Tissue = c(Brain = "#E64B35", Liver = "#4DBBD5"),
  Treatment = c(Control = "gray70", Treated = "gold")
)

# Create annotated heatmap
pheatmap(expr_matrix,
         scale = "row",
         annotation_col = sample_anno,
         annotation_colors = anno_colors,
         show_rownames = FALSE,
         main = "Annotated Gene Expression")

# ComplexHeatmap

# Create color function
col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#1a5b5b", "white", "#d1422f")
)

# Create complex heatmap
Heatmap(scale(expr_matrix),
        name = "Z-score",
        col = col_fun,
        top_annotation = HeatmapAnnotation(
          Tissue = sample_anno$Tissue,
          col = list(Tissue = c(Brain = "#E64B35", Liver = "#4DBBD5"))
        ),
        show_row_names = FALSE,
        column_title = "Complex Heatmap with Annotations")

# Venn diagrams

# Simulate differentially expressed genes
set.seed(456)

# Create gene universe
all_genes <- paste0("Gene_", 1:1000)

# Treatment A DEGs (200 genes)
treatment_A <- sample(all_genes, 200)

# Treatment B DEGs (250 genes, 100 overlap with A)
treatment_B <- c(
  sample(treatment_A, 100),  # 100 overlap
  sample(setdiff(all_genes, treatment_A), 150)  # 150 unique
)

# Treatment C DEGs (180 genes)
treatment_C <- sample(all_genes, 180)

# Create list
gene_lists <- list(
  Treatment_A = treatment_A,
  Treatment_B = treatment_B,
  Treatment_C = treatment_C
)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("Treatment A", "Treatment B", "Treatment C"),
  filename = NULL,
  fill = c("#E64B35", "#4DBBD5", "#00A087"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2
)

# Display
grid.draw(venn.plot)

# UpSet plots

# Simulate 5 different conditions
set.seed(789)
all_genes <- paste0("Gene_", 1:1000)

gene_sets <- list(
  Brain = sample(all_genes, 250),
  Liver = sample(all_genes, 200),
  Heart = sample(all_genes, 180),
  Kidney = sample(all_genes, 220),
  Lung = sample(all_genes, 190)
)

# Create UpSet plot
upset(fromList(gene_sets),
      nsets = 5,
      order.by = "freq")

# UpSet with queries

# Highlight specific intersections
upset(fromList(gene_sets),
      nsets = 5,
      order.by = "freq",
      queries = list(
        list(query = intersects,
             params = list("Brain", "Liver"),
             color = "#d1422f", active = TRUE)
      ))

# Show only specific intersections
upset(fromList(gene_sets),
      nsets = 5,
      nintersects = 20,  # Show top 20 intersections
      mb.ratio = c(0.6, 0.4),  # Adjust matrix/bar ratio
      order.by = "freq",
      text.scale = c(1.5, 1.5, 1.3, 1.3, 1.5, 1.2),
      point.size = 3.5,
      line.size = 1)

# Exercise 1: custom heatmap

# Create cell type-specific expression
set.seed(111)
n_genes <- 100
n_cells <- 30

expr_data <- matrix(
  rnorm(n_genes * n_cells, mean = 5, sd = 2),
  nrow = n_genes,
  dimnames = list(paste0("Gene_", 1:n_genes),
                  paste0("Cell_", 1:n_cells))
)

# Add cell type patterns
expr_data[1:30, 1:10] <- expr_data[1:30, 1:10] + 5   # T cells
expr_data[31:60, 11:20] <- expr_data[31:60, 11:20] + 5  # B cells
expr_data[61:100, 21:30] <- expr_data[61:100, 21:30] + 5  # Monocytes

# Create annotations
cell_anno <- data.frame(
  CellType = c(rep("T cell", 10), rep("B cell", 10), rep("Monocyte", 10)),
  Donor = rep(paste0("Donor", 1:3), length.out = 30),
  row.names = colnames(expr_data)
)

# Plot
library(pheatmap)
pheatmap(expr_data, scale = "row",
         annotation_col = cell_anno,
         show_rownames = FALSE,
         main = "Single-Cell Expression Patterns")

# Exercise 2: correlation heatmap

# Calculate sample correlations
cor_matrix <- cor(expr_matrix)

# Custom colors for correlation
cor_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

pheatmap(cor_matrix,
         color = cor_colors,
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 10,
         main = "Sample Correlation Heatmap")

# Exercise 3: pathway overlap

# Simulate pathway analysis results
set.seed(222)
all_genes <- paste0("Gene_", 1:500)

pathways <- list(
  "Cell Cycle" = sample(all_genes, 80),
  "Apoptosis" = sample(all_genes, 70),
  "Immune Response" = sample(all_genes, 90),
  "Metabolism" = sample(all_genes, 100)
)

# Create Venn for 2 pathways
venn.diagram(
  x = list(CellCycle = pathways[[1]],
           Apoptosis = pathways[[2]]),
  filename = NULL,
  fill = c("lightblue", "pink"),
  alpha = 0.5
)

# For all 4 pathways, use UpSet
upset(fromList(pathways),
      order.by = "freq",
      text.scale = 1.5)

# Exercise 4: Multi-condition comparison

# Simulate DEGs from multiple comparisons
set.seed(333)
all_genes <- paste0("Gene_", 1:800)

deg_sets <- list(
  "Condition1_vs_Control" = sample(all_genes, 150),
  "Condition2_vs_Control" = sample(all_genes, 180),
  "Condition3_vs_Control" = sample(all_genes, 160),
  "Condition4_vs_Control" = sample(all_genes, 140),
  "Condition5_vs_Control" = sample(all_genes, 170)
)

# Create comprehensive UpSet plot
upset(fromList(deg_sets),
      nsets = 5,
      nintersects = 30,
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      text.scale = 1.3,
      point.size = 3,
      main.bar.color = "steelblue")

# Phylogenetic trees ----

# Load or create a tree
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)

# Basic plot
ggtree(tree) +
  geom_tiplab(size = 3) +
  theme_tree2()

# Flipper horizontal
ggtree(tree) +
  coord_flip() +
  geom_tiplab(size = 3)

# Reversed direction
ggtree(tree) +
  coord_flip() +
  scale_x_reverse() +
  geom_tiplab(size = 3)

# Circular layout
set.seed(456)
circ_tree <- rtree(20)

cross_colors <- met.brewer("Cross", 5)

ggtree(circ_tree, layout = "circular") +
  geom_tiplab(size = 2.5, offset = 0.5) +
  labs(title = "Circular Tree")

# Unrooted layout
ggtree(circ_tree, layout = "unrooted") +
  geom_tiplab(size = 2.5) +
  labs(title = "Unrooted Tree")

# Highlighting clades

# Highlight a clade
ggtree(tree, layout = "circular") +
  geom_hilight(node = 23, fill = cross_colors[1], alpha = 0.6) +
  geom_tiplab(size = 3)

# Highlight multiple clades
ggtree(tree, layout = "circular") +
  geom_hilight(node = 5, fill = cross_colors[1], alpha = 0.5, extend = 0.2) +
  geom_hilight(node = 13, fill = cross_colors[3], alpha = 0.5, extend = 0.2) +
  geom_tiplab(size = 3)

set.seed(789)
large_tree <- rtree(50)

cross_colors <- met.brewer("Cross", 5)

# Find some internal nodes for highlighting
p <- ggtree(large_tree, layout = "circular")

# Highlight different clades
p +
  geom_hilight(node = 55, fill = cross_colors[1], alpha = 0.4, extend = 0.15) +
  geom_hilight(node = 70, fill = cross_colors[3], alpha = 0.4, extend = 0.15) +
  geom_hilight(node = 85, fill = cross_colors[5], alpha = 0.4, extend = 0.15) +
  geom_tiplab(size = 2, align = FALSE) +
  labs(title = "Phylogenetic Tree with Highlighted Clades",
       subtitle = "Different colors represent distinct evolutionary groups")

# Adding strips and annotations

ggtree(tree, layout = "circular") +
  geom_strip(
    taxa1 = "A",
    taxa2 = "m",
    color = cross_colors[1],
    barsize = 3,
    offset = 0.05
  ) +
  geom_tiplab(size = 3)

# Fancy tree: multiple annotations

set.seed(101)
fancy_tree <- rtree(40)

# Create some example data for the tips
tip_data <- data.frame(
  taxa = fancy_tree$tip.label,
  group = sample(c("Group A", "Group B", "Group C", "Group D"),
                 40, replace = TRUE),
  value1 = rnorm(40, mean = 50, sd = 15),
  value2 = runif(40, 0, 100)
)

cross_colors <- met.brewer("Cross", 5)

p <- ggtree(fancy_tree, layout = "circular", size = 0.8)

# Add highlighted clades
p2 <- p +
  geom_hilight(node = 45, fill = cross_colors[1], alpha = 0.3, extend = 0.1) +
  geom_hilight(node = 60, fill = cross_colors[2], alpha = 0.3, extend = 0.1) +
  geom_hilight(node = 75, fill = cross_colors[4], alpha = 0.3, extend = 0.1)

# Add tip labels
p3 <- p2 +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.3, offset = 0.5) +
  geom_tippoint(aes(color = tip_data$group[match(label, tip_data$taxa)]),
                size = 2, alpha = 0.8) +
  scale_color_manual(values = cross_colors[1:4], name = "Group")

# Add title
p3 +
  labs(title = "Annotated Phylogenetic Tree",
       subtitle = "With highlighted clades and tip metadata") +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

# Adding external data with ggtreeExtra

set.seed(202)
heat_tree <- rtree(25)

# Create heatmap data
heat_data <- expand.grid(
  taxa = heat_tree$tip.label,
  trait = paste0("Trait", 1:5)
)
heat_data$value <- rnorm(nrow(heat_data), mean = 50, sd = 20)

cross_colors <- met.brewer("Cross", 5)

ggtree(heat_tree, layout = "circular") +
  geom_tiplab(size = 2.5, offset = 1.5) +
  geom_fruit(
    data = heat_data,
    geom = geom_tile,
    mapping = aes(y = taxa, x = trait, fill = value),
    offset = 0.1,
    pwidth = 0.25,
    color = "white",
    size = 0.5
  ) +
  scale_fill_gradientn(
    colors = c(cross_colors[1], "white", cross_colors[3]),
    name = "Value"
  ) +
  labs(title = "Phylogenetic Tree with Trait Heatmap",
       subtitle = "Integrating evolutionary relationships with trait data") +
  theme(legend.position = "right",
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5))

# Interactive trees with ggtree

# Create tree
p <- ggtree(tree) + geom_tiplab(size = 3)
p

# Convert to interactive
plotly(p)
