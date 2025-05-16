# Comprehensive PTM Analysis Project in R
# Author: [Your Name]
# Date: May 1, 2025

# This project demonstrates a comprehensive analysis workflow for post-translational
# modification (PTM) data, with a focus on phosphoproteomics. The workflow includes:
# 1. Data preprocessing and normalization
# 2. Differential PTM site identification
# 3. Functional enrichment analysis
# 4. Motif analysis
# 5. Integration with protein-protein interaction networks
# 6. Visualization of results

# Load required libraries
library(tidyverse)      # For data manipulation and visualization
library(limma)          # For differential analysis
library(PhosR)       # For phosphoproteomics analysis
library(clusterProfiler) # For enrichment analysis
library(ggplot2)        # For visualization
library(pheatmap)       # For heatmaps
library(RColorBrewer)   # For color palettes
library(igraph)         # For network analysis
library(rmotifx)        # For motif analysis
library(org.Hs.eg.db)   # For human gene annotation
library(biomaRt)        # For gene ID conversion
library(ggrepel)        # For improved text labels in plots


# 1. DATA LOADING AND PREPROCESSING -------------------------------------------

# Load phosphoproteomics data
# Sample data format: rows are phosphosites, columns are samples
# Data includes intensity values, site localization, and peptide information
phospho_data <- read.csv("data/phospho_intensity_data.csv", stringsAsFactors = FALSE)
sample_metadata <- read.csv("data/sample_metadata.csv", stringsAsFactors = FALSE)

# Examine data structure
head(phospho_data)
str(phospho_data)
summary(phospho_data)

# Parse phosphosite information
parse_phosphosites <- function(data) {
  # Extract protein names, positions, and amino acids
  data$Protein <- gsub("(.*?)_.*", "\\1", data$Phosphosite)
  data$Position <- as.numeric(gsub(".*_([0-9]+).*", "\\1", data$Phosphosite))
  data$Residue <- gsub(".*_[0-9]+([STY]).*", "\\1", data$Phosphosite)
  
  # Create standardized site identifier (Protein_ResiduePosition)
  data$Site_ID <- paste0(data$Protein, "_", data$Residue, data$Position)
  
  return(data)
}

phospho_data <- parse_phosphosites(phospho_data)

# Filter out low-quality phosphosites
# Remove sites with localization probability < 0.75
filtered_data <- phospho_data %>%
  filter(Localization_prob >= 0.75)

# Log2 transform intensity values
intensity_cols <- grep("intensity_", names(filtered_data), value = TRUE)
filtered_data[intensity_cols] <- log2(filtered_data[intensity_cols])

# Handle missing values (e.g., imputation or filtering)
# Option 1: Remove sites with >50% missing values
missing_threshold <- 0.5
filtered_data <- filtered_data %>%
  mutate(missing_pct = rowSums(is.na(.[intensity_cols])) / length(intensity_cols)) %>%
  filter(missing_pct <= missing_threshold) %>%
  select(-missing_pct)

# Option 2: Impute missing values (minimal imputation example)
impute_minimal <- function(x) {
  if(sum(is.na(x)) / length(x) <= 0.3) {  # Only impute if < 30% missing
    x[is.na(x)] <- min(x, na.rm = TRUE) / 2  # Impute with half the minimum value
  }
  return(x)
}

filtered_data[intensity_cols] <- t(apply(filtered_data[intensity_cols], 1, impute_minimal))

# Normalize data
# Median normalization
normalized_data <- filtered_data
normalized_data[intensity_cols] <- t(apply(filtered_data[intensity_cols], 1, function(x) {
  x - median(x, na.rm = TRUE)
}))

# 2. DIFFERENTIAL PTM ANALYSIS ------------------------------------------------

# Create design matrix for differential analysis
# Example: Treatment vs Control comparison
design <- model.matrix(~ 0 + sample_metadata$Condition)
colnames(design) <- levels(factor(sample_metadata$Condition))

# Create contrast matrix
contrasts <- makeContrasts(
  Treatment_vs_Control = Treatment - Control,
  levels = design
)

# Perform differential analysis using limma
intensity_matrix <- as.matrix(normalized_data[, intensity_cols])
rownames(intensity_matrix) <- normalized_data$Site_ID

fit <- lmFit(intensity_matrix, design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, coef = "Treatment_vs_Control", number = Inf)
results$Site_ID <- rownames(results)

# Merge results with site information
diff_sites <- merge(results, 
                    normalized_data %>% select(Site_ID, Protein, Position, Residue),
                    by = "Site_ID")

# Filter for significantly regulated sites
# adj.P.Val < 0.05 and |logFC| > 1
sig_sites <- diff_sites %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

# Save results
write.csv(diff_sites, "results/differential_phosphosites.csv", row.names = FALSE)
write.csv(sig_sites, "results/significant_phosphosites.csv", row.names = FALSE)

# 3. FUNCTIONAL ENRICHMENT ANALYSIS -------------------------------------------

# Extract significant proteins
sig_proteins <- unique(sig_sites$Protein)

# Map to Entrez Gene IDs
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_map <- getBM(attributes = c("uniprot_swissprot", "entrezgene_id", "external_gene_name"),
                     filters = "uniprot_swissprot",
                     values = sig_proteins,
                     mart = mart)

# Perform GO enrichment analysis
ego <- enrichGO(gene = protein_map$entrezgene_id,
                universe = NULL,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05)

# Simplify GO terms to reduce redundancy
ego_simplified <- simplify(ego)

# Perform KEGG pathway enrichment
ekegg <- enrichKEGG(gene = protein_map$entrezgene_id,
                    organism = "hsa",
                    pvalueCutoff = 0.05)

# Save enrichment results
write.csv(as.data.frame(ego_simplified), "results/GO_enrichment.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "results/KEGG_enrichment.csv", row.names = FALSE)

# 4. MOTIF ANALYSIS ----------------------------------------------------------

# Prepare sequences for motif analysis
# We need the surrounding sequence for each phosphosite

# Load a dataframe with protein sequences
protein_sequences <- read.csv("data/protein_sequences.csv", stringsAsFactors = FALSE)

# Extract sequence windows around phosphosites
extract_motif_window <- function(protein_id, position, residue, window_size = 7) {
  # Get protein sequence
  seq <- protein_sequences$Sequence[protein_sequences$Protein == protein_id]
  
  if(length(seq) == 0) return(NA)
  
  # Calculate window boundaries
  start_pos <- max(1, position - window_size)
  end_pos <- min(nchar(seq), position + window_size)
  
  # Extract window
  window <- substring(seq, start_pos, end_pos)
  
  # Center the window if needed
  if(position - window_size < 1) {
    prefix <- paste(rep("_", abs(position - window_size - 1)), collapse = "")
    window <- paste0(prefix, window)
  }
  
  if(position + window_size > nchar(seq)) {
    suffix <- paste(rep("_", position + window_size - nchar(seq)), collapse = "")
    window <- paste0(window, suffix)
  }
  
  return(window)
}

# Apply to significant sites
sig_sites$Motif <- mapply(extract_motif_window, 
                          sig_sites$Protein, 
                          sig_sites$Position, 
                          sig_sites$Residue)

# Split into upregulated and downregulated sites
up_sites <- sig_sites %>% filter(logFC > 0)
down_sites <- sig_sites %>% filter(logFC < 0)

# Perform motif enrichment using rmotifx
# For upregulated sites
up_motifs <- rmotifx(
  foreground = up_sites$Motif,
  background = diff_sites$Motif,
  central_residue = "S|T|Y",
  min_occurrences = 5,
  min_pct_foreground = 0.05
)

# For downregulated sites
down_motifs <- rmotifx(
  foreground = down_sites$Motif,
  background = diff_sites$Motif,
  central_residue = "S|T|Y",
  min_occurrences = 5,
  min_pct_foreground = 0.05
)

# Save motif results
write.csv(up_motifs, "results/upregulated_motifs.csv", row.names = FALSE)
write.csv(down_motifs, "results/downregulated_motifs.csv", row.names = FALSE)

# 5. INTEGRATION WITH PPI NETWORKS -------------------------------------------

# Load protein-protein interaction data
# Can use STRING or BioGRID data
ppi <- read.csv("data/ppi_data.csv", stringsAsFactors = FALSE)

# Create graph from PPI data
ppi_graph <- graph_from_data_frame(ppi, directed = FALSE)

# Map significant proteins to the network
sig_proteins_in_network <- intersect(sig_proteins, names(V(ppi_graph)))

# Extract subnetwork of significant proteins
sig_subgraph <- induced_subgraph(ppi_graph, sig_proteins_in_network)

# Calculate node metrics
node_degree <- degree(sig_subgraph)
node_betweenness <- betweenness(sig_subgraph)
node_closeness <- closeness(sig_subgraph)

# Create node attribute dataframe
node_attributes <- data.frame(
  Protein = names(node_degree),
  Degree = node_degree,
  Betweenness = node_betweenness,
  Closeness = node_closeness
)

# Identify kinases in the network
# Load kinase annotation data
kinase_data <- read.csv("data/kinase_annotations.csv", stringsAsFactors = FALSE)

# Identify kinases in the significant proteins
sig_kinases <- intersect(sig_proteins, kinase_data$Kinase)

# Map to the network
node_attributes$Is_Kinase <- node_attributes$Protein %in% sig_kinases

# Save network analysis results
write.csv(node_attributes, "results/network_attributes.csv", row.names = FALSE)

# 6. VISUALIZATION ----------------------------------------------------------

# 1. Volcano plot of phosphosites
volcano_plot <- ggplot(diff_sites, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(data = subset(diff_sites, adj.P.Val < 0.01 & abs(logFC) > 2),
                  aes(label = Site_ID), max.overlaps = 20) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 P-value", 
       title = "Volcano Plot of Differential Phosphosites",
       subtitle = "Treatment vs Control")

ggsave("figures/volcano_plot.png", volcano_plot, width = 10, height = 8)

# 2. Heatmap of significant sites
# Prepare data
heatmap_data <- intensity_matrix[rownames(intensity_matrix) %in% sig_sites$Site_ID, ]
heatmap_data <- heatmap_data[order(sig_sites$logFC[match(rownames(heatmap_data), sig_sites$Site_ID)]), ]

# Create annotation
annotation_col <- data.frame(
  Condition = sample_metadata$Condition,
  row.names = colnames(heatmap_data)
)

# Generate heatmap
pheatmap(
  heatmap_data,
  annotation_col = annotation_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Heatmap of Significant Phosphosites",
  filename = "figures/heatmap_phosphosites.png",
  width = 10,
  height = 12
)

# 3. Enrichment plot
dotplot(ego_simplified, showCategory = 15, title = "GO Biological Process Enrichment") +
  theme(axis.text.y = element_text(size = 8))
ggsave("figures/go_enrichment.png", width = 10, height = 8)

# 4. Network visualization
# Prepare network for visualization in igraph
V(sig_subgraph)$size <- node_attributes$Degree * 2
V(sig_subgraph)$color <- ifelse(node_attributes$Is_Kinase, "red", "lightblue")
V(sig_subgraph)$label <- V(sig_subgraph)$name
V(sig_subgraph)$label.cex <- 0.8

# Save network plot
png("figures/ppi_network.png", width = 1200, height = 1000, res = 100)
plot(sig_subgraph, 
     layout = layout_with_fr(sig_subgraph),
     vertex.label.dist = 0.5,
     edge.width = 0.5,
     main = "PPI Network of Significant Proteins")
dev.off()

# 5. Motif visualization
# Use ggseqlogo for motif visualization if available
if (requireNamespace("ggseqlogo", quietly = TRUE)) {
  library(ggseqlogo)
  
  # Prepare motif data
  up_motif_seqs <- up_motifs$motif
  
  # Plot logos
  p <- ggseqlogo(up_motif_seqs) + 
    theme_minimal() +
    ggtitle("Enriched Motifs in Upregulated Phosphosites")
  
  ggsave("figures/motif_logos.png", p, width = 8, height = 6)
}

# 7. SUMMARY REPORT ----------------------------------------------------------

# Generate an HTML report using R Markdown
# This requires the rmarkdown package
if (requireNamespace("rmarkdown", quietly = TRUE)) {
  # Create an R Markdown file
  cat('---
title: "PTM Analysis Report"
author: "[Your Name]"
date: "', format(Sys.Date(), "%B %d, %Y"), '"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(tidyverse)
library(knitr)
library(DT)
```

## Project Overview

This report summarizes the results of a post-translational modification (PTM) analysis
comparing Treatment vs Control conditions.

## Differential Analysis Results

```{r}
sig_sites <- read.csv("results/significant_phosphosites.csv")
DT::datatable(sig_sites)
```

## Functional Enrichment

```{r}
go_results <- read.csv("results/GO_enrichment.csv")
DT::datatable(go_results)
```

## Motif Analysis

```{r}
knitr::include_graphics("figures/motif_logos.png")
```

## Network Analysis

```{r}
knitr::include_graphics("figures/ppi_network.png")
```

## Conclusion

Summary of key findings and biological implications.

', file = "PTM_Analysis_Report.Rmd")
  
  # Render the R Markdown file
  rmarkdown::render("PTM_Analysis_Report.Rmd", output_file = "results/PTM_Analysis_Report.html")
}

# Final message
cat("PTM analysis completed successfully!\n")
cat("Results saved in the 'results' directory.\n")
cat("Figures saved in the 'figures' directory.\n")
cat("Summary report available at 'results/PTM_Analysis_Report.html'.\n")