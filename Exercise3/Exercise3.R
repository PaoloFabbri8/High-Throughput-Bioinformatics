# Load and install necessary libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install and load necessary libraries ------------------------------------

if (!require("limma", quietly = TRUE))
  BiocManager::install("limma")

if (!require("piano", quietly = TRUE))
  BiocManager::install("piano")

if (!require("randomForest", quietly = TRUE))
  install.packages("randomForest")


library(limma)
library(piano)
library(randomForest)
library(ggplot2) # Optional - for plotting


# Set working directory with setwd()
setwd("C:/Users/paolo/Desktop/Università/Magistrale/2° Anno/Secondo Semestre/High-throughput Bioinformatics/Lab/Exercise3/Exercise3")
getwd()


# 1. Drug combination response prediction ---------------------------------

# Read drug combination response data from data/drug_combination_responses.csv
data_drug_combination <- read.csv("data/drug_combination_responses.csv")
data_drug_combination

# Controlling the structure of the data to apply one hot encoding
str(data_drug_combination)

# Hint: factor() and model.matrix() for one-hot encoding
# Hint: randomForest() and predict() for training and predicting with random forest models
# Hint: cor() for computing the Pearson correlations on the test sets
# One hot encoding of drugs, drug concentrations, and cell line IDs
data_drug_combination$drug_1 <- factor(data_drug_combination$drug_1)
data_drug_combination$drug_2 <- factor(data_drug_combination$drug_2)
data_drug_combination$drug_1_concentration <- factor(data_drug_combination$drug_1_concentration)
data_drug_combination$drug_2_concentration <- factor(data_drug_combination$drug_2_concentration)
data_drug_combination$cell_line <- factor(data_drug_combination$cell_line)

# Building the model matrix for the random forest model
X <- model.matrix(response ~ ., data = data_drug_combination)[, -1]
y <- data_drug_combination$response

# Setting the seed for reproducibility
set.seed(123)

n <- nrow(data_drug_combination)
train_idx1 <- sample(1:n, size = 0.8 * n)

X_train1 <- X[train_idx1, ]
y_train1 <- y[train_idx1]

X_test1 <- X[-train_idx1, ]
y_test1 <- y[-train_idx1]

# Train Random Forest
rf1 <- randomForest(x = X_train1,
                    y = y_train1,
                    ntree = 500)

# Prediction
pred1 <- predict(rf1, X_test1)

# Pearson correlation
cor1 <- cor(pred1, y_test1)
cat("Pearson correlation (All data):", cor1, "\n")



cell_lines <- unique(data_drug_combination$cell_line)

train_cells <- sample(cell_lines, size = 0.8 * length(cell_lines))

train_idx2 <- which(data_drug_combination$cell_line %in% train_cells)
test_idx2  <- which(!data_drug_combination$cell_line %in% train_cells)

X_train2 <- X[train_idx2, ]
y_train2 <- y[train_idx2]

X_test2 <- X[test_idx2, ]
y_test2 <- y[test_idx2]

# Train Random Forest
rf2 <- randomForest(x = X_train2,
                    y = y_train2,
                    ntree = 500)

# Prediction
pred2 <- predict(rf2, X_test2)

# Pearson correlation
cor2 <- cor(pred2, y_test2)
cat("Pearson correlation (Cell-line split):", cor2, "\n")

sink("Figures/Pearson_correlation.txt")
cat("Pearson correlation (All data):", cor1, "\n")
cat("Pearson correlation (Cell-line split):", cor2, "\n")
sink()

# 2. Gene expression data exploration -------------------------------------


# Read the gene expression data and metadata from data/gene_expression.csv and
# data/cell_line_metadata.csv
data_gene_expression <- read.csv("data/gene_expression.csv")
metadata_cell_line <- read.csv("data/cell_line_metadata.csv")
data_gene_expression
metadata_cell_line

str(data_gene_expression)
str(metadata_cell_line)

# controlling the missing values
sum(is.na(data_gene_expression))


# Copy of the original Dataframe (just for safety)
data_gene_expression_backup <- data_gene_expression

# Imputation using the mean
for (j in 2:ncol(data_gene_expression)) {
  gene_mean <- mean(data_gene_expression_backup[, j], na.rm = TRUE)
  data_gene_expression_backup[is.na(data_gene_expression_backup[, j]), j] <- gene_mean
}

data_gene_expression <- data_gene_expression_backup

# Control again the missing values after imputation
sum(is.na(data_gene_expression))

# Performing PCA analysis on the gene expression data
# Hint: prcomp() for PCA
gene_expression_matrix <- data_gene_expression[, -1]

# Scale = True standa
pca <- prcomp(gene_expression_matrix)

var_explained <- pca$sdev^2 / sum(pca$sdev^2)
pc1_var <- var_explained[1]
pc2_var <- var_explained[2]

sink("Figures/Explained_variance.txt")
cat("Variance explained by PC1:", pc1_var*100, "%\n")
cat("Variance explained by PC2:", pc2_var*100, "%\n")
sink()

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  cell_line = data_gene_expression$cell_line
)
# Merging the data to create the dataframe with the tissue type information for coloring the points in the PCA plot
pca_df <- merge(pca_df, metadata_cell_line, by = "cell_line")

ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue_type)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of NCI-60 Gene Expression",
    x = paste0("PC1"),
    y = paste0("PC2")
  )

ggsave("Figures/Plot_PCA.png", width = 8, height = 6, dpi = 300)

# 3-5. Differential gene expression analysis ------------------------------


# For performing the differential gene expression analysis using limma, 
# you can refer to the examples shown in Section 9.2 of the limma user guide:
# https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

TP53_status <- factor(metadata_cell_line$TP53_status, levels = c("Wildtype", "Mutated"))
TP53_status

# Design matrix
design <- model.matrix(~ 0 + TP53_status)
colnames(design) <- c("Wildtype", "Mutated")
design

# Defining the contrast matrix for the comparison
contrast_matrix <- makeContrasts(
  TP53_mutated_vs_wildtype = Mutated - Wildtype,
  levels = design
)
contrast_matrix


# Extract expression matrix 
matrix_expression <- t(as.matrix(data_gene_expression[, -1]))

# Fit linear model
fit <- lmFit(matrix_expression, design)

# Apply contrast
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# List of genes that were identified to be statistically significantly differentially expressed
significant_genes <- topTable(fit2, number = Inf, adjust.method = "BH", p.value = 0.05)
significant_genes$Gene <- rownames(significant_genes)
report_genes <- significant_genes[, c("Gene", "logFC", "adj.P.Val")]

head(report_genes)

# Extract all genes with adjusted p-value < 0.05
all_genes <- topTable(fit2, number = Inf, adjust.method = "BH")
all_genes$Gene <- rownames(all_genes)
all_genes$significant <- ifelse(all_genes$adj.P.Val < 0.05, "Significant", "Not Significant")


ggplot(all_genes, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed Genes",
    x = "log2 Fold Change (Mutated vs Wildtype)",
    y = "-log10 Adjusted P-value",
    color = "Gene Status"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")

ggsave("Figures/Volcano_plot.png", width = 8, height = 6, dpi = 300)

write.csv(report_genes, "Figures/Significant_Genes_TP53.csv", row.names = FALSE)

# 6-10. Gene set enrichment analysis --------------------------------------


# Hypergeometric test -----------------------------------------------------

# Load the gene set collection (gsc) from data/gene_set_collection.RData
gene_sets <- readRDS("data/gene_set_collection.rData")
head(gene_sets)
str(gene_sets)
names(gene_sets)
gene_sets$gsc

# Loop over the gene sets and perform hypergeometric test using phyper function,
# see ?phyper for documentation on the input parameters
DE_genes <- report_genes$Gene

background_genes <- rownames(matrix_expression)

actual_gene_sets <- gene_sets$gsc

pvals <- numeric(length(actual_gene_sets))
names(pvals) <- names(actual_gene_sets)

for(gs_name in names(actual_gene_sets)) {
  genes_in_set <- actual_gene_sets[[gs_name]]
  
  m <- sum(background_genes %in% genes_in_set)  # genes of GO term in the background
  n <- length(background_genes) - m              # genes not in the GO term
  K <- length(DE_genes)                          # total DE genes
  k <- sum(DE_genes %in% genes_in_set)           # DE genes in the GO term
  
  pvals[gs_name] <- phyper(k - 1, m, n, K, lower.tail = FALSE)
}

pvals
pvals_sorted <- sort(pvals)
head(pvals_sorted,20)

# Adjusting p-values for multiple testing using Benjamini-Hochberg (BH) method
pvals_adj <- p.adjust(pvals, method = "BH")
pvals_adj_sorted <- sort(pvals_adj)

enrichment_results <- data.frame(
  GO_Term = names(pvals_adj_sorted),
  Adj_P_value = pvals_adj_sorted,
  stringsAsFactors = FALSE
)

significant_GO_terms <- enrichment_results[enrichment_results$Adj_P_value < 0.01, ]
head(significant_GO_terms, 20)
cat("Number of significant GO terms:", nrow(significant_GO_terms), "\n")

write.csv(significant_GO_terms, "Figures/Significant_GO_Terms.csv", row.names = FALSE)


# GSEA --------------------------------------------------------------------

# Hint: runGSA for performing GSEA. Use 1000 permutations and set geneSetStat="gsea" and adjMethod="fdr". 
# Hint: GSAsummaryTable for summarizing results

# Prepare the logFC vector
logFC_vector <- all_genes$logFC
names(logFC_vector) <- all_genes$Gene
logFC_vector <- sort(logFC_vector, decreasing = TRUE)
head(logFC_vector)

# GSEA
gsea_results <- runGSA(
  geneLevelStats = logFC_vector,
  gsc = gene_sets,
  nPerm = 1000,
  geneSetStat = "gsea",
  adjMethod = "fdr",
  verbose = TRUE
)

gsea_summary <- GSAsummaryTable(gsea_results)
head(gsea_summary)

# Genes UP-regulated: p-value "p adj (dist.dir.up)"
up_regulated <- gsea_summary[!is.na(gsea_summary$`p adj (dist.dir.up)`) & 
                               gsea_summary$`p adj (dist.dir.up)` < 0.01, ]
cat("GO terms enriched in UP-regulated genes:", nrow(up_regulated), "\n")
print(up_regulated[, c("Name", "p adj (dist.dir.up)")])

# Genes DOWN-regulated: p-value "p adj (dist.dir.dn)"
down_regulated <- gsea_summary[!is.na(gsea_summary$`p adj (dist.dir.dn)`) & 
                                 gsea_summary$`p adj (dist.dir.dn)` < 0.01, ]
cat("GO terms enriched in DOWN-regulated genes:", nrow(down_regulated), "\n")
print(down_regulated[, c("Name", "p adj (dist.dir.dn)")])

# Save result
write.csv(up_regulated, "Figures/GSEA_UP_GO_Terms.csv", row.names = FALSE)
write.csv(down_regulated, "Figures/GSEA_DOWN_GO_Terms.csv", row.names = FALSE)


# GO terms significative from hypergeometric test
hyper_terms <- significant_GO_terms$GO_Term

# GO terms significative from GSEA (both up and down)
gsea_terms <- unique(c(up_regulated$Name, down_regulated$Name))

# Overlap
overlap <- intersect(hyper_terms, gsea_terms)
cat("GO terms in hypergeometric test:", length(hyper_terms), "\n")
cat("GO terms in GSEA:", length(gsea_terms), "\n")
cat("Overlapping GO terms:", length(overlap), "\n")
print(overlap)


