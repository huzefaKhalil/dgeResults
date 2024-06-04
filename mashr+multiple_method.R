library(RSQLite)
library(mashr)
library(dplyr)
library(ggplot2)

con <- dbConnect(RSQLite::SQLite(), dbname = "/Users/cathy/Desktop/hdrfData.sqlite")
dbListTables(con)
results <- dbGetQuery(con, "SELECT * FROM comparisonData")
print(head(results))
dbDisconnect(con)

###1. Filter genes with low variance and low mean expression:###
###Calculate the variance for each gene if there are multiple logFC measurements per gene, otherwise, skip it.
#if ("logFC" %in% colnames(result)) {
  #gene_variances <- var(result$logFC)
  #print(gene_variances)
#} else {
  #stop("The 'logFC' column is not present in the data.")
#}

#threshold <- 0.1
#gene_variances <- apply(result[, "logFC", drop = FALSE], 1, var) 
#filtered_genes <- result[gene_variances > threshold, ]
#print(head(filtered_genes))

#gene_means <- rowMeans(as.matrix(result$logFC)) 
#mean_threshold <- 0.1  
#filtered_genes <- filtered_genes[gene_means > mean_threshold, ]
#print(head(filtered_genes))

###2. Use of mashr for multivariate analysis:###
### Use observations with name "Akil-1" otherwise there will be too many likelihood matrices
result_Akil_1 <- subset(results, name == "Akil-1")
effects <- as.matrix(result_Akil_1[, "logFC"])  
std_errors <- as.matrix(result_Akil_1[, "se"])  

print("Effect Sizes Matrix:")
print(head(effects))
print("Standard Errors Matrix:")
print(head(std_errors))

data = mash_set_data(effects, std_errors)
U.c = cov_canonical(data)
m.1 = mash(data, U.c)

result = get_significant_results(m.1)# a vector of indices for significant results

significant_genes <- result_Akil_1[result, ]
significant_gene_ids <- significant_genes$id  
print("Significant Gene IDs:")
print(head(significant_gene_ids))

result_Akil_1$significant_mashr <- ifelse(result_Akil_1$id %in% significant_gene_ids, "Yes", "No")

#Volcano plot
ggplot(result_Akil_1, aes(x = logFC, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Gene Expression for Akil-1",
       x = "Log Fold Change",
       y = "-Log10(p-value)",
       color = "Significant")

### Akil-5###
#result_Akil_5 <- subset(results, name == "Akil-5")
#effects <- as.matrix(result_Akil_5[, "logFC"])  
#std_errors <- as.matrix(result_Akil_5[, "se"])  

#data = mash_set_data(effects, std_errors)
#U.c = cov_canonical(data)
#m.1 = mash(data, U.c)

#result = get_significant_results(m.1)# a vector of indices for significant results

#significant_genes <- result_Akil_5[result, ]
#significant_gene_ids <- significant_genes$id  
#print("Significant Gene IDs:")
#print(head(significant_gene_ids))

#result_Akil_1$significant <- ifelse(result_Akil_5$id %in% significant_gene_ids, "Yes", "No")

#Volcano plot
#ggplot(result_Akil_5, aes(x = logFC, y = -log10(pvalue), color = significant)) +
  #geom_point(alpha = 0.5) +
  #scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  #theme_minimal() +
  #labs(title = "Volcano Plot of Gene Expression for Akil-5",
       #x = "Log Fold Change",
       #y = "-Log10(p-value)",
       #color = "Significant")



###3. Apply machine learning techniques to further refine your list of significant genes. ###
### Use Akil-1 as example. ###
library(glmnet)
library(randomForest)
library(caret)
library(stats)
library(ggrepel)

# PCA for Akil-1
X <- as.matrix(result_Akil_1[, c("logFC", "se")])  
y <- as.factor(result_Akil_1$significant)
pca_result <- prcomp(X, scale. = TRUE)
print(summary(pca_result))

pca_scores <- pca_result$x
pca_importance <- summary(pca_result)$importance
top_pcs <- pca_scores[, 1:2]  

# Identify significant genes based on PCA
significant_genes_pca <- result_Akil_1[which(abs(top_pcs[,1]) > 2 | abs(top_pcs[,2]) > 2), ]
significant_gene_ids_pca <- significant_genes_pca$id

print("Significant Gene IDs based on PCA:")
print(head(significant_gene_ids_pca))

# Plot PCA
pca_data <- as.data.frame(pca_result$x)
pca_data$Gene <- rownames(result_Akil_1)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "PCA Plot of Gene Expression",
       x = "Principal Component 1",
       y = "Principal Component 2")

result_Akil_1$significant_pca <- ifelse(result_Akil_1$id %in% significant_gene_ids_pca, "Yes", "No")

# Volcano plot for PCA significant genes
ggplot(result_Akil_1, aes(x = logFC, y = -log10(pvalue), color = significant_pca)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Gene Expression for Akil-1 with PCA",
       x = "Log Fold Change",
       y = "-Log10(p-value)",
       color = "Significant (PCA)")

#########################################################################
set.seed(123)

###Lasso Regression for Akil-1
library(caret)
library(glmnet)
library(ggplot2)

trainIndex <- createDataPartition(y, p = .8, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Lasso Regression with cross-validation
lasso_model <- cv.glmnet(X_train, y_train, family = "binomial", alpha = 1)
best_lambda <- lasso_model$lambda.min
final_lasso_model <- glmnet(X_train, y_train, family = "binomial", alpha = 1, lambda = best_lambda)

selected_features <- coef(final_lasso_model)
nonzero_indices <- which(selected_features != 0)
selected_feature_names <- rownames(selected_features)[nonzero_indices]
selected_feature_values <- selected_features[nonzero_indices]

print("Selected features by Lasso:")
for (i in 1:length(selected_feature_names)) {
  cat(selected_feature_names[i], ":", selected_feature_values[i], "\n")
}

nonzero_indices <- nonzero_indices[-1]
significant_gene_ids_lasso <- rownames(X_train)[nonzero_indices]

print("Significant Gene IDs based on Lasso Regression:")
print(significant_gene_ids_lasso)

result_Akil_1$significant_lasso <- ifelse(result_Akil_1$id %in% significant_gene_ids_lasso, "Yes", "No")

# Volcano plot for Lasso significant genes
ggplot(result_Akil_1, aes(x = logFC, y = -log10(pvalue), color = significant_lasso)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Gene Expression for Akil-1 with Lasso Regression",
       x = "Log Fold Change",
       y = "-Log10(p-value)",
       color = "Significant (Lasso)")

#########################################################################
### UMAP for Akil-1
library(umap)

# UMAP for Akil-1
umap_result <- umap(X)
umap_scores <- umap_result$layout

significant_genes_umap <- result_Akil_1[which(abs(umap_scores[,1]) > 2 | abs(umap_scores[,2]) > 2), ]
significant_gene_ids_umap <- significant_genes_umap$id

print("Significant Gene IDs based on UMAP:")
print(head(significant_gene_ids_umap))

# Plot UMAP
umap_data <- as.data.frame(umap_scores)
umap_data$Gene <- rownames(result_Akil_1)
ggplot(umap_data, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "UMAP Plot of Gene Expression",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2")

result_Akil_1$significant_umap <- ifelse(result_Akil_1$id %in% significant_gene_ids_umap, "Yes", "No")

# Volcano plot for UMAP significant genes
ggplot(result_Akil_1, aes(x = logFC, y = -log10(pvalue), color = significant_umap)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Gene Expression for Akil-1 with UMAP",
       x = "Log Fold Change",
       y = "-Log10(p-value)",
       color = "Significant (UMAP)")

#########################################################################
### ICA for Akil-1
library(fastICA)
ica_result <- fastICA(X, n.comp = 2)
ica_scores <- ica_result$S

significant_genes_ica <- result_Akil_1[which(abs(ica_scores[,1]) > 2 | abs(ica_scores[,2]) > 2), ]
significant_gene_ids_ica <- significant_genes_ica$id

print("Significant Gene IDs based on ICA:")
print(head(significant_gene_ids_ica))

# Plot ICA
ica_data <- as.data.frame(ica_scores)
ica_data$Gene <- rownames(result_Akil_1)
ggplot(ica_data, aes(x = V1, y = V2)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "ICA Plot of Gene Expression",
       x = "Independent Component 1",
       y = "Independent Component 2")

result_Akil_1$significant_ica <- ifelse(result_Akil_1$id %in% significant_gene_ids_ica, "Yes", "No")

# Volcano plot for ICA significant genes
ggplot(result_Akil_1, aes(x = logFC, y = -log10(pvalue), color = significant_ica)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Gene Expression for Akil-1 with ICA",
       x = "Log Fold Change",
       y = "-Log10(p-value)",
       color = "Significant (ICA)")




#################################Summary#################################
###4. Summarize the number of significant genes identified by each method
summary_table <- result_Akil_1 %>%
  summarise(
    significant_pca = sum(significant_pca == "Yes"),
    significant_lasso = sum(significant_lasso == "Yes"),
    significant_umap = sum(significant_umap == "Yes"),
    significant_ica = sum(significant_ica == "Yes")
  )

print(summary_table)
library(reshape2)
summary_table_1 <- melt(summary_table)

colnames(summary_table_1) <- c("Method", "Count")

# Plot the summary table
ggplot(summary_table_1, aes(x = Method, y = Count, fill = Method)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Number of Significant Genes Identified by Each Method",
       x = "Method",
       y = "Count of Significant Genes",
       fill = "Method") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
