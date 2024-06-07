### The following codes use starCounts of Akil-5 as example
library(data.table)
library(mashr)
library(readxl)
library(dplyr)
library(readr)
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(tidyr)

### Preprocess the raw dataset
data <- read.csv("/Users/cathy/Desktop/starCounts_Akil-5.csv", header= TRUE)
head(data)
colnames(data) <- as.character(data[1, ])
data <- data[-1, ]
data <- data[, -c(2:6)]
sample_info <- read_delim("/Users/cathy/Desktop/samples.txt", delim = "\t")
colnames(data)[2:ncol(data)] <- sample_info$sample
data[, 2:19] <- lapply(data[, 2:19], function(x) as.numeric(as.character(x)))

if(any(is.na(data[, 2:19])) | any(data[, 2:19] < 0)) {
  stop("There are NA or negative values in the data.")
}

### Normalize data
counts <- as.matrix(data[, -1]) 
rownames(counts) <- data$Geneid
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge, method="TMM")
dge$samples$norm.factors

### Identify and correct for unknown batch effects
# Is it necessary to Perform Surrogate Variable Analysis (SVA)？？？？？？？

### Differential Expression Analysis using edgeR
group <- as.factor(sample_info$group)
design <- model.matrix(~ group)
y <- estimateDisp(dge, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
top_genes <- topTags(lrt, n = Inf, adjust.method = "BH", p.value = 0.05)
significant_genes <- top_genes$table
view(significant_genes)

### Volcano plot
volcano_data <- data.frame(logFC = lrt$table$logFC, 
                           PValue = lrt$table$PValue,
                           Geneid = rownames(lrt$table))

volcano_data <- volcano_data %>%
  mutate(Significance = ifelse(Geneid %in% rownames(significant_genes), "Significant", "Not Significant"))

ggplot(volcano_data, aes(x = logFC, y = -log10(PValue), color = Significance)) +
  geom_point() +
  scale_color_manual(values = c("grey", "blue")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log Fold Change",
       y = "-log10(PValue)") +
  theme(legend.title = element_blank())
