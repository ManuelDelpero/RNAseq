args <- commandArgs(trailingOnly = TRUE)
counts_files <- strsplit(args[1], " ")[[1]]
group1 <- args[2]
group2 <- args[3]
output_dir_base <- dirname(dirname(counts_files[1]))
output_dir <- file.path(output_dir_base, "dex")
print(output_dir)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load required packages
#install.packages('BiocManager', repos = "https://cran.rstudio.com/")
#BiocManager::install("DESeq2")
library("DESeq2")
#install.packages('data.table', repos = "https://cran.rstudio.com/")
library('data.table')

# Read count datasets
samples_names <- c()
counts_combined <- NA
print(counts_files)
for (file in counts_files){
  counts <- fread(input = file, skip = "ENS", sep = '\t')
  counts <- data.frame(counts)
  rownames(counts) <- counts[,1]
  count <- data.frame(counts[,-1])
  rownames(count) <- counts[,1]
  counts_combined <- cbind(counts_combined, count)
  sample_name <- gsub("_counts.txt", "", basename(file))
  samples_names <- c(samples_names, sample_name)
}
counts_combined <- counts_combined[,-1]
colnames(counts_combined) <- samples_names

# Define conditions
conditions <- rep(c("group1", "group2"), each = length(samples_names)/2)
coldata <- data.frame(sample = samples_names, condition = conditions)
coldata$condition <- as.factor(coldata$condition)


# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts_combined, colData = coldata, design = ~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$condition <- relevel(dds$condition, ref = "group1")
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds)
res <- res[order(res$padj), ]
write.csv(as.data.frame(res), file = paste0(output_dir, "/DE_analysis.csv"))

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = paste0(output_dir, "/Normalized_counts.csv"))

# Generate heatmap with normalized counts
corM <- cor(normalized_counts)
pdf(paste0(output_dir, "/Heatmap.pdf"))
heatmap(corM, main = "Heatmap of Correlation Matrix", xlab = "Samples", ylab = "Samples", cexRow = 0.5, cexCol = 0.5)
dev.off()

# Generate Volcano Plot
pdf(paste0(output_dir, "/Volcano_Plot.pdf"))
plot(main = 'Volcano plot', res$log2FoldChange, -log10(res$pvalue))
abline(h = -log10(0.05 / nrow(res)), col = 'red')
dev.off()
