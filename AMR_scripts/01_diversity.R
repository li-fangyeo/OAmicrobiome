# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsignif)
library(vegan)
library(dunn.test)

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# Distribution ----

# Perform Shapiro-Wilk test on log-transformed values for each tribe
for (tribe in unique(kma_out$tribe)) {
  tribe_data <- kma_out$RPKM[kma_out$tribe == tribe]
  cat("Shapiro-Wilk test for RPKM in", tribe, ":\n")
  print(shapiro.test(tribe_data))
  cat("\n")
}

# kruskal willis
kruskal_result <- kruskal.test(RPKM ~ tribe, data = kma_out)
print(kruskal_result)

# Extracting the chi-squared statistic and p-value
kruskal_chisq <- kruskal_result$statistic
kruskal_pvalue <- kruskal_result$p.value

# dunn test
dunn_test_result <- dunn.test(kma_out$RPKM, kma_out$tribe, method = "bonferroni")

# re-arrange tribe according to least urban to most urban
kma_out$tribe <- factor(kma_out$tribe, levels = c("Jahai", "Temiar", "Temuan", "Malay"))

# visualization
rpkm_boxPlot_tribes <- ggplot(kma_out, aes(x = tribe, y = RPKM, fill = tribe)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_classic() +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  theme(axis.text.x = element_text()) +
  labs(title = "Boxplot of RPKM (log scale) by Tribe", y = "RPKM (log scale)", x = "Sample ID") +
  ggsignif::geom_signif(
    comparisons = list(c("Temiar", "Temuan"), 
                       c("Jahai", "Temiar"), 
                       c("Temiar", "Malay")),
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(7, 6, 6)
  )

plot1 <- rpkm_boxPlot_tribes + 
  annotate("text", x = 1, y = 1e09, 
           label = paste("Kruskal-Wallis chi-squared = ", round(kruskal_chisq, 2),
                         "\np-value = ", format(kruskal_pvalue, scientific = TRUE)), 
           hjust = 0, vjust = 1.5, size = 3.2, color = "red")
plot1

# Save as PNG
# ggsave("resfinder_plot1_distributionBoxPlot.png", plot = plot1, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot1_distributionBoxPlot.pdf", plot = plot1, width = 10, height = 6, units = "in")

# Shannon Index ----

## data preparation
# subset the data
df <- kma_out[c("sampleID", "refSequence", "RPKM")]

# Ensure the RPKM column is numeric
df$RPKM <- as.numeric(df$RPKM)
data_pivot <- pivot_wider(df, names_from = sampleID, values_from = RPKM, values_fill = list(RPKM = 0))

# Convert tibble to data frame
data_pivot <- as.data.frame(data_pivot)

# Set row names to the values in the first column (`refSequence`)
rownames(data_pivot) <- data_pivot[, 1]

# Remove the `refSequence` column since it's now the row names
data_pivot <- data_pivot[, -1]
data_matrix <- as.matrix(data_pivot)

# Transpose the matrix to calculate Shannon index per sample
data_matrix_t <- t(data_matrix)

# Calculate Shannon index per sample
shannon_index <- vegan::diversity(data_matrix_t)

# Add sample names back to the results
names(shannon_index) <- colnames(data_matrix)

# Convert to dataframe
shannon_index_df <- data.frame(
  sampleID = names(shannon_index),
  Shannon_Index = shannon_index
)

# Function to determine tribe based on sampleID prefix
get_tribe <- function(sampleID) {
  prefix <- substring(sampleID, 1, regexpr("[0-9]", sampleID) - 1)
  tribe <- switch(prefix,
                  "J"   = "Jahai", "TRk" = "Temiar",
                  "TM"  = "Temuan", "MLY" = "Malay")
  return(tribe)
}

# add tribe column
shannon_index_df$tribe <- sapply(shannon_index_df$sampleID, get_tribe)

## Normality test
# Q-Q Plot
qqnorm(shannon_index_df$Shannon_Index)
qqline(shannon_index_df$Shannon_Index, col = "red")

# shapiro
shapiro.test(shannon_index_df$Shannon_Index)

## statistical test
kruskal_result <- kruskal.test(Shannon_Index ~ tribe, data = shannon_index_df)

# Extracting the chi-squared statistic and p-value
kruskal_chisq <- kruskal_result$statistic
kruskal_pvalue <- kruskal_result$p.value

# pairwise comparison
dunn.test(shannon_index_df$Shannon_Index, shannon_index_df$tribe, kw = TRUE)

# re-arrange tribe according to least urban to most urban
shannon_index_df$tribe <- factor(shannon_index_df$tribe, 
                                 levels = c("Jahai", "Temiar", "Temuan", "Malay"))

# Create the boxplot
boxplot <- ggplot(shannon_index_df, aes(x = tribe, y = Shannon_Index, fill = tribe)) +
  geom_boxplot(outlier.shape = NA) +
  labs(title = "Shannon Index Distribution by Tribe",
       x = "Tribe",
       y = "Shannon Index") +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  geom_jitter(width = 0.2, color = "blue") +
  theme_classic() +
  theme(axis.text.x = element_text()) +
  ggsignif::geom_signif(
    comparisons = list(c("Jahai", "Malay"), c("Jahai", "Temuan")),
    map_signif_level = TRUE,
    textsize = 3,
    y_position = c(3, 3.2)
  )

plot2 <- boxplot + annotate("text", x = 0.8, y = 4, 
                            label = paste("Kruskal-Wallis chi-squared = ", round(kruskal_chisq, 2),
                                          "\np-value = ", format(kruskal_pvalue, scientific = TRUE)), 
                            hjust = 0, vjust = 1.5, size = 4, color = "red")
plot2

# Save as PNG
# ggsave("resfinder_plot2_shannonBoxPlot.png", plot = plot2, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinder_plot2_shannonBoxPlot.pdf", plot = plot2, width = 10, height = 6, units = "in")