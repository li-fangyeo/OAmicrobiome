# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(tidyverse)
library(ggplot2)

# Data Preparation ----
# Create a list of gene mapping data files in the current directory
fileList <- list.files(pattern = "*.mapstat")

# Initialize an empty list to store all processed gene mapping data
mapstat_conc <- lapply(fileList, function(file) {
  # Extract sampleID from the file name
  sampleID <- sub("\\.mapstat$", "", file)
  
  # Read the gene mapping data from file into a data frame
  data <- read.delim(file, header = TRUE, check.names = FALSE, skip = 6)
  
  # Add sampleID column
  data <- cbind(sampleID, data)
})

# Combine all data frames into a single data frame
mapstat_conc <- do.call(rbind, mapstat_conc)

# rename column
colnames(mapstat_conc)[colnames(mapstat_conc) == "# refSequence"] <- "refSequence"

# add gene length
resfinderFG_meta <- read.csv("resfinderFG_db.csv")
mapstat_conc <- left_join(mapstat_conc, resfinderFG_meta, by = "refSequence")

# Split the first column into four new columns
kma_resfinderFG <- mapstat_conc %>%
  separate(refSequence, into = c("feature", "accession", "source", "antibiotic"), sep = "\\|")

# select columns
kma_resfinderFG <- kma_resfinderFG %>% select(sampleID, feature, accession, source, antibiotic, readCount, Length)

# remove "_FG" suffix at sampleID
kma_resfinderFG$sampleID <- gsub("_FG$", "", kma_resfinderFG$sampleID)

# Calculate total reads per sample
sample_totalReads <- kma_resfinderFG %>%
  group_by(sampleID) %>%
  summarise(total_reads = sum(readCount, na.rm = TRUE)) %>%
  ungroup()

# Join the total reads with the original dataframe
kma_resfinderFG <- kma_resfinderFG %>%
  left_join(sample_totalReads, by = "sampleID")

# Calculate RPKM
kma_resfinderFG <- kma_resfinderFG %>%
  mutate(RPKM = (readCount * 1e9) / (Length * total_reads))

# Function to determine tribe based on sampleID prefix
get_tribe <- function(sampleID) {
  prefix <- substring(sampleID, 1, regexpr("[0-9]", sampleID) - 1)
  tribe <- switch(prefix,
                  "J"   = "Jahai",
                  "TRk" = "Temiar",
                  "TM"  = "Temuan",
                  "MLY" = "Malay")
  return(tribe)
}

# Add tribe column to the merged data frame
kma_resfinderFG$tribe <- sapply(kma_resfinderFG$sampleID, get_tribe)
kma_resfinderFG <- kma_resfinderFG %>% relocate("tribe", .after="sampleID")

head(kma_resfinderFG) %>% knitr::kable() 

# export to csv
# write.csv(kma_resfinderFG, file = "kma_resfinderFG.csv", row.names = FALSE)

# Stacked Bar Chart ----

# calculate percentage within each antibiotic family
amr_summary <- kma_resfinderFG %>%
  group_by(antibiotic, tribe) %>%  
  summarise(
    total_RPKM = sum(RPKM, na.rm = TRUE)
  ) %>%
  group_by(antibiotic) %>%
  mutate(
    percentage = (total_RPKM / sum(total_RPKM)) * 100
  ) %>%
  ungroup()

amr_summary$tribe <- factor(amr_summary$tribe, levels = c("Malay", "Temiar", 
                                                          "Temuan", "Jahai"))

# Create a named vector for the labels
antibiotic_labels <- c(
  "AMC" = "AMC (Amoxicillin/Clavulanic acid)", "AMP" = "AMP (Ampicillin)",
  "AMX" = "AMX (Amoxicillin)", "ATM" = "ATM (Aztreonam)",
  "CAR" = "CAR (Carbenicillin)", "CAZ" = "CAZ (Ceftazidime)",
  "CHL" = "CHL (Chloramphenicol)", "CIP" = "CIP (Ciprofloxacin)",
  "CTX" = "CTX (Cefotaxime)", "CYC" = "CYC (Cycloserine)",
  "FEP" = "FEP (Cefepime)", "GEN" = "GEN (Gentamicin)",
  "KAN" = "KAN (Kanamycin)", "MIN" = "MIN (Minocycline)",
  "OXY" = "OXY (Oxytetracycline)", "PEN" = "PEN (Penicillin)",
  "PIP" = "PIP (Piperacillin)", "SIS" = "SIS (Sisomicin)",
  "SMZ" = "SMZ (Sulfamethoxazole)", "SPT" = "SPT (Spectinomycin)",
  "SXT" = "SXT (Trimethoprim/Sulfamethoxazole)", "TET" = "TET (Tetracycline)",
  "TGC" = "TGC (Tigecycline)", "TMP" = "TMP (Trimethoprim)"
)

# re-arrange tribe according to most urban to least urban
amr_summary$tribe <- factor(amr_summary$tribe, levels = c("Malay", "Temuan", "Temiar", "Jahai"))

resfinderFG_stackedBar <- ggplot(amr_summary, aes(x = antibiotic, 
                                                  y = percentage, 
                                                  fill = tribe)) + 
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(l = 60, r = 30, 
                             t = 20, b = 20,
                             unit = "pt")) +
  scale_fill_manual(values = c("Jahai" = "darkgreen", "Temiar" = "skyblue",
                               "Temuan" = "orange", "Malay" = "pink")) +
  scale_x_discrete(labels = antibiotic_labels) +
  labs(x = "Antibiotic Family", y = "Proportion in Total Pool",
       fill = "Tribe", title = "Distribution of functional genes in each tribes (in %)")

resfinderFG_stackedBar

# Save as PNG
# ggsave("resfinderFG_plot1_stackedBar.png", plot = resfinderFG_stackedBar, width = 10, height = 6, units = "in")

# Save as PDF
# ggsave("resfinderFG_plot1_stackedBar.pdf", plot = resfinderFG_stackedBar, width = 10, height = 8, units = "in")