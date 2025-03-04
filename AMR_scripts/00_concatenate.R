# remove previous session
rm(list = ls()); closeAllConnections()

# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(dplyr)

# Create a list of gene mapping data files in the current directory
fileList <- list.files(pattern = "*.mapstat")

# Initialize an empty list to store all processed gene mapping data
mapstat_conc <- lapply(fileList, function(file) {
  # Extract sampleID from the file name
  sampleID <- sub("\\.mapstat$", "", file)

  # Read the gene mapping data from file into a data frame
  data <- read.delim(file, header = TRUE, check.names = FALSE, skip = 6)

  # Add sampleID column to the left of the data frame without modifying column names
  data <- cbind(sampleID, data)
})

# Combine all data frames into a single data frame
mapstat_conc <- do.call(rbind, mapstat_conc)

# rename column
colnames(mapstat_conc)[colnames(mapstat_conc) == "# refSequence"] <- "refSequence"

# add gene length
resfinder_meta <- read.csv("resfinder_db.csv")
mapstat_conc <- left_join(mapstat_conc, resfinder_meta, by = "refSequence")

# select columns
kma_out <- mapstat_conc %>% select(sampleID, refSequence, geneLength, readCount)

# Calculate total reads per sample
sample_totalReads <- kma_out %>%
  group_by(sampleID) %>%
  summarise(total_reads = sum(readCount, na.rm = TRUE)) %>%
  ungroup()

# Join the total reads with the original dataframe
kma_out <- kma_out %>%
  left_join(sample_totalReads, by = "sampleID")

# Calculate RPKM
kma_out <- kma_out %>%
  mutate(RPKM = (readCount * 1e9) / (geneLength * total_reads))

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
kma_out$tribe <- sapply(kma_out$sampleID, get_tribe)
kma_out <- kma_out %>% relocate("tribe", .after="sampleID")

# export to csv
# write.csv(kma_out, file = "kma_out.csv", row.names = FALSE)
