# Set current directory to source file location
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path) 
setwd(current_working_dir)

# Load necessary libraries
library(readxl)
library(tidyverse) 
library(ggvenn)

# Function ----

# Function to get all pairwise intersections
get_intersections <- function(gene_list) {
  all_combinations <- unlist(lapply(seq_along(gene_list), function(n) combn(names(gene_list), n, simplify = FALSE)), recursive = FALSE)
  intersection_list <- lapply(all_combinations, function(tribe_names) {
    intersected_genes <- Reduce(intersect, gene_list[tribe_names])
    data.frame(Tribes = paste(tribe_names, collapse = " & "), Genes = paste(intersected_genes, collapse = ", "))
  })
  do.call(rbind, intersection_list)
}

# Function to get exact intersections between two sets
exact_intersection <- function(set1, set2, others) {
  # Intersection of the two sets
  intersect_sets <- intersect(set1, set2)
  
  # Subtract genes that are in the intersection with any of the other sets
  for (other_set in others) {
    intersect_sets <- setdiff(intersect_sets, intersect(intersect_sets, other_set))
  }
  
  return(intersect_sets)
}

# Function to get intersections among three sets
exact_intersection_three <- function(set1, set2, set3, other) {
  # Intersection of the three sets
  intersect_sets <- intersect(intersect(set1, set2), set3)
  
  # Subtract genes that are in the intersection with any of the other set
  intersect_sets <- setdiff(intersect_sets, intersect(intersect(intersect(set1, set2), set3), other))
  
  return(intersect_sets)
}

# Function to get intersections among four sets
exact_intersection_four <- function(set1, set2, set3, set4) {
  # Intersection of all four sets
  intersect_sets <- intersect(intersect(intersect(set1, set2), set3), set4)
  
  return(intersect_sets)
}

# Data Preparation ----

# import data
kma_out <- read.csv("kma_out.csv", check.names = FALSE)

# extract gene of interest
refSequence_meta <- read_excel("unique_refSeq.xlsx", sheet = "Sheet2")

gene_of_interest <- refSequence_meta %>% 
  filter(includeHeatmap == "yes") %>% 
  pull(refSequence)

kma_geneOfInterest <- kma_out %>% 
  filter(refSequence %in% gene_of_interest)

# Venn ----

# Split data by tribe and extract unique Gene terms
tribe_list <- kma_out %>%
  group_by(tribe) %>%
  reframe(refSequence = list(unique(refSequence))) %>%
  # Turn the reframe result into a named list
  reframe(tribe_list = setNames(refSequence, unique(tribe))) %>%
  pull(tribe_list)

venn <- ggvenn(tribe_list, fill_color = c("darkgreen", "pink", "skyblue", "orange")) +
  ggtitle("Venn Diagram of Unique AMR Gene Among Tribes")
venn

# Save as PNG
# ggsave("resfinder_plot4_vennDiagram.png", plot = venn, width = 10, height = 8, units = "in", bg = "white")

# save as pdf
# ggsave("resfinder_plot4_vennDiagram.pdf", plot = venn, width = 10, height = 8, units = "in", bg = "white")

# Unshared Gene within Each Tribe ----

unsharedGene_list <- lapply(names(tribe_list), function(tribe) {
  # Find genes in the current tribe
  current_tribe_genes <- tribe_list[[tribe]]
  
  # Combine genes from all other tribes
  other_tribes_genes <- unlist(tribe_list[names(tribe_list) != tribe])
  
  # Find genes unique to the current tribe
  unique_to_tribe <- setdiff(current_tribe_genes, other_tribes_genes)
  
  return(unique_to_tribe)
})

names(unsharedGene_list) <- names(tribe_list)

# Convert list to long format using stack()
unsharedGene <- stack(unsharedGene_list)

# Rename the columns
colnames(unsharedGene) <- c("refSequence", "tribe")
unsharedGene <- unsharedGene[, c("tribe", "refSequence")]

# Collapse the data frame with count of genes
unsharedGene <- unsharedGene %>%
  group_by(tribe) %>%
  summarise(
    Gene_List = paste(sort(refSequence), collapse = ", "),
    Count = n(),
    .groups = 'drop'
  )

head(unsharedGene) %>% knitr::kable() 

# write.csv(unsharedGene, file = "unsharedgenes_allGene.csv")

# Shared Genes (Inclusive Sets) ----

# Get the intersections
sharedGenes <- get_intersections(tribe_list)

# Sort the genes within the Genes column and update the column
sharedGenes$Genes <- sapply(strsplit(as.character(sharedGenes$Genes), ","), function(genes) {
  sorted_genes <- sort(genes)  
  return(paste(sorted_genes, collapse = ", ")) 
})

# Now calculate the gene count for the sorted Genes column
sharedGenes$gene_count <- sapply(strsplit(as.character(sharedGenes$Genes), ","), length)

head(sharedGenes) %>% knitr::kable() 

# write.csv(sharedGenes, file = "sharedGenes_allGene.csv")

# Shared Genes (Exclusive Sets) ----

# Obtain exact intersections
only_Jahai_Malay <- exact_intersection(
  tribe_list$Jahai,
  tribe_list$Malay,
  list(tribe_list$Temuan, tribe_list$Temiar)
)

only_Jahai_Temuan <- exact_intersection(
  tribe_list$Jahai,
  tribe_list$Temuan,
  list(tribe_list$Malay, tribe_list$Temiar)
)

only_Jahai_Temiar <- exact_intersection(
  tribe_list$Jahai,
  tribe_list$Temiar,
  list(tribe_list$Malay, tribe_list$Temuan)
)

only_Malay_Temuan <- exact_intersection(
  tribe_list$Malay,
  tribe_list$Temuan,
  list(tribe_list$Jahai, tribe_list$Temiar)
)

only_Malay_Temiar <- exact_intersection(
  tribe_list$Malay,
  tribe_list$Temiar,
  list(tribe_list$Jahai, tribe_list$Temuan)
)

only_Temuan_Temiar <- exact_intersection(
  tribe_list$Temuan,
  tribe_list$Temiar,
  list(tribe_list$Jahai, tribe_list$Malay)
)

only_Jahai_Malay_Temuan <- exact_intersection_three(
  tribe_list$Jahai,
  tribe_list$Malay,
  tribe_list$Temuan,
  tribe_list$Temiar
)

only_Jahai_Malay_Temiar <- exact_intersection_three(
  tribe_list$Jahai,
  tribe_list$Malay,
  tribe_list$Temiar,
  tribe_list$Temuan
)

only_Jahai_Temuan_Temiar <- exact_intersection_three(
  tribe_list$Jahai,
  tribe_list$Temuan,
  tribe_list$Temiar,
  tribe_list$Malay
)

only_Malay_Temuan_Temiar <- exact_intersection_three(
  tribe_list$Malay,
  tribe_list$Temuan,
  tribe_list$Temiar,
  tribe_list$Jahai
)

only_Jahai_Malay_Temuan_Temiar <- exact_intersection_four(
  tribe_list$Jahai,
  tribe_list$Malay,
  tribe_list$Temuan,
  tribe_list$Temiar
)

# Output results
shared_refGene_exact_list <- list(
  Jahai_Malay = only_Jahai_Malay,
  Jahai_Temuan = only_Jahai_Temuan,
  Jahai_Temiar = only_Jahai_Temiar,
  Malay_Temuan = only_Malay_Temuan,
  Malay_Temiar = only_Malay_Temiar,
  Temuan_Temiar = only_Temuan_Temiar,
  Jahai_Malay_Temuan = only_Jahai_Malay_Temuan,
  Jahai_Malay_Temiar = only_Jahai_Malay_Temiar,
  Jahai_Temuan_Temiar = only_Jahai_Temuan_Temiar,
  Malay_Temuan_Temiar = only_Malay_Temuan_Temiar,
  Jahai_Malay_Temuan_Temiar = only_Jahai_Malay_Temuan_Temiar
)

# Convert list to data frame
shared_refGene_exact <- stack(shared_refGene_exact_list)
names(shared_refGene_exact) <- c("Gene", "Shared_Between")

# Collapse the data frame with count of genes
shared_refGene_exact <- shared_refGene_exact %>%
  group_by(Shared_Between) %>%
  summarise(
    Gene_List = paste(sort(Gene), collapse = ", "),
    Count = n(),
    .groups = 'drop'
  )

# write.csv(shared_refGene_exact, file = "sharedGene_allGene_exact.csv")