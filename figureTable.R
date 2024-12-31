pheatmap(mat, annotation_row = taxa_clusters, 
         annotation_col = sample_data,
         annotation_colors = ann_colors,
         breaks = breaks,
         color = colors,
         fontsize = 12,
         show_colnames = FALSE)
#Jahai as ref
#2nd now Malay as ref group
malay <- df %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  select(taxon, contains("Malay")) %>%
  filter(if_any(contains("q"), ~ . < 0.05)) 
temuan
malay
temiar
jahai

#rename the columns
my <- function(x, y) {
  x %>%
    rename("Log-fold change" = !!sym(paste0("lfc_Group", y)), 
           "SE" = !!sym(paste0("se_Group", y)), 
           "W" = !!sym(paste0("W_Group", y)),
           "FDR p-value" = !!sym(paste0("q_Group", y)),
           "Sensitivity test" = !!sym(paste0("passed_ss_Group", y)))
}

malay %>%
  #drop these columns
  select(-c(diff_GroupMalay, p_GroupMalay))%>%
  #arrange by alphabetic reversed
  arrange(desc(passed_ss_GroupMalay)) %>%
  #rename columns
  my("Malay") %>%
  #remove _ in the taxon
  mutate(taxon = str_replace_all(taxon, "_", " ")) %>%
  write.csv("./ancombc/temiar-ref/Malay-ancombc.csv")


