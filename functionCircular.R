##Circular heatmap
library(ggplot2)
library(dplyr)
library(tibble)

sigfunc <- read.csv("C:/Users/lifyeo/oa-microbiome/pathway-IRN-full.csv")

e <-sigfunc %>%
  group_by(LevelOne) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  mutate(pathway2 = ifelse(count < 3, "Others", LevelOne)) %>%
  select(-count)%>% 
  arrange(pathway2, Pathway) %>%
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)))


# Create the circular bar plot
A<- ggplot(e, aes(x = factor(Pathway), y = estimate, fill = pathway2)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(start = 0) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_fill_brewer(palette = "Spectral", name = "Pathway class")

A
##Top 20 functional pathway dotplot
et<- e %>% 
  arrange(qval_fdr) %>%
  slice_head(n=15) %>%
  select(Pathway, pathway2, estimate, qval_fdr)

funcolor<- c("Amino acid biosynthesis"="#9E0142" ,"Biosynthesis" = "#D53E4F", 
             "Carbohydrate biosynthesis" = "#F46D43", 
             "Carrier biosynthesis" = "#FDAE61", "Cell wall biosynthesis" = "#FEE08B",
             "CO2 Fixation" = "#FFFFBF" , "Fatty acid and lipid biosynthesis" = "#E6F598",
             "Generation of precursor metabolites and energy" = "#ABDDA4", 
             "Nucleoside and nucleotide biosynthesis" = "#66C2A5", 
             "Others" = "#3288BD", "Vitamin biosynthesis" = "#5E4FA2")
# Create the dot plot
B<- ggplot(et, aes(x = estimate, y = Pathway, fill = pathway2)) +
  geom_col() +
  labs(
    x = "Estimate",
    y = "Pathway",
    color = "Class"
  ) +
  theme_classic() +
  scale_fill_manual(values = funcolor) +
  theme(
    text = element_text(color = "black"),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16))

B
library(patchwork)
nested <- (A / B) +
  plot_annotation(tag_levels = 'A')
nested
ggplot2::ggsave(filename = "nested.pdf", 
                plot = nested,
                #dpi = 300,
                width = 25,
                height = 20,
                units = "in" )

