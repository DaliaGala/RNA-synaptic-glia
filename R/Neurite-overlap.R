
# Overlap with Neurite data -----------------------------------------------

library(tidyverse)
library(readxl)


# Prepare data ------------------------------------------------------------

## Neurite data
neurite_enriched <- read_excel("./RNAseq/external_data/Supplementary online tables_neurite-enriched.xlsx",
                               skip = 1) %>% janitor::clean_names()

neurite_present <- read_excel("./RNAseq/external_data/Supplementary online tables_Chekulaeva_most-abundant.xlsx",
                              skip = 1) %>% janitor::clean_names()

## Glial data 
df <- readRDS("./summary_table.RDS")

## Summarise 
neurite_rna_present_n8 <- neurite_present %>% 
  filter(datasets_with_neurite_tpm_10 >= 8) %>% pull(gene_id)
neurite_translation_present_n2 <- neurite_present %>%
  filter(studies_with_neurite_ribosome_association >= 2) %>% pull(gene_id)
neurite_rna_enriched_n6 <- neurite_enriched %>%
  filter(datasets_with_significant_neurite_enrichment_p_0_1 >= 6) %>% pull(gene_id)

df_with_neurite <- df %>%
  mutate(neurite_rna_present = if_else(gene_id %in% neurite_rna_present_n8, TRUE, FALSE)) %>%
  mutate(neurite_translation_present = if_else(gene_id %in% neurite_translation_present_n2, TRUE, FALSE)) %>%
  mutate(neurite_rna_enriched = if_else(gene_id %in% neurite_rna_enriched_n6, TRUE, FALSE))

df_with_neurite_summary <- df_with_neurite %>%
  mutate(RNA_present = case_when(
    RNA_in_protrusion >= 8 & neurite_rna_present == TRUE ~ "Both",
    RNA_in_protrusion >= 8 ~ "Glia only",
    neurite_rna_present == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(translation_present = case_when(
    translation_in_protrusion >= 4 & neurite_translation_present == TRUE ~ "Both",
    translation_in_protrusion >= 4 ~ "Glia only",
    neurite_translation_present == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  mutate(RNA_enriched = case_when(
    enriched_in_protrusion >= 3 & neurite_rna_enriched == TRUE ~ "Both",
    enriched_in_protrusion >= 3 ~ "Glia only",
    neurite_rna_enriched == TRUE ~ "Neurite only",
    TRUE ~ as.character(NA)
  )) %>%
  dplyr::select(RNA_present, translation_present, RNA_enriched) %>%
  pivot_longer(cols = everything(),
               names_to = "type",
               values_to = "celltype") %>%
  group_by(type, celltype) %>%
  summarise(count = n()) %>% ungroup() %>%
  filter(!is.na(celltype)) %>%
  group_by(type) %>%
  mutate(total_count = sum(count)) %>% ungroup() %>%
  mutate(percentage = count/total_count * 100) %>%
  mutate(label = paste0(round(percentage, digits = 1), "%\n", "(", count, ")")) %>%
  mutate(celltype = str_replace_all(celltype, " ", "\n")) %>%
  mutate(type = str_replace_all(type, "_", "\n")) %>%
  mutate(type = str_replace_all(type, "trans", "Trans")) %>%
  mutate(type = fct_relevel(type, c("RNA\npresent", "RNA\nenriched", "Translation\npresent"))) 


# Plot --------------------------------------------------------------------
fill_colour <- c("#d18975", "#8fd175", "#c9d175")

df_with_neurite_summary %>%
  ggplot(aes(
    x = type,
    y = count,
    fill = celltype
  )) +
  geom_col(width = 0.5, alpha = 1) +
  geom_text(aes(
    x = 1.38,
    label = celltype
  ), position = position_stack(vjust = 0.5), cex = 2) +
  geom_text(aes(
    x = 1.25,
    label = "-"
  ), position = position_stack(vjust = 0.5), cex = 2) +
  geom_label(aes(
    x = type,
    # y = count,
    label = label,
    group = celltype
  ), position = position_stack(vjust = 0.5), cex = 2.25, fill = "white", alpha = 0.5, fontface = "bold") + 
  labs(title = "Overlap between Glial protrusion & Neurite transcriptomes",
       subtitle = "",
       x = "",
       y = "Gene count",
       fill = "") + 
  scale_fill_manual(values = fill_colour) + 
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ type, scales = "free_x") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10))

ggsave("./RNAseq/plots/overlap_with_neurite_data_bar_free-x.pdf",
       width = 4.6, height = 5.2)
  























