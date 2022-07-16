
# Create glial protrusion library table -------------------------------------------------------

library(tidyverse)
library(readxl)
library(qs)


# Get plotting dataframe ----------------------------------------------------------------------

## Glial data
avgTPM_wide <- readRDS("./RNAseq/quant_results/all-glia_avgTPM_wide.RDS")
protrusion_RNA_library_info <- read_csv("./RNAseq/protrusion_RNA_library_info.csv") %>%
  unite(col = "index", c(Model.system, Species, Separation.method, Data.type, Reference), sep = "\n\n")
protrusion_RNA_libraries <- protrusion_RNA_library_info$library

commonly_expressed_genes <- avgTPM %>%
  group_by(gene_id) %>%
  summarise(n_exp_lib = sum(TPM > 10)) %>% ungroup() %>%
  filter(n_exp_lib >= 3) %>%
  pull(gene_id)

## Neurite data
neurite_enriched <- read_excel("./RNAseq/external_data/Supplementary online tables_neurite-enriched.xlsx",
                               skip = 1) %>% janitor::clean_names()

neurite_present <- read_excel("./RNAseq/external_data/Supplementary online tables_Chekulaeva_most-abundant.xlsx",
                              skip = 1) %>% janitor::clean_names()

## Summarise
avgTPM <- avgTPM_wide %>%
  dplyr::select(gene_id, protrusion_RNA_libraries) %>%
  pivot_longer(cols = contains(c("trap", "txn")),
               names_to = "library",
               values_to = "TPM") %>%
  mutate(protrusion_RNA = if_else(TPM > 10, TRUE, FALSE)) %>%
  mutate(is_common = if_else(gene_id %in% commonly_expressed_genes, "Detected in ≥ 3 datasets", "Detected in < 3 datasets")) %>%
  mutate(is_neurite_present = if_else(gene_id %in% neurite_present$gene_id, "Shared with neurites", "Glia only"))

# avgTPM %>% qsave("./RNAseq/misc/avgTPM_for_plotting.qs")

avgTPM_summary_glia_dataset <- avgTPM %>%
  filter(protrusion_RNA == TRUE) %>%
  group_by(library, is_common) %>%
  summarise(count = n()) %>% ungroup() %>%
  left_join(protrusion_RNA_library_info)

avgTPM_summary_glia_neurite <- avgTPM %>%
  filter(protrusion_RNA == TRUE) %>%
  group_by(library, is_neurite_present) %>%
  summarise(count = n()) %>% ungroup() %>%
  left_join(protrusion_RNA_library_info)

expression("">=3)
# Plot overall study summary ------------------------------------------------------------------

avgTPM_summary_glia_dataset %>%
  mutate(index = str_replace(index, "trans", "\ntrans")) %>%
  mutate(index = str_replace(index, "al,", "al,\n")) %>%
  # mutate(is_common, str_replace(is_common, "≥", "\u2265")) %>%
  ggplot(aes(
    x = index,
    y = count, 
    fill = is_common
  )) +
  labs(title = "Glial protrusion transcriptome meta-analysis overview",
       subtitle = "", 
       fill = "",
       x = "",
       y = "Detected gene count") + 
  geom_col(width = 0.6, alpha = 0.9) + 
  scale_fill_manual(values = c("gray75", "lightsalmon3")) +
  scale_y_continuous(expand = c(0, 100)) + 
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(size = 17),
        axis.text.y = element_text(size = 10))

ggsave("./RNAseq/plots/meta-analysis-overview.pdf", 
       width = 13.5, height = 6.5,
       device = cairo_pdf)
 























