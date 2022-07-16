
# Perform gene ontology analysis and make plots -----------------------------------------------

library(tidyverse)
source("./src/runTopGO.R")


# Run TopGO -----------------------------------------------------------------------------------

## Background brain transcriptome
Mmus_brain_glia_txome <- readRDS("./RNAseq/external_data/Mmus_brain_glia_txome.RDS")

## Summary dataframe (mouse)
df <- readRDS("./summary_table.RDS")

## Ontologies to test
ontologies <- c("BP", "MF", "CC") %>% set_names()

## Run analyses
topgo_rna_present <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, RNA_in_protrusion >= 8), gene_id))
}, .id = "ontology")

topgo_translation_present <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, translation_in_protrusion >= 4), gene_id))
}, .id = "ontology")

topgo_rna_and_translation_present <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, RNA_in_protrusion >= 8 & translation_in_protrusion >= 4), gene_id))
}, .id = "ontology")

topgo_rna_enriched <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, enriched_in_protrusion >= 3), gene_id))
}, .id = "ontology")

topgo_translation_enriched <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, enhanced_translation_in_protrusion >= 3), gene_id))
}, .id = "ontology")

topgo_rna_and_translation_enriched <- map_dfr(ontologies, function(ontology){
  runTopGO(ontology = ontology,
           species = "mouse",
           backgroundGenes = Mmus_brain_glia_txome,
           targetedGenes = pull(filter(df, enriched_in_protrusion >= 3 & enhanced_translation_in_protrusion >= 3), gene_id))
}, .id = "ontology")


# Plots ---------------------------------------------------------------------------------------

## Plot helpers
ontology_order <- c("BP", "MF", "CC")
ontology_colours <- c("#d1ab75", "#8fd175", "#d175b8")
plot_topgo <- function(df, bonferroni_th, foldEnrichment_th, cc_sig, other_sig, title){
  df %>%
    filter(bonferroni < bonferroni_th) %>%
    filter(foldEnrichment > foldEnrichment_th) %>%
    filter(case_when(
      ontology == "CC" ~ Significant > cc_sig,
      TRUE ~ Significant > other_sig
    )) %>%
    mutate(plot_FDR = -log10(bonferroni)) %>%
    mutate(ontology = fct_relevel(ontology, ontology_order)) %>%
    ggplot(aes(x = reorder(Term, plot_FDR), y = plot_FDR)) +
    geom_segment(aes(xend = reorder(Term, plot_FDR), y = 0, yend = plot_FDR), colour = "gray70") +
    geom_point(aes(colour = ontology), size = 4) +
    geom_label(aes(y = plot_FDR + 2, label = Significant), cex = 2) + 
    scale_colour_manual(values = ontology_colours,
                        labels = c("Biological process", "Molecular function", "Cellular component")) +
    labs(title = title,
         subtitle = "",
         y = "-log10 Bonferroni p-value",
         colour = "") + 
    coord_flip() +
    theme_classic(base_size = 9) + 
    theme(axis.title.y = element_blank(),
          legend.position = "bottom",
          axis.text.y = element_text(size = 10))
}

## Plot
plot_topgo(topgo_rna_present, 0.05, 2, 189, 100, "GO: RNA present in glial protrusions")
ggsave("./RNAseq/plots/GO_rna_present.pdf", width = 7, height = 5)

plot_topgo(topgo_translation_present, 0.05, 2, 200, 100, "GO: Translation present in glial protrusions")
ggsave("./RNAseq/plots/GO_translation_present.pdf", width = 7, height = 5)

plot_topgo(topgo_rna_and_translation_present, 0.05, 2, 200, 100, "GO: RNA & Translation present in glial protrusions")
ggsave("./RNAseq/plots/GO_rna_and_translation_present.pdf", width = 7, height = 5)

plot_topgo(topgo_rna_enriched, 0.05, 2, 40, 30, "GO: RNA enriched in glial protrusions")
ggsave("./RNAseq/plots/GO_rna_enriched.pdf", width = 7, height = 5)

plot_topgo(topgo_translation_enriched, 0.05, 1.5, 20, 20, "GO: Translation enriched in glial protrusions")
ggsave("./RNAseq/plots/GO_translation_enriched.pdf", width = 7, height = 5)

plot_topgo(topgo_rna_and_translation_enriched, 0.05, 1, 5, 5, "GO: RNA & Translation enriched in glial protrusions")
ggsave("./RNAseq/plots/GO_rna_and_translation_enriched.pdf", width = 7, height = 5)































