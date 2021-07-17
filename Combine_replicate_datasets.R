#!/usr/bin/env Rscript


### Script to combine replicates from CLIP data. Takes the outputs for each replicate
# from the iCLIP_Genialis_feature_annotation.R as inputs, together with the number
# of replicates in which a crosslink site must have been identified to be retained

# Can either provide only the significant crosslink sites (without needing an FDR column),
# or all crosslink sites in which case they will be filtered for the FDR threshold provided

# To run:
# Ensure tidyverse package is installed (available from CRAN)
# Amend dataset list, output file name, number of replicates and FDR threshold for the data to be analysed
# save and source script

# Provide list of annotated replicate datasets:
datasets <- list("iCLIP_rep1_chr17_100kb_annotated.tsv", "iCLIP_rep2_chr17_100kb_annotated.tsv", "iCLIP_rep3_chr17_100kb_annotated.tsv")

# Provide output filename for the combined data
output_file <- "iCLIP_combined_replicates.txt"

# Number of replicates in which crosslink site must be detected:
min_Replicates <- 2

# FDR threshold to filter data (if 'FDR' column present)
FDRthresh <- 0.05


library(tidyverse)

lapply(datasets, read_tsv, col_types = cols(chromosome = col_character())) %>%
  bind_rows(.id = "ID") %>%
  distinct(ID, chromosome, start, strand, .keep_all = T) -> all_data
# for all crosslink sites from iMaps, some sites can be duplicated with different gene
# names, so distinct() is necessary to removed these. Gene names/IDs/Features that were
# added by the iCLIP_Genialis_feature_annotation script will be identical for the duplicate
# sites.

if("FDR" %in% colnames(all_data)) {
  all_data %>%
    filter(FDR <= FDRthresh) -> all_data
}

all_data %>%
  group_by(chromosome, start, strand, Gene_name, Gene_ID, Feature) %>%
  summarise(nReplicates = n(), mean_score = mean(score), sum_score = sum(score)) %>%
  ungroup() %>%
  rename(position = start) %>%
  filter(nReplicates >= min_Replicates) -> filtered_data

filtered_data %>% write_tsv(output_file)
