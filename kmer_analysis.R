#!/usr/bin/env Rscript

# Best run on cluster (it will probably take > 1day run locally..!) and one file per job.

# GenomicRanges, GenomicFeatures, rtracklayer and BSgenome.Mmusculus.UCSC.mm10 are available through Bioconductor (BiocManager::install())
# tidyverse and parallel packages available through CRAN (install.packages())

# Example command to submit to Babraham cluster:
# ssub -c 4 --mem=50G --email -o kmer_analysis.out Rscript --vanilla kmer_analysis.R output_dir=kmer_results/ ~/Annotations/Mus_musculus.GRCm38.96.gtf.gz ptbp1-cd8-24h-rep3_mapped_to_genome_single_peaks_annotated.tsv.gz k=8 p=4 window_size=31 sample_features=T sample_dataset=10000

# Except for GTF file and file to analyse which must be provided, all other arguments have defaults (see below)

# Based on method from Wang, Ule et al, PLoS Biol, 2010, iCLIP predicts the dual splicing effects of TIA-RNA interactions, PMID: 21048981:
# (i) iCLIP reads were associated with expressed genomic regions as defined by ENSEMBL (version Hg18/NCBI36). Each coding or non-coding gene was defined as its own region (in case of overlapping genes, the shorter gene always had the priority). Introns, 5' UTR, ORF and 3' UTR were considered as separate regions. (ii) iCLIP reads antisense to the transcriptional direction of the associated gene and reads that mapped to non-annotated genomic regions were removed before proceeding to further analysis.  (iii) The control files were generated 100 times with randomised iCLIP positions. iv) Both in iCLIP and control files, the positions were extended by 10 nt in both directions, such that 21 nt long sequences were used for analysis. v) The occurrence (pentamer frequency) was calculated for each pentamer in each file. vi) The z-score was calculated for each pentamer as:
#(occurrence in iCLIP sequences - average occurrence in control sequences) / standard deviation of occurrence in control sequences

# One issue I came across for larger values of k was that very occasionally, no occurrences were found in the randomised sites.
# This was when running tests with only a few sample sets, but is theoretically possible
# To avoid this being an issue and giving infinite z-scores, we add 1 to the standard deviation before dividing

# Takes annotated crosslink site files, which should include gene/feature annotations, and assesses enrichment of kmers relative to random positions
# Positions are randomised across all genes with at least one xlink site, with the same distribution across features as in the iCLIP dataset
# ie, the random sites will include the same number of sites within the 3'UTR as there are 3'UTR xlink sites, etc. 

# Input files must include a GTF file (ideally the same as was used for iCLIP data annotation)

# One or more iCLIP xlink files can be analysed - but note that it takes quite a long time to run so may be more advisable to run one at once (esp for large k). 
# They should be tab-delimited text files containing significant crosslink sites - ie width 1, FDR < 0.05, not non-significant or clusters
# Column headers must include chromosome, either start or position, strand, Gene_name, Feature
# BED files directly from iMaps are not suitable as they are not annotated with genes and features - first run iCLIP_Genialis_feature_annotation.R 
# On the command line, argument with extension .gtf(.gz) will be recognised as the GTF file
# Anything else that doesn't correspond to another argument will be treated as a file to analyse

# Users can define the size of the kmers as [k] (default 5) by including the argument k=[k]
# Users can define the size of the window assessed [w] (default 21, ie 10 bases up/downstream of site)
# by including the argument window_size=[w]

# By default will analyse all crosslink sites in the file.
# If instead you want to analyse a random subset of size [s], include the argument sample_dataset=[s]
# This would be useful eg for more consistent comparison between datasets. The individual feature zscores would still not be identical in size, but for the 
# same type of data would likely be similar. If [s] cannot be coerced to numeric, it will be ignored.
# The sampling is done after removing non-coding/intergenic sites and unlocalised chromosome scaffolds (so need to take this into account when choosing the size), 
# therefore for the whole gene, the number of sites considered will be identical if the same value of [s] is used across different datasets, as long as there are enough. 
# Note that if there are not sufficient sites in the dataset, all sites will be retained; this will be flagged in the R output file.
# Including the argument sample_features=T will result in an additional column per feature, where all are downsampled to the smallest number of crosslink sites.
# NB if used together, the dataset sample will be taken first, so feature samples could end up quite small!
# Alternatively, sample_features=[f] will attempt to downsample all features to size [f], or take all crosslink sites if there are fewer than [f]
# If you use this, check the R output file as this will report if any features did not have sufficient crosslink sites.

# Intergenic xlink sites and those in ncRNAs will not be included in the analysis. 
# Other feature types will be analysed separately (CDS, 3'UTR, 5'UTR, Intron). 
# Additionally, zscore for all together will also be calculated.


# For each input file, an output file corresponding to the input file name will be saved in the output directory - 
# defined using the argument output_dir=path/to/directory (default is current working directory - which will be the 
# directory from which you submit the job if on the cluster).
# File extension of input filename will be replaced in output with eg "_5mer_zscores.txt" (where 5 is replaced by value of k)
# Contains kmers in column 1, then additional columns with z scores for all feature types together, or each separately

# To use parallel processing (will speed it up), set the number of threads to parallelise over [p] as p=[p]
# This is not possible on Windows where p must equal 1 (default). 
# On the cluster, p can equal 2x number of cores specified to ssub (as each can handle 2 threads)

# kmers containing 'N' are removed from the output (sometimes found in introns)

# Note that occasionally the cluster job appears to 'fail' with an error message at the end of the log file - but having 
# already finished processing and writing out the output.
# There is therefore a line printed when processing of each file is started and finished - if 'Finished processing [filename]'
# appears in the log file, the output for that file was successfully written out before it 'failed'.

library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)
library(tidyverse)

#########################################################################################
# Code chunk 1a: Parse arguments - comment out 1b if reading arguments from command line
#########################################################################################

Args = commandArgs(trailingOnly=TRUE)
Args <- c("iCLIP_combined_replicates.txt", "Mouse_GRCm38.90_chr17_100kb.gtf")
gtf_file <- Args[grepl("gtf$|gtf.gz$", Args)]
if(length(gtf_file) > 1) {
  warning("More than one GTF file provided; using only first listed")
  gtf_file <- gtf_file[1]
}
if(length(gtf_file) == 0) {
  stop("No GTF file provided")
} 

files_to_analyse <- Args[!(grepl("^k=|^output_dir=|^window_size=|^p=|^sample_features=|^sample_dataset=", Args) | grepl("gtf$|gtf.gz$", Args))]
if(length(files_to_analyse) == 0) {
  stop("No files provided for analysis")
}

if(sum(grepl("^k=", Args)) > 0) {
  k <- as.numeric(sub("^k=", "", Args[grepl("k=", Args)]))
} else { k <- 5 }

if(sum(grepl("^p=", Args)) > 0) {
  p <- as.numeric(sub("^p=", "", Args[grepl("p=", Args)]))
} else { p <- 1 }

if(sum(grepl("^window_size=", Args)) > 0) {
  window_size <- as.numeric(sub("window_size=", "", Args[grepl("^window_size=", Args)]))
} else { window_size <- 21 }

if(sum(grepl("^output_dir=", Args)) > 0) {
  output_dir <- sub("output_dir=", "", Args[grepl("^output_dir=", Args)])
} else { output_dir <- "." }

if(sum(grepl("^sample_features=", Args)) > 0) {
  sample_features <- sub("^sample_features=", "", Args[grepl("sample_features=", Args)])
  if(sample_features == "T") {
    sample_features <- T
  } else {
    sample_features <- as.numeric(sample_features)
  }
} else { sample_features <- NA }

if(sum(grepl("^sample_dataset=", Args)) > 0) {
  sample_dataset <- as.numeric(sub("sample_dataset=", "", Args[grepl("^sample_dataset=", Args)]))
} else { sample_dataset <- NA }


########################################################################################
# Code chunk 1b: Provide arguments - comment out 1a if providing the arguments yourself
########################################################################################

# Uncomment lines below and modify as appropriate:


# gtf_file <- "path/to/gtffile.gtf.gz"

# files_to_analyse <- c("path/to/textfile") # more than one can be provided, but will take a long time to run!!

# k <- 5 # size of kmer to search for

# window_size <- 21 # size of extended window around each site

# output_dir <- "." # change if you want output saved elsewhere

# p <- 1 # Number of threads to parallelise over - on Windows must equal 1, parallel processing not possible

# sample_dataset <- NA

# sample_features <- NA


###################################################################################################################
# Code chunk 2: Define all possible kmers of length k
###################################################################################################################

tibble(var = c("A", "C","G","T")) -> all_kmers
for(position in 2:k) { all_kmers %>% add_column(var = c("A","C","G","T"), .name_repair = "unique") -> all_kmers}
all_kmers %>% 
  expand.grid() %>% 
  unite("kmer", everything(), sep = "") %>%
  deframe() -> 
  all_kmers

###################################################################################################################
# Code chunk 3: Import GTF file, extract features into separate GRanges objects for each type (excluding low TSL)
###################################################################################################################

gtf <- import(gtf_file)
seqlevelsStyle(gtf) <- "UCSC"
gtf <- gtf[seqnames(gtf) %in% seqnames(BSgenome.Mmusculus.UCSC.mm10)]
# compatible with either/both/neither having chr in names.
# won't work for MT but I think for this it's fine to exclude

# Some GTF file eg from Gencode have 'gene_type' and 'transcript_type' instead of 'gene/transcript_biotype' 
# Also have 'ccdsid' rather than 'ccds_id'
# So we need to rename these if that's the case to ensure the rest of the code works
names(mcols(gtf))[names(mcols(gtf)) == "transcript_type"] <- "transcript_biotype"
names(mcols(gtf))[names(mcols(gtf)) == "gene_type"] <- "gene_biotype"
names(mcols(gtf))[names(mcols(gtf)) == "ccdsid"] <- "ccds_id"

# First, eliminate all transcript isoforms that have a TSL > 3 (or NA)
# transcript_support_level in the gtf file (at least v101 which I'm using to test..!) contains many annotations where the TSL is followed by a string in brackets, eg:
# "1 (assigned to previous version 13)". So first need to strip these out, before converting to numeric and selecting only those <=3
gtf$transcript_support_level <- as.numeric(sub(" .*$", "", gtf$transcript_support_level))
gtf <- gtf[!is.na(gtf$transcript_support_level)]
gtf <- gtf[gtf$transcript_support_level <=3]


# Include immunoglobulin & Tcr genes in protein coding ones. Exclude anything with pseudogene in transcript_biotype, but include any isoforms with protein_coding as transcript_biotype (can be designated pseudogenes if you just look at gene_biotype)
pc_genes <- gtf[(gtf$gene_biotype=="protein_coding"| grepl("^[IT][GR]_.*_gene$",gtf$gene_biotype) | gtf$transcript_biotype %in% "protein_coding") & !(grepl("pseudogene", gtf$transcript_biotype))]

## so this is taking everything that has a GENE biotype of protein_coding, and all transcripts from that gene regardless of whether they are coding or not (unless transcript biotype includes pseudogene). 
# When you extract CDS/UTRs, this will ignore the non-coding isoforms - but introns would be taken from both when extracted from TxDb. 
# I think it's best to remove non-coding isoforms, based on all transcript IDs for which there is >=1 CDS annotation (NB coding will include eg NMD transcript biotype)

all_coding_isoforms <- na.omit(unique(pc_genes$transcript_id[pc_genes$type == "CDS"]))

# Create TxDb object (just coding isoforms of pc_genes), then CDS/UTRs, introns can be extracted.

TxDb_pc <- makeTxDbFromGRanges(pc_genes[pc_genes$transcript_id %in% all_coding_isoforms]) # this will not include gene body annotation since no transcript ID, but does not matter for extracting gene segments
# Warning message:
# In .get_cds_IDX(type, phase) :
#  The "phase" metadata column contains non-NA values for features of type stop_codon. This information was ignored.


# extract the different features in different GR objects from TxDb objects

CDS <- cdsBy(TxDb_pc, use.names = T) 
CDS <- unlist(CDS)
mcols(CDS) <- NULL
CDS$Transcript_ID <- names(CDS)
CDS$Gene_ID <- gtf$gene_id[match(CDS$Transcript_ID, gtf$transcript_id)]
CDS$Gene_name <- gtf$gene_name[match(CDS$Transcript_ID, gtf$transcript_id)]

ThreeUTR <- threeUTRsByTranscript(TxDb_pc, use.names=T)
ThreeUTR <- unlist(ThreeUTR) # then names(ThreeUTR) can access the transcript IDs
mcols(ThreeUTR) <- NULL
ThreeUTR$Transcript_ID <- names(ThreeUTR)
ThreeUTR$Gene_ID <- gtf$gene_id[match(ThreeUTR$Transcript_ID, gtf$transcript_id)]
ThreeUTR$Gene_name <- gtf$gene_name[match(ThreeUTR$Transcript_ID, gtf$transcript_id)]

FiveUTR <- fiveUTRsByTranscript(TxDb_pc, use.names=T)
FiveUTR <- unlist(FiveUTR)
mcols(FiveUTR) <- NULL
FiveUTR$Transcript_ID <- names(FiveUTR)
FiveUTR$Gene_ID <- gtf$gene_id[match(FiveUTR$Transcript_ID, gtf$transcript_id)]
FiveUTR$Gene_name <- gtf$gene_name[match(FiveUTR$Transcript_ID, gtf$transcript_id)]

Intron <- intronsByTranscript(TxDb_pc,use.names=T)
Intron <- unlist(Intron)
mcols(Intron) <- NULL
Intron$Transcript_ID <- names(Intron)
Intron$Gene_ID <- gtf$gene_id[match(Intron$Transcript_ID, gtf$transcript_id)]
Intron$Gene_name <- gtf$gene_name[match(Intron$Transcript_ID, gtf$transcript_id)]


# I think we want to only search within genes with at least one xlink (surrogate for expression - if you had expression data could use that instead)
# Then want to get same relative proportions from the different feature types. 
# Elisa had excluded the extended regions around sites from being allowed to be used in the random selection - but I think just leave them in. 
# This gives a truly 'random' selection of sites, as opposed to a comparison of bound vs non-bound - which I think is reasonable

####################################################################################################################################
# Code chunk 4: Import and analyse kmer frequency/zscores for each xlink file in turn
####################################################################################################################################

window_size <- if_else(window_size %% 2 == 0, window_size + 1, window_size)
half_window <- window_size %/% 2

substring_startL <- 1+window_size-k

for(f in files_to_analyse) {
  print(paste("Starting analysis of", f))
  set.seed(123)
  
  dataset_name <- sub("^.*/", "", sub("....$", "", sub(".gz", "", f)))
  
  all_sites <- list()
  
  ##############################################################################################################
  # 4a: Import xlink files for analysis; filter feature GRanges to only include genes with xlink sites
  ##############################################################################################################
  
  read_tsv(f) %>%
    rename_with(~ tolower(.x)) %>%
    filter(!feature %in% c("ncExon", "ncIntron", "Intergenic")) %>%
    filter(chromosome %in% seqnames(BSgenome.Mmusculus.UCSC.mm10)) -> xlink_sites
  if(!is.na(sample_dataset)) {
    if(sample_dataset <= nrow(xlink_sites)){
      xlink_sites %>%
        slice(sample(1:nrow(xlink_sites), sample_dataset)) -> xlink_sites
    } else {
      print("\n")
      print(paste("*** File", f, "has", nrow(xlink_sites), "xlink sites; insufficient to downsample to size", sample_dataset, "; Retaining all ***"))
      print("\n")
    }
  }
  xlink_sites %>%
    mutate(feature = factor(feature, levels = c("CDS", "UTR3", "UTR5", "Intron"))) %>%
    arrange(feature) -> xlink_sites
  feature_types <- as.character(xlink_sites$feature)
  
  if(!is.na(sample_features)) {
    sapply(unique(feature_types), function(x) sum(feature_types == x)) -> feature_counts
    if(is.numeric(sample_features)) {
      names(feature_counts)[which(feature_counts < sample_features)] -> too_short
      if(length(too_short) > 0) {
        print("\n")
        print(paste("***", paste(too_short, collapse = ", "), "features have less than", sample_features, "xlink sites in file", f, "; downsampling only performed for other features ***"))
        print("\n")
      }
    } else {
      min(feature_counts) -> sample_features
    }
  }
      
  if("start" %in% colnames(xlink_sites)) {
    if("end" %in% colnames(xlink_sites)) {
      all_sites[["XLsites"]] <- with(xlink_sites, GRanges(chromosome, IRanges(start, end), strand = strand))
    } else {
      all_sites[["XLsites"]] <- with(xlink_sites, GRanges(chromosome, IRanges(start, start), strand = strand))
    }
  } else if("position" %in%  colnames(xlink_sites)) {
    all_sites[["XLsites"]] <- with(xlink_sites, GRanges(chromosome, IRanges(position, position), strand = strand))
  } else {
    stop("Xlink site position column not identified - header should contain 'start' or 'position'")
  }
  
  CDS_subset <- GPos(CDS[CDS$Gene_name %in% xlink_sites$gene_name])
  ThreeUTR_subset <- GPos(ThreeUTR[ThreeUTR$Gene_name %in% xlink_sites$gene_name])
  FiveUTR_subset <- GPos(FiveUTR[FiveUTR$Gene_name %in% xlink_sites$gene_name])
  Intron_subset <- GPos(Intron[Intron$Gene_name %in% xlink_sites$gene_name])
  
  ######################################################################################################################
  # 4b: Generate 100 sets of randomly sampled sites with the same distribution between feature types as the xlink sites
  ######################################################################################################################
  
  for(i in 1:100) {
    as(c(
      sample(CDS_subset, sum(xlink_sites$feature == "CDS")), 
      sample(ThreeUTR_subset, sum(xlink_sites$feature == "UTR3")), 
      sample(FiveUTR_subset, sum(xlink_sites$feature == "UTR5")),
      sample(Intron_subset, sum(xlink_sites$feature == "Intron"))
    ), "GRanges") -> all_sites[[paste0("sample", i)]]
  }  
    
  ######################################################################################################################
  # 4c: Create extended windows around xl/sampled sites and calculate frequency of kmers
  ######################################################################################################################
  
  mclapply(all_sites, function(x) {
    x %>%
      flank(half_window, both = T) %>%
      resize(window_size) -> sites
    getSeq(BSgenome.Mmusculus.UCSC.mm10, sites, as.character = T) %>%
      lapply(function(w) as_tibble_col(unique(substring(w, 1:substring_startL, k:window_size)), column_name = "kmer")) -> seq_substrings
    seq_substrings %>%
      split(feature_types) -> seq_substrings
    if(!is.na(sample_features)) {
      seq_substrings %>%
        lapply(function(x) sample(x, min(sample_features, length(x)))) -> downsampled_substrings
      seq_substrings %>%
        unlist(recursive = F) %>%
        sample(min(nrow(xlink_sites), sample_features)) -> downsampled_substrings$Gene
      names(downsampled_substrings) <- sapply(names(downsampled_substrings), function(x) paste0(x, "_downsampled"))
      seq_substrings <- c(seq_substrings, downsampled_substrings)
    }
    
    seq_substrings %>%
      lapply(bind_rows) %>%
      bind_rows(.id = "feature")  %>% 
      group_by(kmer, feature) %>% 
        summarise(frequency = n()) %>%
        ungroup()
  }, mc.cores = p) %>%
      bind_rows(.id = "Sample") -> kmer_usage
  
   
  
  ##############################################################################################################
  # 4d: Calculate z-scores and write out
  ##############################################################################################################
  
  kmer_usage %>% 
    pivot_wider(names_from = feature, values_from = frequency) %>% 
    mutate(across(where(is.numeric), ~ if_else(is.na(.x), 0, as.numeric(.x)))) %>%
    # If a kmer is not there in any sample or the xlinks, it will be ignored, otherwise add 0s where not found
    mutate(Gene = CDS+Intron+UTR3+UTR5) %>%
    pivot_longer(-c(Sample, kmer), names_to = "feature", values_to = "Frequency") %>%
    mutate(Type = if_else(Sample == "XLsites", "XLsites", "Random")) %>%
    group_by(Type, kmer, feature) %>%
    summarise(Mean = mean(Frequency), StDev = sd(Frequency)) %>%
    ungroup() %>%
    #mutate(StDev = StDev + 1) %>% # Avoids issue of having s.d. = 0 if no samples have the kmer
    pivot_wider(id_cols = kmer, names_from = c(feature, Type), values_from = c(Mean, StDev)) %>%
    mutate(
      zscore_gene = (Mean_Gene_XLsites - Mean_Gene_Random)/StDev_Gene_Random,
      zscore_3UTR = (Mean_UTR3_XLsites - Mean_UTR3_Random)/StDev_UTR3_Random,
      zscore_CDS = (Mean_CDS_XLsites - Mean_CDS_Random)/StDev_CDS_Random,
      zscore_5UTR = (Mean_UTR5_XLsites - Mean_UTR5_Random)/StDev_UTR5_Random,
      zscore_Intron = (Mean_Intron_XLsites - Mean_Intron_Random)/StDev_Intron_Random
    ) -> kmer_usage
  if(!is.na(sample_features)) {
    kmer_usage %>%
      mutate(
        zscore_gene_downsampled = (Mean_Gene_downsampled_XLsites - Mean_Gene_downsampled_Random)/StDev_Gene_downsampled_Random,
        zscore_3UTR_downsampled = (Mean_UTR3_downsampled_XLsites - Mean_UTR3_downsampled_Random)/StDev_UTR3_downsampled_Random,
        zscore_CDS_downsampled = (Mean_CDS_downsampled_XLsites - Mean_CDS_downsampled_Random)/StDev_CDS_downsampled_Random,
        zscore_5UTR_downsampled = (Mean_UTR5_downsampled_XLsites - Mean_UTR5_downsampled_Random)/StDev_UTR5_downsampled_Random,
        zscore_Intron_downsampled = (Mean_Intron_downsampled_XLsites - Mean_Intron_downsampled_Random)/StDev_Intron_downsampled_Random
      ) -> kmer_usage
    colnames(kmer_usage)[colnames(kmer_usage) == "zscore_gene_downsampled"] <- paste0("zscore_gene_downsampled_", min(sample_features, sum(feature_counts)))
    colnames(kmer_usage)[colnames(kmer_usage) == "zscore_3UTR_downsampled"] <- paste0("zscore_3UTR_downsampled_", min(sample_features, feature_counts["UTR3"]))
    colnames(kmer_usage)[colnames(kmer_usage) == "zscore_CDS_downsampled"] <- paste0("zscore_CDS_downsampled_", min(sample_features, feature_counts["CDS"]))
    colnames(kmer_usage)[colnames(kmer_usage) == "zscore_5UTR_downsampled"] <- paste0("zscore_5UTR_downsampled_", min(sample_features, feature_counts["UTR5"]))
    colnames(kmer_usage)[colnames(kmer_usage) == "zscore_Intron_downsampled"] <- paste0("zscore_Intron_downsampled_", min(sample_features, feature_counts["Intron"]))
  }
  kmer_usage %>%
    select(kmer, starts_with("zscore")) %>%
    arrange(desc(zscore_gene))  %>%
    mutate(across(where(is.numeric), ~ if_else(is.infinite(.x), NA_real_, .x))) %>%
    filter(!grepl("N", kmer)) %>%
    write_tsv(paste0(output_dir, "/", dataset_name, "_", k, "mer_zscores.txt"))
  print(paste("Finished analysis of", f))
}
