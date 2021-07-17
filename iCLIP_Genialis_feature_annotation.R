#!/usr/bin/env Rscript

######################################################################################
# Annotation of iCLIP data from Genialis, to identify the genes and regions bound
# LM, Oct 2020
######################################################################################

# This script takes as input output files from iMaps (run on Genialis platform, soon moving to Goodwright), and a GTF file (.gtf or .gtf.gz extension) providing the annotation. 

# Files from Genialis can be BED files containing crosslinks, peaks, clusters, paraclu output (.bed or .bed.gz extension). It can also take the 
# .tsv.gz file containing all crosslinks, output from iCount peaks (requires 'chrom', 'position' and 'strand' columns - this could be made more 
# adaptable to take files with start/end).

# Annotated files will be written out as tsv files, with '_annotated.tsv' replacing the original file extension. They will be written to the same directory as the input files.
# One or more files for annotation can be provided, and the script will loop through and annotate each

# For each crosslink site/peak/cluster, the 'most important' gene/feature with which it overlaps is identied and added as columns (see below for full details of hierarchy of genes/features). 
# By default, all transcript isoforms are considered, but those with TSL <= 3 are prioritised.
# --removeTSLNA and --removeTSLover3 flags can be used to exclude isoforms that don't have a TSL (I think includes single exon genes), or have TSL 4 or 5, respectively.
# By default, the importance of features is:
# 3'UTR > CDS > 5'UTR > Intron > non-coding exon > non-coding intron > intergenic
# 'non-coding' category includes pseudogene/non-coding gene biotypes, as well as non-coding isoforms of protein-coding genes (any isoform with a pseudogene transcript biotype, or that does not have a CDS)
# But - for non-coding isoforms of protein coding genes, these will only be reported if none of the protein coding isoforms for that gene overlap with the crosslink/cluster (including in the detailed output)
# The flag --prioritiseCDS can be used to swap CDS and 3'UTR in this hierarchy, so that CDS is considered the most important feature. 

# By default, the output will include only the most 'important' gene/feature, based on the hierarchy of gene/transcript biotypes and feature types.
# Alternatively, the --detailed flag can be used so that in addition to the columns showing the most important gene and feature, additional columns show all genes that overlap with that site/cluster for each type of feature.
# Eg, 'Genes_with_UTR3' shows all the genes with a 3'UTR overlapping the site/Cluster. Both gene ID and name are included separated by '_', whilst different genes are separated by ';'
# This detailed analysis is much more time consuming so is best run on files separately, and with more memory (depending how large the files you have to annotate are).  

# The script can be run from the command line on linux, eg for Babraham cluster:
# module load ssub  ## NOTE this is now necessary before running jobs with ssub
# module load R
# ssub -c 1 --mem 16G --email -o iCLIP_annotation.log Rscript --vanilla path/to/script/iCLIP_Genialis_feature_annotation.R path/to/annotation.gtf.gz bedfile1.bed annotationfile2.tsv.gz subfolder/*.bed.gz --detailed --prioritiseCDS

# If running with --detailed and the files are large, use more memory and only annotate one file at once!
# Sometimes with large files in detailed mode, the script appears to 'fail' - though it seems to do this right at the end, after finishing processing the files. 
# This seems less likely to happen if given more memory, but not entirely predictable! All the tests I ran completed without any error if given 150G memory - for most files less is required.
# If the annotated file has been written out, and the log contains the line 'Finished processing {filename}' - I think it can be assumed to have completed correctly, despite the error message at the end.

# If a GTF file is not provided, or one or more files for annotation (anything not recognised as GTF or a flag), the script will fail with appropriate error message
# If more than one GTF file is provided, only the first will be used (warning message will be generated in the log)
# If one or more file for annotation cannot be imported to generate a GRanges object, it will be skipped, and result in a warning message in the log
# This is likely to arise if it is not in BED format, and does not include column headers 'chrom', 'position' and 'strand' which are currently what is used to convert a tab-delimited text file without bed extension to a GRanges object.


# Alternatively, to run interactively, modify code chunk 1 below to provide the arguments within the script. 

# Requires GenomicRanges, GenomicFeatures, rtracklayer (all available using BiocManager::install) and tidyverse (CRAN install.packages) packages. 

### Hierarchy of genes/features
# Features in the GTF file are ordered based on the type of feature, and the gene/transcript isoform, in the following way:
# (1 done first, 2 used to break ties in that, etc...)
# 1. Based on functional_coding category - all isoforms with protein_coding or IG/TR gene transcript_biotype will appear before other transcript biotypes
# 2. Based on whether the TSL is 3 or below - these isoforms are prioritised over those with TSL 4, 5 or NA
# 3. Based on feature type - Type column factorised based on hierarchy shown above (either UTR3 or CDS first, depending with --prioritiseCDS included as argument)
# 4. Based on TSL - any isoforms with TSL > 3 already excluded; if there are still multiple isoforms, all with the same feature type, the lowest TSL will be prioritised
# 5. Based on CCDS - those with a CCDS will be ranked first



library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(tidyverse)

###########################################################################################
# Code chunk 1: Parse arguments and check both GTF file and file(s) to annotate provided
###########################################################################################

Args <- commandArgs(trailingOnly=TRUE)

### If running interactively, comment out the commandArgs line above.
### The below line will then find files with .gtf or .gtf.gz extension, and files with .bed or .bed.gz extension, in the working directory, and use these as arguments. 
### If another format eg .txt or .tsv is used for the files to be annotated, this can be modified, or the files simply listed.
### Alternatively, custom file paths for the gtf file and files to be annotated can be assigned to Args.
### Also need to add the --prioritiseCDS or --detailed arguments to Args if these are required.
# Args <- c(list.files(pattern = "\\.gtf$"), list.files(pattern = "\\.gtf\\.gz$"), list.files(pattern = "\\.bed$"), list.files(pattern = "\\.bed\\.gz$"))

gtf_file <- Args[grepl("gtf$|gtf.gz$", Args)]
if(length(gtf_file) > 1) {
  warning("More than one GTF file provided; using only first listed")
  gtf_file <- gtf_file[1]
}
if(length(gtf_file) == 0) {
  stop("No GTF file provided")
} 

print("GTF file used for annotation:")
print(gtf_file)

if("--prioritiseCDS" %in% Args) {
  featureHierarchy <- c("CDS", "UTR3", "UTR5", "Intron", "ncExon", "ncIntron")
} else{
  featureHierarchy <- c("UTR3", "CDS", "UTR5", "Intron", "ncExon", "ncIntron")
}

print("Features will be prioritised in the following order:")
print(paste(c(featureHierarchy, "Intergenic"), collapse = " > "))

files_to_annotate <- Args[!(Args %in% c("--prioritiseCDS", "--detailed") | grepl("gtf$|gtf.gz$", Args))]

if(length(files_to_annotate) == 0) {
  stop("No files provided for annotation")
}

print("The following files will be annotated:")
print(files_to_annotate)

###################################################################################################################
# Code chunk 2: Import GTF file, extract features and make overall GRanges with features organised hierarchically
###################################################################################################################

gtf <- import(gtf_file)

# Some GTF file eg from Gencode have 'gene_type' and 'transcript_type' instead of 'gene/transcript_biotype' 
# Also have 'ccdsid' rather than 'ccds_id'
# So we need to rename these if that's the case to ensure the rest of the code works
names(mcols(gtf))[names(mcols(gtf)) == "transcript_type"] <- "transcript_biotype"
names(mcols(gtf))[names(mcols(gtf)) == "gene_type"] <- "gene_biotype"
names(mcols(gtf))[names(mcols(gtf)) == "ccdsid"] <- "ccds_id"

# For this version of the script, I am not going to eliminate things with TSL > 3
# transcript_support_level in the gtf file (at least v101 which I'm using to test..!) contains many annotations where the TSL is followed by a string in brackets, eg:
# "1 (assigned to previous version 13)". So first need to strip these out, before converting to numeric and selecting only those <=3
gtf$transcript_support_level <- as.numeric(sub(" .*$", "", gtf$transcript_support_level))

if("--removeTSLNA" %in% Args) {
  gtf <- gtf[!is.na(gtf$transcript_support_level)]
}
if("--removeTSLover3" %in% Args) {
  gtf <- gtf[gtf$transcript_support_level <=3]
}

# Including anything that is not a protein coding gene in the ncRNA category, this includes pseudogenes (NB, some pseudogene biotypes include protein coding isoforms - these are excluded here and included in the protein coding category instead. Also, protein coding genes include pseudogene transcript isoforms which are included here rather than in protein coding category)
ncRNA_genes <- gtf[(gtf$gene_biotype %in% c("3prime_overlapping_ncRNA","bidirectional_promoter_lncRNA","lincRNA","macro_lncRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","rRNA","scaRNA","scRNA","snoRNA","snRNA","sRNA", "ribozyme", "antisense", "processed_transcript", "sense_intronic", "sense_overlapping") | grepl("pseudogene", gtf$gene_biotype) | grepl("pseudogene", gtf$transcript_biotype)) & !(gtf$transcript_biotype %in% "protein_coding")]

# Include immunoglobulin & Tcr genes in protein coding ones. Exclude anything with pseudogene in transcript_biotype, but include any isoforms with protein_coding as transcript_biotype (can be designated pseudogenes if you just look at gene_biotype)
pc_genes <- gtf[(gtf$gene_biotype=="protein_coding"| grepl("^[IT][GR]_.*_gene$",gtf$gene_biotype) | gtf$transcript_biotype %in% "protein_coding") & !(grepl("pseudogene", gtf$transcript_biotype))]

## so this is taking everything that has a GENE biotype of protein_coding, and all transcripts from that gene regardless of whether they are coding or not (unless transcript biotype includes pseudogene). 
# When you extract CDS/UTRs, this will ignore the non-coding isoforms - but introns would be taken from both when extracted from TxDb. 
# I think it's best to separate these into coding and non-coding isoforms, based on all transcript IDs for which there is >=1 CDS annotation (NB coding will include eg NMD transcript biotype)

all_coding_isoforms <- na.omit(unique(pc_genes$transcript_id[pc_genes$type == "CDS"]))

# Create TxDb objects for each (splitting pc_genes based on coding isoforms above), then CDS/UTRs, exons, introns can be extracted from each.

TxDb_pc <- makeTxDbFromGRanges(pc_genes[pc_genes$transcript_id %in% all_coding_isoforms]) # this will not include gene body annotation since no transcript ID, but does not matter for extracting gene segments
# Warning message:
# In .get_cds_IDX(type, phase) :
#  The "phase" metadata column contains non-NA values for features of type stop_codon. This information was ignored.

TxDb_nc <- makeTxDbFromGRanges(c(ncRNA_genes, pc_genes[!pc_genes$transcript_id %in% all_coding_isoforms]))
# this now combines both ncRNAs and non-coding isoforms of protein coding genes 

# extract the different features in different GR objects from TxDb objects

CDS <- cdsBy(TxDb_pc, use.names = T) 
CDS <- unlist(CDS)
mcols(CDS) <- NULL
mcols(CDS)$Type = "CDS"
mcols(CDS)$Transcript_ID <- names(CDS)

ThreeUTR <- threeUTRsByTranscript(TxDb_pc, use.names=T)
ThreeUTR <- unlist(ThreeUTR) # then names(ThreeUTR) can access the transcript IDs
mcols(ThreeUTR) <- NULL
mcols(ThreeUTR)$Type <- "UTR3"
mcols(ThreeUTR)$Transcript_ID <- names(ThreeUTR)

FiveUTR <- fiveUTRsByTranscript(TxDb_pc, use.names=T)
FiveUTR <- unlist(FiveUTR)
mcols(FiveUTR) <- NULL
mcols(FiveUTR)$Type <- "UTR5"
mcols(FiveUTR)$Transcript_ID <- names(FiveUTR)

Intron <- intronsByTranscript(TxDb_pc,use.names=T)
Intron <- unlist(Intron)
mcols(Intron) <- NULL
mcols(Intron)$Type <- "Intron"
mcols(Intron)$Transcript_ID <- names(Intron)

# For ncRNA and pseudogenes, extract just exons and introns (a few pgs have CDS/UTRs, but ignoring these and just taking exons)
ncExon <- exonsBy(TxDb_nc, use.names = T)
ncExon <- unlist(ncExon)
mcols(ncExon) <- NULL
mcols(ncExon)$Type <- "ncExon"
mcols(ncExon)$Transcript_ID <- names(ncExon)


ncIntron <- intronsByTranscript(TxDb_nc, use.names = T)
ncIntron <- unlist(ncIntron)
mcols(ncIntron) <- NULL
mcols(ncIntron)$Type <- "ncIntron"
mcols(ncIntron)$Transcript_ID <- names(ncIntron)


# Join all of these together
all_features <- c(CDS, ThreeUTR, FiveUTR, Intron, ncExon, ncIntron)

## Merge in other info we want to include for ranking features, allowing the GRanges object to be ranked hierarchically
# This will mean the first overlap found will be the one we want to list as the main hit in file we are annotating

# CCDS - genes which have one of these are typically more confident isoforms. Make a column that is TRUE if CCDS ID is NA, so when ordered by this those with a CCDS ID (ie F) end up at the top
all_features$no_ccds <- is.na(gtf$ccds_id)[match(all_features$Transcript_ID, gtf$transcript_id)]

# Transcript support level - 1 is most confident, 5 least; these have already been converted to numeric, but none excluded in this version of the script
all_features$TSL <- gtf$transcript_support_level[match(all_features$Transcript_ID, gtf$transcript_id)]
# since we're not excluding TSL > 3 here, I want to add a column that will allow to order based on whether or not TSL > 3, so those with low (high confidence) TSL are prioritised early on
all_features$TSL45NA <- all_features$TSL > 3 | is.na(all_features$TSL)

# Transcript biotype - essentially here I just want to know whether it is coding and functional. Ie, want to decrease the importance of things like NMD isoforms (otherwise these may have a much longer 3'UTR than other more 'normal' isoforms, so it may look like 3'UTR binding when it should be classified as CDS)
# So test whether it is either protein_coding, or one of the IG/TR gene segments
# but simplest to first add in the transcript biotype to all_features
all_features$transcript_biotype <- gtf$transcript_biotype[match(all_features$Transcript_ID, gtf$transcript_id)]
all_features$functional_coding <- if_else(all_features$transcript_biotype %in% "protein_coding" | grepl("^[IT][GR]_.*_gene$", all_features$transcript_biotype), "coding", "noncoding")

# Also merge in the gene name and ID, since I these need to be added into the file we are annotating once we've found the overlaps
all_features$Gene_ID <- gtf$gene_id[match(all_features$Transcript_ID, gtf$transcript_id)]
all_features$Gene_name <- gtf$gene_name[match(all_features$Transcript_ID, gtf$transcript_id)]

# Order in the following way (1 done first, 2 used to break ties in that, etc...):
# 1. Based on functional_coding category - so all isoforms with protein_coding or IG/TR gene transcript_biotype will appear before other transcript biotypes
# 2. TSL_over_3 column so those with high (low confidence) TSL or NA are moved to the end
# 2. Based on feature type - Type column factorised based on hierarchy established above (either UTR3 or CDS first, depending with --prioritiseCDS included as argument)
all_features$Type <- factor(all_features$Type, levels = featureHierarchy)
# 3. Based on TSL - the lowest TSL will be prioritised
# 4. Based on CCDS - those with a CCDS will be ranked first

all_features_ordered <- all_features[order(all_features$functional_coding, all_features$TSL45NA, all_features$Type, all_features$TSL, all_features$no_ccds)]


###############################################################################################################################
# Code chunk 3: Import file to be annotated - this and later chunks done in a loop through all provided files for annotation
###############################################################################################################################

for(f in files_to_annotate) {
  print(paste("Processing", f))
  if(grepl("bed$|bed.gz$", f)) {
    try(ann_file <- import(f))
  } else {
    try(read_tsv(f) %>% with(GRanges(chrom, IRanges(position, position), strand = strand, mcols = dplyr::select(., name, group_id, score, score_extended, FDR))) -> ann_file)
  }
  if(!exists("ann_file")) {
    warning(paste(f, "file for annotation could not be read; skipping this file. Please provide either bed/bed.gz format, or a tab delimited text format file with column names 'chrom', 'position' and 'strand' designating the chromosome, crosslink position and strand respectively (eg, the tsv file for all crosslinks following iCount peak calling)")) 
    next }
  
  # Output from Genialis generally has 'chr' preceding chromosome names, and 'chrM' rather than 'MT' for mitochondrial DNA, whereas GTF does not
  # But might be safest to assume this could be the case for one or both, and remove chr from both (which doesn't matter if it's not there), and convert '^M$' to 'MT' after that.
  # (NB realised this could be simplified using seqlevelsStyle())
  seqlevels(all_features_ordered) <- sub('chr','',seqlevels(all_features_ordered)) 
  seqlevels(all_features_ordered) <- sub('^M$','MT',seqlevels(all_features_ordered))
  
  ann_file2 <- ann_file # create new version used for finding overlaps where we have made chromosome names consistent, but leave them as they were in the orignal file which is what will be exported
  seqlevels(ann_file2) <- sub('chr','',seqlevels(ann_file2)) 
  seqlevels(ann_file2) <- sub('^M$','MT',seqlevels(ann_file2))
  
  dropSeqlevels(ann_file, c(paste0("chr", seqlevels(ann_file2)[!seqlevels(ann_file2) %in% seqlevels(all_features_ordered)]), sub("^MT$", "M", seqlevels(ann_file2)[!seqlevels(ann_file2) %in% seqlevels(all_features_ordered)])), pruning.mode = "coarse") -> ann_file # must do this before ann_file2 or the level(s) to remove will no longer be there!
  dropSeqlevels(ann_file2, seqlevels(ann_file2)[!seqlevels(ann_file2) %in% seqlevels(all_features_ordered)], pruning.mode = "coarse") -> ann_file2
  
  
#######################################################################################################################
# Code chunk 3: Overlap crosslinks/clusters from file to be annotated, with the ranked features derived from the GTF
#######################################################################################################################

  # By default this will ONLY output the 'most likely' feature (based on the hierarchy in the all_features_ordered GRanges object)
  # Instead can use the --detailed flag to give both the 'most likely' feature, and lists of all possible genes for each feature type - this will take much longer!
  
  
  RegionsBound <- findOverlaps(ann_file2, all_features_ordered, select = "first")
  ann_file$Gene_name <- all_features_ordered$Gene_name[RegionsBound]
  ann_file$Gene_ID <- all_features_ordered$Gene_ID[RegionsBound]
  ann_file$Feature <- as.character(all_features_ordered$Type[RegionsBound])
  ann_file$Feature[is.na(ann_file$Feature)] <- "Intergenic"
  if("--detailed" %in% Args) {
    ann_file$Genes_with_CDS <- NA_character_
    ann_file$Genes_with_UTR3 <- NA_character_
    ann_file$Genes_with_UTR5 <- NA_character_
    ann_file$Genes_with_Intron <- NA_character_
    ann_file$Nc_genes_with_exon <- NA_character_
    ann_file$Nc_genes_with_intron <- NA_character_
    all_features_ordered[all_features_ordered$Type == "CDS"] -> CDS_features
    as_tibble(CDS_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> CDS_tb
    CDSBound <- findOverlaps(ann_file2, CDS_features)
    for(i in unique(CDSBound@from)) {
      CDS_tb %>%
        slice(CDSBound@to[CDSBound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        deframe() %>%
        paste(collapse = ";") -> ann_file$Genes_with_CDS[i]
    }
    all_features_ordered[all_features_ordered$Type == "UTR3"] -> UTR3_features
    as_tibble(UTR3_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> UTR3_tb
    UTR3Bound <- findOverlaps(ann_file2, UTR3_features)
    for(i in unique(UTR3Bound@from)) {
      UTR3_tb %>%
        slice(UTR3Bound@to[UTR3Bound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        deframe() %>%
        paste(collapse = ";") -> ann_file$Genes_with_UTR3[i]
    }
    all_features_ordered[all_features_ordered$Type == "UTR5"] -> UTR5_features
    as_tibble(UTR5_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> UTR5_tb
    UTR5Bound <- findOverlaps(ann_file2, UTR5_features)
    for(i in unique(UTR5Bound@from)) {
      UTR5_tb %>%
        slice(UTR5Bound@to[UTR5Bound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        deframe() %>%
        paste(collapse = ";") -> ann_file$Genes_with_UTR5[i]
    }
    all_features_ordered[all_features_ordered$Type == "Intron"] -> Intron_features
    as_tibble(Intron_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> Intron_tb
    IntronBound <- findOverlaps(ann_file2, Intron_features)
    for(i in unique(IntronBound@from)) {
      Intron_tb %>%
        slice(IntronBound@to[IntronBound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        deframe() %>%
        paste(collapse = ";") -> ann_file$Genes_with_Intron[i]
    }
    all_features_ordered[all_features_ordered$Type == "ncExon"] -> ncExon_features
    as_tibble(ncExon_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> ncExon_tb
    ncExonBound <- findOverlaps(ann_file2, ncExon_features)
    for(i in unique(ncExonBound@from)) {
      overlapping_coding_isoforms <- na.omit(c(ann_file$Genes_with_CDS[i], ann_file$Genes_with_UTR3[i], ann_file$Genes_with_UTR5[i], ann_file$Genes_with_Intron[i]))
      if(length(overlapping_coding_isoforms) > 0) {
        paste(overlapping_coding_isoforms, collapse = ";") -> overlapping_coding_isoforms
        unique(unlist(strsplit(overlapping_coding_isoforms, ";"))) -> overlapping_coding_isoforms
      }
      ncExon_tb %>%
        slice(ncExonBound@to[ncExonBound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        filter(!Gene_name_ID %in% overlapping_coding_isoforms) %>%
        deframe() %>%
        paste(collapse = ";") -> ncExon_genes
      if(ncExon_genes != "") {
        ncExon_genes -> ann_file$Nc_genes_with_exon[i]
      }
    }
    all_features_ordered[all_features_ordered$Type == "ncIntron"] -> ncIntron_features
    as_tibble(ncIntron_features) %>%
      unite("Gene_name_ID", Gene_ID, Gene_name) -> ncIntron_tb
    ncIntronBound <- findOverlaps(ann_file2, ncIntron_features)
    for(i in unique(ncIntronBound@from)) {
      overlapping_coding_isoforms <- na.omit(c(ann_file$Genes_with_CDS[i], ann_file$Genes_with_UTR3[i], ann_file$Genes_with_UTR5[i], ann_file$Genes_with_Intron[i]))
      if(length(overlapping_coding_isoforms) > 0) {
        paste(overlapping_coding_isoforms, collapse = ";") -> overlapping_coding_isoforms
        unique(unlist(strsplit(overlapping_coding_isoforms, ";"))) -> overlapping_coding_isoforms
      }
      ncIntron_tb %>%
        slice(ncIntronBound@to[ncIntronBound@from == i]) %>%
        distinct(Gene_name_ID) %>%
        filter(!Gene_name_ID %in% overlapping_coding_isoforms) %>%
        deframe() %>%
        paste(collapse = ";") -> ncIntron_genes
      if(ncIntron_genes != "") {
        ncIntron_genes -> ann_file$Nc_genes_with_intron[i]
      }
    } 
  }
  
#################################################################################################
# Code chunk 4: Write out annotated file in tsv format
#################################################################################################
  
  output_filename <- sub("\\.[a-zA-Z]+?$", "_annotated.tsv", sub(".gz$","", f))
  
  ann_file %>%
    as_tibble()  %>% 
    rename(chromosome = seqnames) %>%
    rename_all(function(x) sub("^mcols.", "", x)) %>%
    write_tsv(output_filename)
  print(paste("Finished processing", f))
}


#################################################################################################
# Code chunk 5: sessionInfo()
#################################################################################################

sessionInfo()

