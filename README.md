# Process_CLIP_data
A collection of R scripts used to process iCLIP or other types of CLIP data

## iCLIP_Genialis_feature_annotation.R

This script takes data that has been through the iCount pipeline (available on the iMaps platform from Jernej Ule's lab).

It requires the tidyverse package (available through CRAN with install.packages("tidyverse")), and the rtracklayer, GenomicRanges and GenomicFeatures packages, available through Bioconductor (BiocManager::install("package_name")).

Whilst some of the output files from iMaps include the genes to which each crosslink site maps, the features are not annotated. This script reannotates either the BED files containing only significant crosslink sites, or the tsv files containing all crosslink sites (with FDR annotated), together with a GTF file, and assigns the most likely gene/feature bound at that site. This is done in a hierarchical way - full details provided in annotation of the script. 

The script gives the flexibility to prioritise either 3'UTR (default) or CDS (--prioritise_CDS flag) in the hierarchy.

It also gives the option of running in a more detailed mode (--detailed), which will identify and report all possible genes and features to which the crosslink site could be assigned - note this will take a lot longer to run! 

The script can be tested using the Mouse_GRCm38.90_chr3_200kb.gtf and three iCLIP...chr3_200kb.txt files as input; these files include data covering a 200kb region of chromosome 3 for replicate iCLIP datasets. The test can be run with default parameters by setting the working directory to test_data and replacing the line:
```
Args <- commandArgs(trailingOnly=TRUE)
```
with:
```
Args <- c(list.files(pattern = "\\.gtf$"), list.files(pattern = "\\.txt$"))
```

The expected outputs can be found in the iCLIP...chr3_200kb_annotated.tsv files. 



## Combine_replicate_datasets.R

This script takes as input a list of replicate datasets, and outputs the crosslink sites that identified as significant in at least a user-specified number of replicates. 

It requires the tidyverse package, available through CRAN (install.packages("tidyverse"))

If an 'FDR' column is present in the data, it will be filtered to include only sites with FDR <= 0.05 (unless this value is amended by the user), otherwise, all sites in the dataset are assumed to be significant. 

The input data is assumed to come from the iCLIP_Genialis_feature_annotation.R script, and so the column names must match the output from this. 

The 'start' and 'end' columns are summarised into a single 'position' column (since they are always identical for crosslink sites), and sum and mean of the score column are calculated (although only including replicates for which the site was significant).

The filtered data including only the sites that are significant in at least the specified number of replicates is saved out as a tab-delimited text file with the user-specified name.

As written, the inputs specified are to take the three iCLIP...chr3_200kb_annotated.tsv files (in test_data) as input (with working directory set to test_data), and report only sites that are significant (FDR <= 0.05) in at least 2 replicates. The expected output file (iCLIP_combined_replicates.txt) is provided.


## kmer_analysis.R

This script takes a gtf file together with annotated iCLIP datasets (that have been through iCLIP_Genialis_feature_annotation.R above), and calculates z-scores for kmers of length specified by the user. Ideally, crosslink sites should be limited to only those that are significant, or can be the output from Combine_replicate_datasets.R. 

It requires the tidyverse and parallel packages (available through CRAN with install.packages("package_name")), and the rtracklayer, BSgenome.Mmusculus.UCSC.mm10, GenomicRanges and GenomicFeatures packages, available through Bioconductor (BiocManager::install("package_name")).

Z-scores are calculated for the whole gene and for each feature type (3'UTR, CDS, 5'UTR and intron). Sites mapping to non-coding genes or intergenic regions are excluded. By default all sites are analysed, but subsampling to ensure either that an equal number of sites are analysed between different datasets, or between different features within the same datasets, can be done - see details in the script annotation.

The value of k is provided by the user (default 5).

For full datasets, it is recommended to run with parallel processing on a computer cluster. However, the script can be tested on the iCLIP_combined_replicates.txt file, with Mouse_GRCm38.90_chr3_200kb.gtf, by setting the working directory to test_data and replacing the line:

```
Args <- commandArgs(trailingOnly=TRUE)
```
with:
```
Args <- c("iCLIP_combined_replicates.txt", "Mouse_GRCm38.90_chr3_200kb.gtf")
```

The expected output when run with default parameters (k = 5 and no subsampling) can be found in the iCLIP_combined_replicates_5mer_zscores.txt file. Note that due to the sparsity of the data, this test results in relatively low z scores and many NA values.
