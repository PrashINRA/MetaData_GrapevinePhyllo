# Microbiome Data Analysis MetaData_GrapevinePhyllo
This repository consit of the meta-data file concerning the informations regarding the samples collected from different grapvevine cultivars
(9 cultivatrs, 3 genetic pools and two organ types), published here https://pubmed.ncbi.nlm.nih.gov/30248973/.
Also illustrate the whole data analysis procedure from FASTQ files to Unsupervised clustering and Hypothesis testing in R.

# Load Packages:
```{r}
pkgs <- c("ggplot2", "gridExtra", "gdata", "ape", "phyloseq", "dada2")
sapply(pkgs, require, character.only = TRUE)
theme_set(theme_bw())
set.seed(100)
```
# Read in your data (FastQ files with forward(R1) and reverse(R2) reads) from you working directory:
```{r}
#Set paths like:
S_path <- file.path("data1", "data_16S")
filt_path <- file.path("data1", "filtered_16S") 
fns <- sort(list.files(ITS_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
# Quality check
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }
#This will generate quality score plote for forward & reverse reads
```
