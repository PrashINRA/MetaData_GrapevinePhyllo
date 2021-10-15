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
S_path <- file.path("data1", "data_16S") #Path to the folder ofthe raw data
filt_path <- file.path("data1", "filtered_16S") # Path to the folder of filtered data
fns <- sort(list.files(ITS_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
# Quality check
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }
#This will generate quality score plote for forward & reverse reads
```
# Trimming & Filtering
For filtering/trimming various tools like Trimmomatic or Cutadapt can also be used but DADA2 also provides an handy function for it
```{r}
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  out[i] = fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]), trimLeft = 10,
                    truncLen=c(240,230),
                    maxN=0, maxEE=c(2,2), truncQ=11, rm.phix=TRUE,
                    n = 2e+06, compress=TRUE, verbose=TRUE)
}

out = data.frame(t(out))
colnames(out) = c("reads.in", "reads.out")
rownames(out)<- basename(filtFs)
#If too few reads are passing the filter, consider increasing maxEE and/or reduce truncQ. If quality drops sharply at the end of your reads, reduce truncLen. If your reads are high quality and you want to reduce computation time in the sample inference step, reduce maxEE.

#if filtering is already done and the filtered files are in folder then reach them by this:
 
 filtS <- sort(list.files(filt_path, full.names = TRUE))
filtFs <- filtS[grepl("R1", filtS)]
filtRs <- filtS[grepl("R2", filtS)]
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 2)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 2)
if(!identical(sample.names, sample.namesR)) stop
names(filtFs) <- sample.names
names(filtRs) <- sample.names

```

