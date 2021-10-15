# Microbiome Data Analysis MetaData_GrapevinePhyllo
This repository consit of the meta-data file concerning the informations regarding the samples collected from different grapvevine cultivars
(9 cultivatrs, 3 genetic pools and two organ types), published here https://pubmed.ncbi.nlm.nih.gov/30248973/.
Also illustrate the whole upstream data analysis procedure from FASTQ files to generating phyloseq object in R. Datasets consist 2x250 16S_V4/ITS sequence data.

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
# Learn Error rates
The DADA2 algorithm uses a parametric error model and every amplicon dataset has a different set of error rates. The learnErrors method learns this error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

```{r}
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#To visualize estimated error rates
plotErrors(errF, nominalQ=TRUE) 
#This will show The error rates for each possible transition (A→C, A→G etc.).
If the learned error rates are merging with observed error rates, it is reasonable to go ahead.
```
# Sample Inference

```{r}
dadaFs <- dada(filtFs, err=errF, multithread=T)
dadaFs <- dada(filtRs, err=errR, multithread=T)
#Now we canmerge the forward and reverse reads together to obtain the full denoised sequences. 
#Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, 
#and then constructing the merged “contig” sequences. 
#By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap #region (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

# Construct sequence table
Now, we can construct an amplicon sequence variant table (ASV) table, a higher-resolution version of the OTU table.

```{r}
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=T, verbose=T) #Remove chimeras
```
# Assign taxonomy
 DADA2 uses naive Bayesian classifier method for this purpose. The assignTaxonomy function takes as input a set of sequences to be classified and a training set of reference sequences with known taxonomy, and outputs taxonomic assignments with at least minBoot bootstrap confidence.

```{r}
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/rdp_train_set_14.fa.gz", multithread=T)  
#Dowload a training set from RDP/UNITE databases(reference data) to your directory
```
# Fitting Data to get Phylogenetic tree in data object

```{r}

```

```{r}
library(ape)
fseqs <- getSequences(fseqtab)
names(fseqs) <- fseqs
alignment <- AlignSeqs(DNAStringSet(fseqs), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=T,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#Convert the seqtab into Phyloseq object for data visualization

rownames(taxa)<- paste0("OTU", 1:nrow(taxa))
colnames(seqtab)<- paste0("OTU", 1:ncol(seqtab))
mapfile = read.xls("16s.xlsx", sheet = 1, header = TRUE) # Load Metadata to map your results
all(rownames(seqtab) %in% mapfile$SampleId)
keep.cols <- c("SampleId", "Genetic_Pool", "Region", "Sample_Type", "Cultivar_Name", "Block", "Sampling_Time")
rownames(mapfile) <- mapfile$SampleId
mapfile <- mapfile[rownames(seqtab), keep.cols]

#Now get the phyloseq object-
ps <- phyloseq(tax_table(taxa), sample_data(mapfile),
               otu_table(seqtab, taxa_are_rows = F))
```
# Follow phyloseq tutorial for downstream analysis.
