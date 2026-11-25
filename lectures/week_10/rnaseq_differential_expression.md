# RNA-seq Part 3: Differential Expression on the Cluster

This week we finish the expression portion of the [RNA-seq Bioinformatics course](https://rnabio.org/course/) (rnabio.org). Starting from last week's StringTie estimates and htseq-count read counts, we'll test for genes that are differentially expressed (DE) between the UHR and HBR samples using three different tools (Ballgown, edgeR, and DESeq2), compare their results, and visualize them.

By the end of this exercise you should be able to:

- Run a differential expression analysis with Ballgown (from StringTie FPKM estimates) and with edgeR and DESeq2 (from raw read counts)
- Explain why log fold change shrinkage is used, and what it does
- Compare the DE gene lists produced by different methods, and check them against the known ERCC spike-in fold changes
- Make and interpret common DE plots: MA plots, sample distance heatmaps, MDS plots, and volcano plots

> [!NOTE]
> As always, we strongly suggest typing the shell commands out yourself rather than copying and pasting. Most of this week is R code, though, and much of it is long. Feel free to paste the R blocks, but read them first so you know what each piece does, and run them a few lines at a time so you can look at the output as you go.

This page builds directly on [last week's assignment](../week_09/rnaseq_expression.md). All of the commands have already been adjusted to work on the cluster, so **follow the commands on this page, rather than the ones on rnabio.org.** Each section heading links to the original lesson if you'd like to see it.

## Contents

1. [Before you start](#1-before-you-start)
2. [Differential expression with Ballgown](#2-differential-expression-with-ballgown)
3. [Differential expression with edgeR](#3-differential-expression-with-edger)
4. [Differential expression with DESeq2](#4-differential-expression-with-deseq2)
5. [Compare DE genes across methods](#5-compare-de-genes-across-methods)
6. [DE visualization with Ballgown](#6-de-visualization-with-ballgown)
7. [DE visualization with DESeq2](#7-de-visualization-with-deseq2)
8. [DE visualization advanced](#8-de-visualization-advanced)
9. [For-credit deliverable](#9-for-credit-deliverable)
10. [Cleaning up](#10-cleaning-up)
11. [Cheat sheet](#11-cheat-sheet)

---

## 1. Before you start

### Launch the container

Log in to Compute2, and start an interactive job in the same RNA-seq toolkit container as the last two weeks, with your workspace mounted at `/workspace`:

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu

srun -A compute2-workshop -p workshop -c 2 --mem=8G --pty \
  --container-image=griffithlab/rnaseq-toolkit:1.1.1 \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch/$USER/rnabio-workspace:/workspace" \
  --container-workdir=/workspace \
  /bin/bash
```

Then load the course environment variables:

```bash
source /workspace/rnabio_env.sh
echo $RNA_HOME
```

That should print `/workspace/rnaseq`. If it prints nothing, or `source` complains that the file doesn't exist, go back and revisit [Section 4 of the week 8 assignment](../week_08/rnaseq_alignment.md#4-set-up-the-course-environment-variables).

> [!IMPORTANT]
> Unless otherwise noted, **all of the commands for the rest of this assignment are run inside this container**. If you get disconnected or `exit`, re-run the `srun` command above and `source /workspace/rnabio_env.sh` again.

### Check that last week's results are there

Everything this week starts from the expression results we made last week. Make sure they exist:

```bash
ls -d $RNA_HOME/expression/stringtie/ref_only/*_Rep*
ls -l $RNA_HOME/expression/stringtie/ref_only/gene_tpm_all_samples.tsv
ls -l $RNA_HOME/expression/htseq_counts/gene_read_counts_table_all_final.tsv
ls -l $RNA_HOME/de/htseq_counts/ENSG_ID2Name.txt
```

You should see six StringTie directories (`HBR_Rep1` through `UHR_Rep3`) and the three files. If any of these are missing, finish [last week's assignment](../week_09/rnaseq_expression.md) first.

### Running R inside the container

Almost everything this week happens in R, which is installed in the container along with all of the packages we need. You start it by typing `R` at the container's shell prompt. Your prompt will change to a `>`, which means you're now typing R commands rather than shell commands. To get back to the shell, type `quit(save = "no")`.

### Viewing PDFs and downloading files

Most of this week's results are plots saved as PDF files, and a few steps ask you to open a results table in a spreadsheet or paste a gene list into a website. To do this, mount your storage on your laptop the same way we did in week 8 (see [Viewing HTML reports](../week_08/rnaseq_alignment.md#viewing-html-reports)), and open the files from the `rnabio-workspace/rnaseq/` folder. You can copy files from there to your laptop, too.

---

## 2. [Differential expression with Ballgown](https://rnabio.org/module-03-expression/0003/03/01/Differential_Expression-Ballgown/)

![RNA-seq_Flowchart4](https://rnabio.org/assets/module_3/RNA-seq_Flowchart4.png)

### Ballgown DE analysis

First, let's use Ballgown to compare the UHR and HBR conditions. Refer to the Ballgown manual for a more detailed explanation:

- [https://www.bioconductor.org/packages/release/bioc/html/ballgown.html](https://www.bioconductor.org/packages/release/bioc/html/ballgown.html)

Create and change to the Ballgown ref-only results directory:

```bash
mkdir -p $RNA_HOME/de/ballgown/ref_only/
cd $RNA_HOME/de/ballgown/ref_only/
```

Perform the UHR vs. HBR comparison, using all replicates, for known (reference only mode) transcripts.

First, start an R session:

```bash
R
```

Run the following R commands in your R session.

```R
# load the required libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# Create phenotype data needed for ballgown analysis
ids = c("UHR_Rep1", "UHR_Rep2", "UHR_Rep3", "HBR_Rep1", "HBR_Rep2", "HBR_Rep3")
condition = c("UHR", "UHR", "UHR", "HBR", "HBR", "HBR")
inputs = "/workspace/rnaseq/expression/stringtie/ref_only/"
path = paste(inputs, ids, sep="")
pheno_data = data.frame(ids, condition, path)

# Load ballgown data structure and save it to a variable "bg"
bg = ballgown(samples = as.vector(pheno_data$path), pData = pheno_data)

# View the current order of the comparison: UHR/HBR (fold-changes > 1 means higher expression in UHR than HBR)
factor(pData(bg)$condition) # "Levels: HBR UHR", alphabetically by default determines the comparison

# Note that if you want to reverse the comparison of HBR and UHR. We would achieve this by changing the order of the "levels":
pData(bg)$condition <- factor(pData(bg)$condition, levels = c("UHR","HBR"))
factor(pData(bg)$condition) # "Levels: UHR HBR"

# Lets change the condition order back for consistency
pData(bg)$condition <- factor(pData(bg)$condition, levels = c("HBR","UHR"))
factor(pData(bg)$condition) # "Levels: HBR UHR"

# Display a description of this object
bg

# Load all attributes including gene name
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])
bg_transcript_names = unique(bg_table[, c(1, 6)])

# Save the ballgown object to a file for later use
save(bg, file = 'bg.rda')

# Pull the gene and transcript expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))
transcript_expression = as.data.frame(texpr(bg))

# Perform differential expression (DE) analysis with no filtering, at the transcript level
results_transcripts = stattest(bg, feature = "transcript", covariate = "condition", getFC = TRUE, meas = "FPKM")

# Then add on transcript/gene names and sample level fpkm values for export
results_transcripts = merge(results_transcripts, bg_transcript_names, by.x = c("id"), by.y = c("t_id"))
results_transcripts = merge(results_transcripts, transcript_expression, by.x = c("id"), by.y = c("row.names"))

# Perform differential expression (DE) analysis with no filtering, at the gene level
results_genes = stattest(bg, feature = "gene", covariate = "condition", getFC = TRUE, meas = "FPKM")

# Then add on transcript/gene names and sample level fpkm values for export
results_genes = merge(results_genes, bg_gene_names, by.x = c("id"), by.y = c("gene_id"))
results_genes = merge(results_genes, gene_expression, by.x = c("id"), by.y = c("row.names"))

# Note that in the statistical tests above, the design is very simple with condition as the single covariate.
# More complex design are possible including time series ("timecourse" parameter), confounders (using adjustvars parameter), etc.
# The ballgown vignettes describe these options in detail.

# Save a tab delimited file for both the transcript and gene results
write.table(results_transcripts, "UHR_vs_HBR_transcript_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset (bg, "SparseArray::rowVars(texpr(bg)) > 1", genomesubset = TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])
bg_filt_transcript_names = unique(bg_filt_table[, c(1,6)])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature = "transcript", covariate = "condition", getFC = TRUE, meas = "FPKM")
results_transcripts = merge(results_transcripts, bg_filt_transcript_names, by.x = c("id"), by.y = c("t_id"))
results_transcripts = merge(results_transcripts, transcript_expression, by.x = c("id"), by.y = c("row.names"))

results_genes = stattest(bg_filt, feature = "gene", covariate = "condition", getFC = TRUE, meas = "FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x = c("id"), by.y = c("gene_id"))
results_genes = merge(results_genes, gene_expression, by.x = c("id"), by.y = c("row.names"))

# Output the filtered list of genes and transcripts and save to tab delimited files
write.table(results_transcripts, "UHR_vs_HBR_transcript_results_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_genes, "UHR_vs_HBR_gene_results_filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)

# Output the significant gene results to a pair of tab delimited files
write.table(sig_transcripts, "UHR_vs_HBR_transcript_results_sig.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig_genes, "UHR_vs_HBR_gene_results_sig.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# Exit the R session
quit(save = "no")
```

Once you have completed the Ballgown analysis in R, exit the R session and continue with the steps below. A copy of the above R code is located [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_Part1_ballgown.R) (but note that it uses the course's `/home/ubuntu/...` paths, not ours).

What does the raw output from Ballgown look like?

```bash
head UHR_vs_HBR_gene_results.tsv
```

How many genes are there on this chromosome?

```bash
grep -v feature UHR_vs_HBR_gene_results.tsv | wc -l
```

How many passed filter in UHR or HBR?

```bash
grep -v feature UHR_vs_HBR_gene_results_filtered.tsv | wc -l
```

How many differentially expressed genes were found on this chromosome (p-value < 0.05)?

```bash
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | wc -l
```

Display the top 20 DE genes. In the following commands we use `grep -v feature` to remove lines that contain "feature". Then we use `sort` to sort the data in various ways. The `k` option specifies that we want to sort on a particular column (`3` in this case, which has the DE fold change values). The `n` option tells `sort` to sort numerically. The `r` option tells `sort` to reverse the sort.

```bash
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -rnk 3 | head -n 20 | column -t #Higher abundance in UHR
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | sort -nk 3 | head -n 20 | column -t #Higher abundance in HBR
```

Save all genes with P<0.05 to a new file.

```bash
grep -v feature UHR_vs_HBR_gene_results_sig.tsv | cut -f 6 | sed 's/\"//g' > DE_sig_genes_ballgown.tsv
head DE_sig_genes_ballgown.tsv
```

### Practical exercise 9

Assignment: Use Ballgown to identify differentially expressed genes from the StringTie expression estimates (i.e., Ballgown table files) which you created in Practical Exercise 8 last week.

- Hint: Follow the example R code above.
- Hint: You will need to change how the `pheno_data` object is created to point to the correct sample ids, condition, and path to your inputs (the StringTie results files). Remember to use `/workspace/rnaseq/...` paths.
- Hint: Make sure to save your Ballgown data object to file (e.g., `bg.rda`) for use in subsequent practical exercises.
- Hint: You may wish to save both a complete list of genes with differential expression results as well as a subset which are filtered and pass a significance test.

Solution: When you are ready you can check your approach against the [Solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-9---differential-expression). (The solution uses `/home/ubuntu/workspace/rnaseq/...`; change that to `/workspace/rnaseq/...`.)

### ERCC DE analysis

This section will compare the differential expression estimates obtained by the RNA-seq analysis against the expected differential expression results for the ERCC spike-in RNAs (mix1-UHR vs mix2-HBR).

First set up a directory to store the results of this analysis:

```bash
mkdir $RNA_HOME/de/ercc_spikein_analysis/
cd $RNA_HOME/de/ercc_spikein_analysis/
wget https://genomedata.org/rnaseq-tutorial/ERCC_Controls_Analysis.txt
cat ERCC_Controls_Analysis.txt
```

Using R, load the expected and observed ERCC DE results and produce a visualization. Make sure you're still in the `ercc_spikein_analysis` directory when you start R, because the code below reads and writes files there:

```bash
cd $RNA_HOME/de/ercc_spikein_analysis/
R
```

Work through the following R commands:

```R
library(ggplot2)

#load the ERCC expected fold change values for mix1 vs mix2
ercc_ref = read.table("ERCC_Controls_Analysis.txt", header=TRUE, sep="\t")
names(ercc_ref) = c("id", "ercc_id", "subgroup", "ref_conc_mix_1", "ref_conc_mix_2", "ref_fc_mix1_vs_mix2", "ref_log2_mix1_vs_mix2")
head(ercc_ref)
dim(ercc_ref)

#load the observed fold change values determined by our RNA-seq analysis
rna_de_file = "/workspace/rnaseq/de/ballgown/ref_only/UHR_vs_HBR_gene_results.tsv";
rna_de = read.table(rna_de_file, header=TRUE, sep="\t")
tail(rna_de)
dim(rna_de)

#combine the expected and observed data into a single data table
ercc_ref_de = merge(x = ercc_ref, y = rna_de, by.x = "ercc_id", by.y = "id", all.x = TRUE)
head(ercc_ref_de)
dim(ercc_ref_de)

#convert fold change values to log2 scale
ercc_ref_de$observed_log2_fc = log2(ercc_ref_de$fc)
ercc_ref_de$expected_log2_fc = ercc_ref_de$ref_log2_mix1_vs_mix2

#fit a linear model and calculate R squared between the observed and expected fold change values
model = lm(observed_log2_fc ~ expected_log2_fc, data=ercc_ref_de)
r_squared = summary(model)[["r.squared"]]

#create a scatterplot to compare the observed and expected fold change values
p = ggplot(ercc_ref_de, aes(x = expected_log2_fc, y = observed_log2_fc))
p = p + geom_point(aes(color = subgroup))
p = p + geom_smooth(method = lm)
p = p + annotate("text", 1, 2, label=paste("R^2 =", r_squared, sep=" "))
p = p + xlab("Expected Fold Change (log2 scale)")
p = p + ylab("Observed Fold Change in RNA-seq data (log2 scale)")

#save the plot to a PDF
pdf("ERCC_Ballgown-DE_vs_SpikeInDE.pdf")
print(p)
dev.off()

# Exit the R session
quit(save="no")
```

To view the plot, open this file from your mounted storage (see [Viewing PDFs and downloading files](#viewing-pdfs-and-downloading-files)):

- `rnabio-workspace/rnaseq/de/ercc_spikein_analysis/ERCC_Ballgown-DE_vs_SpikeInDE.pdf`

---

## 3. [Differential expression with edgeR](https://rnabio.org/module-03-expression/0003/03/02/Differential_Expression-edgeR/)

![RNA-seq_Flowchart4](https://rnabio.org/assets/module_3/RNA-seq_Flowchart4-2.png)

In this tutorial you will:

- Make use of the raw counts you generated previously using htseq-count
- edgeR is a Bioconductor package designed specifically for differential expression of count-based RNA-seq data
- This is an alternative to using StringTie/Ballgown to find differentially expressed genes

First, create a directory for results:

```bash
cd $RNA_HOME/
mkdir -p de/htseq_counts/edgeR
cd de/htseq_counts/edgeR
```

Launch R:

```bash
R
```

R code has been provided below. (If you wish to have a script with all of the code, it can be found [here](https://github.com/griffithlab/rnabio.org/blob/master/assets/scripts/Tutorial_edgeR.R), but again, it uses the course's paths, not ours.) Run the R commands below.

```R
# set working directory where output will go
working_dir = "/workspace/rnaseq/de/htseq_counts/edgeR"
setwd(working_dir)

# read in gene mapping
mapping = read.table("/workspace/rnaseq/de/htseq_counts/ENSG_ID2Name.txt", header = FALSE, stringsAsFactors = FALSE, row.names = 1)

# read in count matrix
rawdata = read.table("/workspace/rnaseq/expression/htseq_counts/gene_read_counts_table_all_final.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# Check dimensions
dim(rawdata)

# Require at least 1/6 of samples to have expressed count >= 10
sample_cutoff = (1/6)
count_cutoff = 10

#Define a function to calculate the fraction of values expressed above the count cutoff
getFE = function(data,count_cutoff){
 FE = (sum(data >= count_cutoff) / length(data))
 return(FE)
}

#Apply the function to all genes, and filter out genes not meeting the sample cutoff
fraction_expressed = apply(rawdata, 1, getFE, count_cutoff)
keep = which(fraction_expressed >= sample_cutoff)
rawdata = rawdata[keep, ]

# Check dimensions again to see effect of filtering
dim(rawdata)

#################
# Running edgeR #
#################

# load edgeR
library("edgeR")

# make class labels
class = c(rep("UHR", 3), rep("HBR", 3))

# Get common gene names
Gene = rownames(rawdata)
Symbol = mapping[Gene, 1]
gene_annotations = cbind(Gene, Symbol)

# Make DGEList object
y = DGEList(counts = rawdata, genes = gene_annotations, group = class)
nrow(y)

# TMM Normalization
y = calcNormFactors(y)

# Estimate dispersion
y = estimateCommonDisp(y, verbose = TRUE)
y = estimateTagwiseDisp(y)

# Differential expression test
et = exactTest(y)

# Extract raw counts to add back onto DE results
counts = getCounts(y)

# Print top genes
topTags(et)

# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTests(et, adjust.method = "BH", p = 0.05))

#Get output with BH-adjusted FDR values - all genes, any p-value, unsorted
out = topTags(et, n = "Inf", adjust.method = "BH", sort.by = "none", p.value = 1)$table

#Add raw counts back onto results for convenience (make sure sort and total number of elements allows proper join)
out2 = cbind(out, counts)

#Limit to significantly DE genes
out3 = out2[as.logical(de), ]

# Order by p-value
o = order(et$table$PValue[as.logical(de)], decreasing=FALSE)
out4 = out3[o, ]

# Save table
write.table(out4, file = "DE_sig_genes_edgeR.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

#To exit R type the following
quit(save = "no")
```

---

## 4. [Differential expression with DESeq2](https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/)

In this tutorial you will:

- Make use of the raw counts you generated previously using htseq-count
- DESeq2 is a Bioconductor package designed specifically for differential expression of count-based RNA-seq data
- This is an alternative to using StringTie/Ballgown to find differentially expressed genes

### Setup

Two pieces of information are required to perform analysis with DESeq2. The first is a matrix of raw counts, such as was generated previously while running [HTSeq](https://htseq.readthedocs.io/en/release_0.9.0/). This is important, as DESeq2 normalizes the data, correcting for differences in library size using this type of data. DESeq2 also requires the experimental design, which can be supplied as a data.frame, detailing the samples and conditions.

Start a new R session:

```bash
cd $RNA_HOME/de/htseq_counts/
R
```

All of the packages we need are already installed, so we will load the libraries, set working directories, and read in the raw read counts data:

```R
# define working dir paths
datadir = "/workspace/rnaseq/expression/htseq_counts/"
outdir = "/workspace/rnaseq/de/htseq_counts/deseq2/"

# load R libraries we will use in this section - suppress some harmless warnings from the data.table package
library(DESeq2)
suppressPackageStartupMessages(library(data.table))
library(ggplot2)

# set working directory to data dir
setwd(datadir)

# create our output directory, if it doesn't already exist
if (!dir.exists(outdir)) dir.create(outdir)

# read in the RNAseq read counts for each gene (produced by htseq-count)
htseqCounts = fread("gene_read_counts_table_all_final.tsv")
```

### Format htseq counts data to work with DESeq2

DESeq2 has a number of options for data import and has a function to read HTSeq output files directly. Here the most universal option is used, reading in raw counts from a matrix in simple TSV format (one row per gene, one column per sample). The HTSeq count data that was read in above is stored as an object of class "data.table"; this can be verified with the `class()` function. Before use in this exercise, it is required to convert this object to an appropriate matrix format with gene names as rows and samples as columns.

It should be noted that while the replicate samples are technical replicates (i.e. the same library was sequenced), herein they are treated as biological replicates for illustrative purposes. DESeq2 does have a function to collapse technical replicates, though.

```R
# view the class of the data
class(htseqCounts)

# convert the data.table to matrix format
htseqCounts = as.matrix(htseqCounts)
class(htseqCounts)

# set the gene ID values to be the row names for the matrix
rownames(htseqCounts) = htseqCounts[, "GeneID"]

# now that the gene IDs are the row names, remove the redundant column that contains them
htseqCounts = htseqCounts[, colnames(htseqCounts) != "GeneID"]

# convert the count values from strings (with spaces) to integers, because originally the gene column contained characters, the entire matrix was set to character
class(htseqCounts) = "integer"

# view the first few lines of the gene count matrix
head(htseqCounts)
```

### Filter raw counts

Before running DESeq2 (or any differential expression analysis) it is useful and common practice to pre-filter the data. There are computational benefits to doing this, as the memory size of the objects within R will decrease and DESeq2 will have less data to work through and will be faster. More importantly, removing "low quality" data also reduces the number of statistical tests that will be performed, which in turn reduces the impact of multiple test correction and can lead to more significant genes.

The amount of pre-filtering is up to the analyst. However, it must be done in an unbiased way (i.e. you can NOT remove genes that don't have a difference between conditions of interest, or that show an unexpected pattern).

DESeq2 recommends removing any gene with less than 10 reads across all samples. Below, we filter out a gene unless at least 1 sample has at least 10 reads. Either way, mostly what is being removed here are genes with very little evidence for expression in any sample (in many cases genes with 0 counts in all samples).

```R
# run a filtering step
# i.e. require that for every gene: at least 1 of 6 samples must have counts greater than 10
# get the index of rows that meet this criterion and use that to create a subset of the matrix
# note the dimensions of the matrix before and after filtering with dim()

# breaking apart the command below to understand its outcome
tail(htseqCounts) # look at the raw counts
tail(htseqCounts >= 10) # determine which cells have counts greater than 10
tail(rowSums(htseqCounts >= 10)) # for each gene row, count how many samples have a count greater than 10
tail(rowSums(htseqCounts >= 10) >= 1) # filter to those genes where at least one sample has a count greater than 10
tail(which(rowSums(htseqCounts >= 10) >= 1)) #obtain the index for rows that meet the above filter criteria

dim(htseqCounts)
htseqCounts = htseqCounts[which(rowSums(htseqCounts >= 10) >= 1), ]
dim(htseqCounts)
```

### Specifying the experimental design

As mentioned above, DESeq2 also needs to know the experimental design. In our experimental design this means which samples belong to which condition. The experimental design for the example dataset is quite simple, as there are 6 samples with two conditions to compare (UHR vs HBR). We can create the experimental design right within R. There is one important thing to note: DESeq2 does not check sample names. It expects that the column names in the counts matrix we created correspond exactly to the row names we specify in the experimental design.

```R
# construct a mapping of the meta data for our experiment (comparing UHR cell lines to HBR brain tissues)
# this is defining the biological condition/label for each experimental replicate
# create a simple one column dataframe to start
metaData <- data.frame("Condition" = c("UHR", "UHR", "UHR", "HBR", "HBR", "HBR"))

# convert the "Condition" column to a factor data type
# the arbitrary order of these factors will determine the direction of log2 fold-changes for the genes (i.e. up or down regulated)
# the condition that is listed second here will correspond to a +ve log2 fold change
# genes with counts higher in UHR will have a fold change > 1, or +ve log2 fold change
# genes with counts higher in HBR will have a fold change < 1 (fractions), or -ve log2 fold change

metaData$Condition = factor(metaData$Condition, levels = c("HBR", "UHR"))

# set the row names of the metaData dataframe to be the names of our sample replicates from the read counts matrix
rownames(metaData) = colnames(htseqCounts)

# view the metadata dataframe
head(metaData)

# check that names of htseq count columns match the names of the meta data rows
# use the "all" function which tests whether an entire logical vector is TRUE
all(rownames(metaData) == colnames(htseqCounts))
```

### Construct the DESeq2 object piecing all the data together

With all the data properly formatted, it is now possible to combine all the information required and run the differential expression analysis in one object. This object will hold the input data and intermediary calculations. It is also where the condition to test is specified.

```R
# make deseq2 data sets
# here we are setting up our experiment by supplying: (1) the gene counts matrix, (2) the sample/replicate for each column, and (3) the biological conditions we wish to compare.
# this is a simple example that works for many experiments but these can also get more complex
# for example, including designs with multiple variables such as "~ group + condition",
# and designs with interactions such as "~ genotype + treatment + genotype:treatment".

dds = DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~Condition)

# the design formula above is often a point of confusion, it is useful to put in words what is happening, when we specify "design = ~Condition" we are saying
# regress gene expression on condition, or put another way, model gene expression on condition
# gene expression is the response variable and condition is the explanatory variable
# you can put words to formulas using this [cheat sheet](https://www.econometrics.blog/post/the-r-formula-cheatsheet/)
```

### Running DESeq2

With all the data now in place, DESeq2 can be run. Calling DESeq2 will perform the following actions:

- Estimation of size factors, i.e. accounting for differences in sequencing depth (or library size) across samples.
- Estimation of dispersion, i.e. estimating the biological variability (over and above the expected variation from sampling) in gene expression across biological replicates. This is needed to assess the significance of differences across conditions. Additional work is performed to correct for higher dispersion seen for genes with low expression.
- Perform "independent filtering" to reduce the number of statistical tests performed (see `?results` and [this paper](https://doi.org/10.1073/pnas.0914005107) for details).
- Negative binomial GLM fitting and performing the Wald statistical test.
- Correct p-values for multiple testing using the Benjamini and Hochberg method.

```R
# run the DESeq2 analysis on the "dds" object
dds = DESeq(dds)

# view the first 5 lines of the DE results
res = results(dds)
head(res, 5)
```

### Log fold change shrinkage

It is good practice to shrink the log fold change values. This does exactly what it sounds like, reducing the amount of log fold change for genes where there are few counts, which create huge variability that is not true biological signal. Consider for example a gene for two samples: one sample has 1 read, and one sample has 6 reads. That is a 6-fold change, which is likely not accurate. There are a number of algorithms that can be used to shrink log2 fold change; here we will use the apeglm algorithm, which requires the apeglm package to be installed (it is, in our container).

```R
# shrink the log2 fold change estimates
#   shrinkage of effect size (log fold change estimates) is useful for visualization and ranking of genes

#   In simplistic terms, the goal of calculating "dispersion estimates" and "shrinkage" is also to account for the problem that
#   genes with low mean counts across replicates tend to have higher variability than those with higher mean counts.
#   Shrinkage attempts to correct for this. For a more detailed discussion of shrinkage refer to the DESeq2 vignette

# first get the name of the coefficient (log fold change) to shrink
resultsNames(dds)

# now apply the shrinkage approach
resLFC = lfcShrink(dds, coef = "Condition_UHR_vs_HBR", type = "apeglm")

# make a copy of the shrinkage results to manipulate
deGeneResult = resLFC

# contrast the values before and after shrinkage
# create a temporary, simplified data structure with the values before/after shrinkage, and create a new column called "group" with this label
res_before = res
res_before$group = "before shrinkage"
res_after = deGeneResult
res_after$group = "after shrinkage"
res_combined = rbind(res_before[,c("log2FoldChange","group")], res_after[,c("log2FoldChange","group")])

# adjust the order so that the legend has "before shrinkage" listed first
res_combined$group = factor(res_combined$group, levels = c("before shrinkage", "after shrinkage"))

# look at the fold change values
head(res_before)
head(res_after)

p <- ggplot(res_combined, aes(x = log2FoldChange, color = group)) + geom_density() + theme_bw() +
  scale_color_manual(values = c("before shrinkage" = "tomato4", "after shrinkage" = "slategray")) +
  labs(color = "Shrinkage status")
ggsave(plot=p, filename = paste0(outdir, "before_after_shrinkage.pdf"), device="pdf", width=6, height=4, units="in")
```

To see the plot, open `rnabio-workspace/rnaseq/de/htseq_counts/deseq2/before_after_shrinkage.pdf` from your mounted storage.

How did the results change before and after shrinkage? What direction is each log2 fold change value moving?

Note that for these data, the impact of shrinkage is very subtle, but the pattern is that fold change values move towards 0.

### Annotate gene symbols onto the DE results

DESeq2 was run with Ensembl gene IDs as identifiers, which is not the most human-friendly way to interpret results. Here gene symbols are merged onto the differentially expressed gene list to make the results a bit more interpretable.

```R
# read in gene ID to name mappings (using "fread" an alternative to "read.table")
gene_id_mapping <- fread("/workspace/rnaseq/de/htseq_counts/ENSG_ID2Name.txt", header = FALSE)

# add names to the columns in the "gene_id_mapping" dataframe
setnames(gene_id_mapping, c("ensemblID", "Symbol"))

# view the first few lines of the gene ID to name mappings
head(gene_id_mapping)

# merge on gene names
deGeneResult$ensemblID = rownames(deGeneResult)
deGeneResult = as.data.table(deGeneResult)
deGeneResult = merge(deGeneResult, gene_id_mapping, by = "ensemblID", all.x = TRUE)

# merge the original raw count values onto this final dataframe to aid interpretation
original_counts = as.data.frame(htseqCounts)
original_counts[,"ensemblID"] = rownames(htseqCounts)
deGeneResult = merge(deGeneResult, original_counts, by = 'ensemblID', all.x = TRUE)

# view the final merged dataframe
# based on the raw counts and fold change values what does a negative fold change mean with respect to our two conditions HBR and UHR?
head(deGeneResult)
```

### Data manipulation

With the DE analysis complete, it is useful to view and filter the data frames to only the relevant genes. Here some basic data manipulation is performed, filtering to significant genes at specific thresholds.

```R
# view the top genes according to adjusted p-value
deGeneResult[order(deGeneResult$padj), ]

# view the top genes according to fold change
deGeneResult[order(deGeneResult$log2FoldChange), ]

# determine the number of up/down significant genes at FDR < 0.05 significance level
dim(deGeneResult)[1] # number of genes tested
dim(deGeneResult[deGeneResult$padj < 0.05])[1] #number of significant genes

# order the DE results by adjusted p-value
deGeneResultSorted = deGeneResult[order(deGeneResult$padj), ]

# create a filtered data frame that limits to only the significant DE genes (adjusted p.value < 0.05)
deGeneResultSignificant = deGeneResultSorted[deGeneResultSorted$padj < 0.05]
```

### Save results to files

The data generated will now be written out as tab-separated files (.tsv). Some of the DESeq2 objects are also saved as serialized R (RDS) objects, which can be quickly read back into R later for visualization.

```R
# set the working directory to the output dir where we will store any results files
setwd(outdir)

# save the final DE result (all genes)  to an output file
fwrite(deGeneResultSorted, file="DE_all_genes_DESeq2.tsv", sep="\t")

# save the final DE result (significant genes only)  to an output file
fwrite(deGeneResultSignificant, file="DE_sig_genes_DESeq2.tsv", sep="\t")

# save the DESeq2 objects for the data visualization section
saveRDS(dds, "dds.rds")
saveRDS(res, "res.rds")
saveRDS(resLFC, "resLFC.rds")
```

### Briefly examine the top over-expressed genes

For both conditions (HBR and UHR), let's take a look at the top n genes, but this time according to fold change instead of p-value.

```R
# use the dplyr library to manipulate the dataframe
library(dplyr)

# create a new copy of the data frame, sorted by log2 fold change
deGeneResultSortedFoldchange = arrange(deGeneResultSignificant, log2FoldChange)

# create a convenient data structure with just the top n genes from each condition
top_bottom = bind_rows(
  head(deGeneResultSortedFoldchange, 10) %>% mutate(Set = "Bottom 10"),
  tail(deGeneResultSortedFoldchange, 10) %>% mutate(Set = "Top 10")
)

# visualize data for the top n genes. Simplify the output a bit
print(top_bottom[,c("log2FoldChange","padj","Symbol","UHR_Rep1","UHR_Rep2","UHR_Rep3","HBR_Rep1","HBR_Rep2","HBR_Rep3","Set")])

quit(save="no")
```

### Perform some preliminary exploration of DE genes using web tools

From your mounted storage, open (or copy to your laptop) `rnabio-workspace/rnaseq/de/htseq_counts/deseq2/DE_sig_genes_DESeq2.tsv`. Open this spreadsheet, sort on the "log2FoldChange" column, and find the top 100 significant genes with higher expression in HBR (brain). Also grab `rnabio-workspace/rnaseq/de/htseq_counts/deseq2/DE_all_genes_DESeq2.tsv` (to be used as a list of background genes, or where we want the fold change value for every gene).

- Try querying with the top 100 HBR over-expressed genes using [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost).
- Try querying with the top 100 HBR over-expressed genes using [TissueEnrich](https://tissueenrich.gdcb.iastate.edu/). Use the Tissue Enrichment tool. This tool also requires a Background Gene List; use all genes in `DE_all_genes_DESeq2.tsv` for this purpose. You can also manually explore some individual genes over-expressed in UHR with the Tissue-Specific Genes tool. For example, try *PRAME* and *SERPIND1*, two of the top UHR genes.

**g:Profiler example result**

![g:Profiler example](https://rnabio.org/assets/module_3/Gprofiler-example.png)

**TissueEnrich summary for top 100 HBR over-expressed genes**

![TissueEnrich HBR summary](https://rnabio.org/assets/module_3/TissueEnrich-HBR-Genes-1.png)

**TissueEnrich result for Cerebral Cortex tissue**

![TissueEnrich cerebral cortex](https://rnabio.org/assets/module_3/TissueEnrich-HBR-Genes-2.png)

**TissueEnrich example for UHR over-expressed gene: *PRAME***

![TissueEnrich PRAME](https://rnabio.org/assets/module_3/TissueEnrich-UHR-PRAME.png)

**TissueEnrich example for UHR over-expressed gene: *SERPIND1***

![TissueEnrich SERPIND1](https://rnabio.org/assets/module_3/TissueEnrich-UHR-SERPIND1.png)

Does all of this make sense when we think about the makeup of the HBR and UHR samples? Refer back to the [description of the samples](../week_08/rnaseq_alignment.md#64-rna-seq-data).

---

## 5. [Compare DE genes across methods](https://rnabio.org/module-03-expression/0003/03/04/Differential_Expression-DE-Venn/)

In this section we will compare the DE gene lists obtained from different DE methods (Ballgown, edgeR, and DESeq2).

### Visualize overlap with a Venn diagram

This can be done with simple web tools like:

- [https://www.biovenn.nl/](https://www.biovenn.nl/)
- [https://bioinfogp.cnb.csic.es/tools/venny/](https://bioinfogp.cnb.csic.es/tools/venny/)

Now that you have run the DESeq2 tutorial, compare the significant DE genes to those saved earlier from Ballgown and edgeR:

```bash
head $RNA_HOME/de/ballgown/ref_only/DE_sig_genes_ballgown.tsv
head $RNA_HOME/de/htseq_counts/edgeR/DE_sig_genes_edgeR.tsv
head $RNA_HOME/de/htseq_counts/deseq2/DE_sig_genes_DESeq2.tsv
```

Pull out the gene symbols:

```bash
cd $RNA_HOME/de/

cut -f 1 $RNA_HOME/de/ballgown/ref_only/DE_sig_genes_ballgown.tsv | sort | uniq > ballgown_DE_gene_symbols.txt
cut -f 2 $RNA_HOME/de/htseq_counts/edgeR/DE_sig_genes_edgeR.tsv | sort | uniq | grep -v Symbol > htseq_counts_edgeR_DE_gene_symbols.txt
cut -f 7 $RNA_HOME/de/htseq_counts/deseq2/DE_sig_genes_DESeq2.tsv | sort | uniq | grep -v Symbol > htseq_counts_DESeq2_DE_gene_symbols.txt
```

To get the gene lists you could use `cat` to print out each list in your terminal and then copy/paste:

```bash
cat ballgown_DE_gene_symbols.txt
cat htseq_counts_edgeR_DE_gene_symbols.txt
cat htseq_counts_DESeq2_DE_gene_symbols.txt
```

Alternatively, you could open the lists from your mounted storage. These files should be here:

- `rnabio-workspace/rnaseq/de/ballgown_DE_gene_symbols.txt`
- `rnabio-workspace/rnaseq/de/htseq_counts_edgeR_DE_gene_symbols.txt`
- `rnabio-workspace/rnaseq/de/htseq_counts_DESeq2_DE_gene_symbols.txt`

**Example Venn diagram (two-way comparison: DE genes from StringTie/Ballgown vs HTSeq/DESeq2)**

![Venn: Ballgown vs DESeq2](https://rnabio.org/assets/module_3/venn-ballgown-vs-deseq2.png)

**Example Venn diagram (three-way comparison: DE genes from StringTie/Ballgown vs HTSeq/DESeq2 vs HTSeq/edgeR)**

![Venn: Ballgown vs DESeq2 vs edgeR](https://rnabio.org/assets/module_3/venn-ballgown-vs-deseq2-vs-edger.png)

### Compare the ERCC genes detected as significant by each DE method

Start an R session and produce a visualization that compares how many of the expected DE results for the ERCC spike-in control data are detected by each method:

```bash
R
```

```R
library(ggplot2)

setwd("/workspace/rnaseq/de/ercc_spikein_analysis")

#load the ERCC expected fold change values for mix1 vs mix2
ercc_ref = read.table("ERCC_Controls_Analysis.txt", header=TRUE, sep="\t")[,c("ERCC.ID","expected.fold.change.ratio")]

#load the observed fold change values determined by our RNA-seq analysis
rna_de_fileBG = "/workspace/rnaseq/de/ballgown/ref_only/UHR_vs_HBR_gene_results.tsv"
rna_de_fileDS = "/workspace/rnaseq/de/htseq_counts/deseq2/DE_sig_genes_DESeq2.tsv"
rna_de_fileER = "/workspace/rnaseq/de/htseq_counts/edgeR/DE_sig_genes_edgeR.tsv"

# Read files and rename relevant columns
rna_deBG = read.table(rna_de_fileBG, header=TRUE, sep="\t")[,c("id", "qval")]
colnames(rna_deBG) = c("ensemblID", "Ballgown")
rna_deDS = read.table(rna_de_fileDS, header=TRUE, sep="\t")[,c("ensemblID", "padj")]
colnames(rna_deDS) = c("ensemblID", "DESeq2")
rna_deER = read.table(rna_de_fileER, header=TRUE, sep="\t")[,c("Gene", "FDR")]
colnames(rna_deER) = c("ensemblID", "edgeR")

#merge and melt results from all DE algorithms
ercc_ref_de = merge(x = ercc_ref, y = rna_deBG, by.x = "ERCC.ID", by.y = "ensemblID", all.x = TRUE)
ercc_ref_de = merge(x = ercc_ref_de, y = rna_deDS, by.x = "ERCC.ID", by.y = "ensemblID", all.x = TRUE)
ercc_ref_de = merge(x = ercc_ref_de, y = rna_deER, by.x = "ERCC.ID", by.y = "ensemblID", all.x = TRUE)
ercc_ref_de_long = reshape2::melt(ercc_ref_de, id.vars = c("ERCC.ID","expected.fold.change.ratio"), variable.name = "algorithm")

# transform adjusted p-values into binary
ercc_ref_de_long$value[ercc_ref_de_long$value > 0.05] = NA
ercc_ref_de_long$is_sig = !(is.na(ercc_ref_de_long$value))
# change to factor, or else ggplot will interpret these numbers as positions
ercc_ref_de_long$expected.fold.change.ratio = as.factor(ercc_ref_de_long$expected.fold.change.ratio)

#create a scatterplot to compare the observed and expected fold change values
p = ggplot(ercc_ref_de_long, aes(x = expected.fold.change.ratio, fill =is_sig))
p = p + geom_bar() + facet_grid(~algorithm)
p

ggsave(plot=p, filename = "ERCC_DE_Gene_Method_Comparison.pdf", device="pdf", width=6, height=4, units="in")

quit(save="no")
```

To view the plot, open `rnabio-workspace/rnaseq/de/ercc_spikein_analysis/ERCC_DE_Gene_Method_Comparison.pdf` from your mounted storage.

---

## 6. [DE visualization with Ballgown](https://rnabio.org/module-03-expression/0003/04/01/DE_Visualization_Ballgown/)

![RNA-seq_Flowchart4](https://rnabio.org/assets/module_3/RNA-seq_Flowchart4.png)

Navigate to the correct directory and then launch R:

```bash
cd $RNA_HOME/de/ballgown/ref_only
R
```

Run the R commands below. All results are directed to PDF files in the `ref_only` directory.

First you'll need to load the libraries needed for this analysis and set the working directory, where the output PDFs will be written.

```R
#load libraries
library(ballgown)
library(genefilter)
library(dplyr)
library(devtools)

# set the working directory
setwd("/workspace/rnaseq/de/ballgown/ref_only")
```

Next we'll load our data into R.

```R
# Define the conditions being compared for use later
condition = c("UHR", "UHR", "UHR", "HBR", "HBR", "HBR")

# Load the ballgown object from file
load("bg.rda")

# The load command, loads an R object from a file into memory in our R session.
# You can use ls() to view the names of variables that have been loaded
ls()

# Print a summary of the ballgown object
bg

# Load gene names for lookup later in the tutorial
bg_table = texpr(bg, "all")
bg_gene_names = unique(bg_table[, 9:10])

# Pull the gene and transcript expression data frame from the ballgown object
gene_expression = as.data.frame(gexpr(bg))
transcript_expression = as.data.frame(texpr(bg))

#View expression values for the transcripts of a particular gene symbol of chromosome 22.  e.g. 'TST'
#First determine the transcript_ids in the data.frame that match 'TST', aka. ENSG00000128311, then display only those rows of the data.frame
i = bg_table[, "gene_name"] == "TST"
bg_table[i,]

# Display the transcript ID for a single row of data
ballgown::transcriptNames(bg)[2763]

# Display the gene name for a single row of data
ballgown::geneNames(bg)[2763]

#What if we want to view values for a list of genes of interest all at once?
genes_of_interest = c("TST", "MMP11", "LGALS2", "ISX")
i = bg_table[, "gene_name"] %in% genes_of_interest

bg_table[i,]

# Load the transcript to gene index from the ballgown object
transcript_gene_table = indexes(bg)$t2g
head(transcript_gene_table)

#Each row of data represents a transcript. Many of these transcripts represent the same gene. Determine the numbers of transcripts and unique genes
length(unique(transcript_gene_table[, "t_id"])) #Transcript count
length(unique(transcript_gene_table[, "g_id"])) #Unique Gene count

# Extract FPKM values from the 'bg' object
fpkm = texpr(bg, meas = "FPKM")

# View the last several rows of the FPKM table
tail(fpkm)

# Transform the FPKM values by adding 1 and convert to a log2 scale
fpkm = log2(fpkm + 1)

# View the last several rows of the transformed FPKM table
tail(fpkm)
```

Now we'll start to generate figures with the following R code.

```R
#### Plot #1 - the number of transcripts per gene.
#Many genes will have only 1 transcript, some genes will have several transcripts
#Use the 'table()' command to count the number of times each gene symbol occurs (i.e. the # of transcripts that have each gene symbol)
#Then use the 'hist' command to create a histogram of these counts
#How many genes have 1 transcript?  More than one transcript?  What is the maximum number of transcripts for a single gene?
pdf(file="TranscriptCountDistribution.pdf")
counts=table(transcript_gene_table[, "g_id"])
c_one = length(which(counts == 1))
c_more_than_one = length(which(counts > 1))
c_max = max(counts)
hist(counts, breaks = 50, col = "bisque4", xlab = "Transcripts per gene", main = "Distribution of transcript count per gene")
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text, lty = NULL)
dev.off()

#### Plot #2 - the distribution of transcript sizes as a histogram
#In this analysis we supplied StringTie with transcript models so the lengths will be those of known transcripts
#However, if we had used a de novo transcript discovery mode, this step would give us some idea of how well transcripts were being assembled
#If we had a low coverage library, or other problems, we might get short "transcripts" that are actually only pieces of real transcripts
pdf(file = "TranscriptLengthDistribution.pdf")
hist(bg_table$length, breaks = 50, xlab = "Transcript length (bp)", main = "Distribution of transcript lengths", col = "steelblue")
dev.off()

#### Plot #3 - distribution of gene expression levels for each sample
# Create boxplots to display summary statistics for the FPKM values for each sample
# set color based on condition which is UHR vs. HBR
# set labels perpendicular to axis (las=2)
# set ylab to indicate that values are log2 transformed
pdf(file = "All_samples_FPKM_boxplots.pdf")
boxplot(fpkm, col = as.numeric(as.factor(condition)) + 1,las = 2,ylab = "log2(FPKM + 1)")
dev.off()

#### Plot 4 - BoxPlot comparing the expression of a single gene for all replicates of both conditions
# set border color for each of the boxplots
# set title (main) to gene : transcript
# set x label to Type
# set ylab to indicate that values are log2 transformed

pdf(file = "TST_ENST00000249042_boxplot.pdf")
transcript = which(ballgown::transcriptNames(bg) == "ENST00000249042")[[1]]
boxplot(fpkm[transcript,] ~ condition, border = c(2, 3), main = paste(ballgown::geneNames(bg)[transcript],": ", ballgown::transcriptNames(bg)[transcript]), pch = 19, xlab = "Type", ylab = "log2(FPKM+1)")

# Add the FPKM values for each sample onto the plot
# set plot symbol to solid circle, default is empty circle
points(fpkm[transcript,] ~ jitter(c(2,2,2,1,1,1)), col = c(2,2,2,1,1,1)+1, pch = 16)
dev.off()


#### Plot 5 - Plot of transcript structures observed and expression level for UHR vs HBR with representative replicate

pdf(file = "TST_transcript_structures_expression.pdf")
par(mfrow = c(2, 1))
plotTranscripts(ballgown::geneIDs(bg)[transcript], bg, main = c("TST in HBR"), sample = c("HBR_Rep1"), labelTranscripts = TRUE)
plotTranscripts(ballgown::geneIDs(bg)[transcript], bg, main = c("TST in UHR"), sample = c("UHR_Rep1"), labelTranscripts = TRUE)
dev.off()

# Exit the R session
quit(save = "no")
```

To view the plots, open these files from your mounted storage (see [Viewing PDFs and downloading files](#viewing-pdfs-and-downloading-files)):

- `rnabio-workspace/rnaseq/de/ballgown/ref_only/TranscriptCountDistribution.pdf`
- `rnabio-workspace/rnaseq/de/ballgown/ref_only/TranscriptLengthDistribution.pdf`
- `rnabio-workspace/rnaseq/de/ballgown/ref_only/All_samples_FPKM_boxplots.pdf`
- `rnabio-workspace/rnaseq/de/ballgown/ref_only/TST_ENST00000249042_boxplot.pdf`
- `rnabio-workspace/rnaseq/de/ballgown/ref_only/TST_transcript_structures_expression.pdf`

---

## 7. [DE visualization with DESeq2](https://rnabio.org/module-03-expression/0003/04/02/DE_Visualization-DESeq2/)

![RNA-seq_Flowchart4](https://rnabio.org/assets/module_3/RNA-seq_Flowchart4-2.png)

In this section we will be going over some basic visualizations of the DESeq2 results generated in the "Differential expression with DESeq2" section above. Our goal is to quickly obtain some interpretable results using built-in visualization functions from DESeq2 or recommended packages. For a very extensive overview of DESeq2 and how to visualize and interpret the results, refer to the [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

### Setup

Start a new R session:

```bash
cd $RNA_HOME/de/htseq_counts/deseq2/
R
```

Then load the DESeqDataSet object and the results tables that we saved in the DESeq2 section:

```R
# set the working directory
setwd("/workspace/rnaseq/de/htseq_counts/deseq2/")

# view the contents of this directory
dir()

# load libs
library(DESeq2)
library(data.table)
library(pheatmap)

# Load in the DESeqDataSet object (https://genomedata.org/cri-workshop/deseq2/dds.rds)
dds = readRDS("dds.rds")

# Load in the results object before shrinkage (https://genomedata.org/cri-workshop/deseq2/res.rds)
res = readRDS("res.rds")

# Load in the results object after shrinkage (https://genomedata.org/cri-workshop/deseq2/resLFC.rds)
resLFC = readRDS("resLFC.rds")

# Load in the final results file with all sorted DE results (https://genomedata.org/cri-workshop/deseq2/DE_all_genes_DESeq2.tsv)
deGeneResultSorted = fread("DE_all_genes_DESeq2.tsv")
```

All of the plots in this section are saved as PDFs in `rnabio-workspace/rnaseq/de/htseq_counts/deseq2/` on your mounted storage.

### MA-plot before LFC shrinkage

MA-plots were originally used to evaluate microarray expression data, where M is the log ratio and A is the mean intensity for a gene (both based on scanned intensity measurements from the microarray).

These types of plots are still useful in RNA-seq DE experiments with two conditions, as they can immediately give us information on the number of significantly differentially expressed genes, the ratio of up- vs down-regulated genes, and any outliers. To interpret these plots it is important to keep a couple of things in mind. The Y axis (M) is the log2 fold change between the two conditions tested; a higher fold change indicates a greater difference between condition A and condition B. The X axis (A) is a measure of read alignment to a gene, so as you go higher on the X axis you are looking at regions which have higher totals of aligned reads. In other words, the gene is "more" expressed overall (with the caveat that gene length is not being taken into account by raw read counts here).

Using the built-in `plotMA` function from DESeq2, we also see that the genes are color-coded by a significance threshold (e.g., adjusted p-value < 0.1). Genes with higher expression values and higher fold changes are more often significant, as one would expect.

```R
# use DESeq2 built in MA-plot function
pdf("maplot_preShrink.pdf")
plotMA(res, alpha = 0.1, ylim=c(-2, 2), cex=1.5, main = "MA-plot before LFC shrinkage")
dev.off()
```

### MA-plot after LFC shrinkage

When we ran DESeq2 we obtained two results, one with and one without log fold change shrinkage. When you have genes with low hits you can get some very large fold changes. For example, 1 hit on a gene vs 6 hits on a gene is a 6x fold change. This high level of variance, though, is probably quantifying noise instead of real biology. Running `plotMA` on the results where we applied an algorithm for log fold change shrinkage, we can see that this "effect" is somewhat controlled for. Note that while shrinking LFC is part of a typical DE workflow, there are cases where you would not want to perform this, namely when there is already low variation amongst samples (i.e. from technical replicates), as most shrinkage algorithms rely on some variability to build a prior distribution.

```R
# use DESeq2 built in MA-plot function
pdf("maplot_postShrink.pdf")
plotMA(resLFC, alpha = 0.1, ylim = c(-2, 2), cex=1.5, main = "MA-plot after LFC shrinkage")
dev.off()
```

The effect is very subtle here due to the focused nature of our dataset (chr22 genes only), but if you toggle between the two plots and look in the upper left and bottom left corners you can see some fold change values are moving closer to 0.

### Viewing individual gene counts between two conditions

Often it may be useful to view the normalized counts for a gene amongst our samples. DESeq2 provides a built-in function for that which works off of the dds object. Here we view SEPT3, which we can see in our DE output is significantly higher in the HBR cohort, and PRAME, which is significantly higher in the UHR cohort. This is useful as we can see the per-sample distribution of our corrected counts; we can immediately determine if there are any outliers within each group and investigate further if need be.

```R
# hint! you defined intgroup when creating the dds object, you can view the name by printing the dds objct in your R session
# dds

pdf("normalized_count_examples.pdf")

# view SEPT3 normalized counts
plotCounts(dds, gene = "ENSG00000100167", intgroup = "Condition", main = "SEPT3")

# view PRAME normalized counts
plotCounts(dds, gene = "ENSG00000185686", intgroup = "Condition", main = "PRAME")

dev.off()
```

### Viewing pairwise sample clustering

It may often be useful to view inter-sample relatedness. In other words, how similar or dissimilar samples are to one another overall. While not part of the DESeq2 package, there is a convenient library that can easily construct a hierarchically clustered heatmap from our DESeq2 data. It should be noted that using "raw count" data when doing a distance calculation is not ideal; the count data should be transformed using `vst()` or `rlog()`, which can be performed directly on the dds object. The reason for this is described in detail in the DESeq2 manuscript. Suffice it to say that transforming gene variance to be more homoskedastic will make inferences of sample relatedness more interpretable.

```R
# note that we use rlog because we don't have a large number of genes, for a typical DE experiment with 1000's of genes use the vst() function
rld <- rlog(dds, blind = FALSE)

# view the structure of this object
rld

# compute sample distances (the dist function uses the euclidean distance metric by default)
# in this command we will pull the rlog transformed data ("regularized" log2 transformed, see ?rlog for details) using "assay"
# then we transpose that data using t()
# then we calculate distance values using dist()
# the distance is calculated for each vector of sample gene values, in a pairwise fashion comparing all samples

# view the first few lines of raw data
head(assay(dds))

# see the rlog transformed data
head(assay(rld))

# see the impact of transposing the matrix
t(assay(rld))[1:6, 1:5]

# see the distance values
dist(t(assay(rld)))

# put it all together and store the result
sampleDists <- dist(t(assay(rld)))

# convert the distance result to a matrix
sampleDistMatrix <- as.matrix(sampleDists)

# view the distance numbers directly in the pairwise distance matrix
head(sampleDistMatrix)

pdf("distance_sample_heatmap.pdf", width = 8, height = 8)

# construct clustered heatmap, important to use the computed sample distances for clustering
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists)

dev.off()
```

Instead of a distance metric we could also use a similarity metric, such as a Pearson correlation.

There are many correlation and distance options:

- Similarity metrics: "pearson", "kendall", "spearman"
- Distance metrics: "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"

```R
sampleCorrs = cor(assay(rld), method = "pearson")
sampleCorrMatrix = as.matrix(sampleCorrs)
head(sampleCorrMatrix)

pdf("similarity_sample_heatmap.pdf", width = 8, height = 8)

pheatmap(sampleCorrMatrix)

dev.off()
```

Instead of boiling all the gene count data for each sample down to a distance metric, you can get a similar sense of the pattern by just visualizing all the genes and their expression at once:

```R
pdf("all_gene_heatmap.pdf", width = 10, height = 10)

# because there are so many gene we choose not to display them

pheatmap(mat = t(assay(rld)), show_colnames = FALSE)

dev.off()
```

Make a heatmap for just significant genes (e.g., FC >= 3 and padj < 0.001). Note, we chose these cutoffs arbitrarily to result in a number of genes that can be legibly labeled.

```R
#First, get a list of just significant genes
deGeneSigResult=deGeneResultSorted[(deGeneResultSorted$padj<0.001 & (deGeneResultSorted$log2FoldChange>=3 | deGeneResultSorted$log2FoldChange<=-3)),]
SigGenes=deGeneSigResult$ensemblID
SigGeneSymbols=deGeneSigResult$Symbol

#extract expression data
expression_data=assay(rld)

#Limit to significant genes
expression_data=expression_data[SigGenes,]

#Create a heatmap as before
pdf("sig_gene_heatmap.pdf", width = 10, height = 10)
pheatmap(mat = t(expression_data), show_colnames = TRUE, labels_col=SigGeneSymbols)
dev.off()

quit(save="no")
```

---

## 8. [DE visualization advanced](https://rnabio.org/module-03-expression/0003/04/03/DE_Visualization_AdvancedR/)

### Supplementary R DE visualization

Occasionally you may wish to reformat and work with expression estimates in R in an ad hoc way. Here, we provide an optional/advanced tutorial on how to visualize your results in R and perform "old school" (non-Ballgown, non-DESeq2) visualization of your data.

In this tutorial you will:

- Learn basic R usage and commands (common plots, and data manipulation tasks)
- Examine the expression estimates
- Create an MDS plot to visualize the differences between/among replicates, library prep methods and UHR versus HBR
- Examine the differential expression estimates
- Visualize the expression estimates and highlight those genes that appear to be differentially expressed
- Ask how reproducible technical replicates are

Expression and differential expression files will be read into R. The R analysis will make use of the gene-level expression estimates from HISAT2/StringTie (TPM values) and differential expression results from HISAT2/htseq-count/DESeq2 (fold changes and p-values).

Navigate to the correct directory and then launch R:

```bash
cd $RNA_HOME
R
```

First you'll load your libraries and your data. All of the plots in this section are written to a new `visualization_advanced` directory.

```R
#Load your libraries
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(ggrepel)

#Set a base working directory
setwd("/workspace/rnaseq/")


#Import expression results (TPM values) from the HISAT2/Stringtie pipeline (https://genomedata.org/cri-workshop/gene_tpm_all_samples.tsv)
gene_expression = read.table("/workspace/rnaseq/expression/stringtie/ref_only/gene_tpm_all_samples.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

#Import gene name mapping file (https://genomedata.org/cri-workshop/ENSG_ID2Name.txt)
gene_names=read.table("/workspace/rnaseq/de/htseq_counts/ENSG_ID2Name.txt", header = FALSE, stringsAsFactors = FALSE)

colnames(gene_names) = c("gene_id", "gene_name")

#Import DE results from the HISAT2/htseq-count/DESeq2 pipeline (https://genomedata.org/cri-workshop/deseq2/DE_all_genes_DESeq2.tsv)
results_genes = read.table("/workspace/rnaseq/de/htseq_counts/deseq2/DE_all_genes_DESeq2.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#Set a directory for the output to go to
output_dir = "/workspace/rnaseq/visualization_advanced/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
setwd(output_dir)
```

Let's briefly explore the imported data:

```R
#### Working with 'dataframes'
#View the first five rows of data (all columns) in the gene_expression (Stringtie TPM) dataframe
head(gene_expression)

#View the column names
colnames(gene_expression)

#View the row names
row.names(gene_expression)

#Determine the dimensions of the dataframe.  'dim()' will return the number of rows and columns
dim(gene_expression)

#Get the first 3 rows of data and a selection of columns
gene_expression[1:3, c(1:3, 6)]

#Do the same thing, but using the column names instead of numbers
gene_expression[1:3, c("HBR_Rep1", "HBR_Rep2", "HBR_Rep3", "UHR_Rep3")]

#Now, explore the differential expression (DESeq2 results)
head(results_genes)
dim(results_genes)

#Assign some colors for use later.  You can specify color by RGB, Hex code, or name
#To get a list of color names:
colours()
data_colors = c("tomato1", "tomato2", "tomato3", "royalblue1", "royalblue2", "royalblue3")
```

The following code blocks generate various plots using the above data set.

```R
#### Plot #1 - View the range of values and general distribution of TPM values for all 6 libraries
#Create boxplots for this purpose
#Display on a log2 scale and set a minimum non-zero value to avoid log2(0)
min_nonzero = 1

# Set the columns for finding TPM and create shorter names for figures
data_columns = c(1:6)
short_names = c("HBR_1", "HBR_2", "HBR_3", "UHR_1", "UHR_2", "UHR_3")

pdf(file = "All_samples_TPM_boxplots.pdf")
boxplot(log2(gene_expression[, data_columns] + min_nonzero), col = data_colors, names = short_names, las = 2, ylab = "log2(TPM)", main = "Distribution of TPMs for all 6 libraries")
#Note that the bold horizontal line on each boxplot is the median
dev.off()

#### Plot #2 - plot a pair of replicates to assess reproducibility of technical replicates
#Transform the data by converting to log2 scale after adding an arbitrary small value to avoid log2(0)
x = gene_expression[, "UHR_Rep1"]
y = gene_expression[, "UHR_Rep2"]
pdf(file = "UHR_Rep1_vs_Rep2_scatter.pdf")
plot(x = log2(x + min_nonzero), y = log2(y + min_nonzero), pch = 16, col = "blue", cex = 0.25, xlab = "TPM (UHR, Replicate 1)", ylab = "TPM (UHR, Replicate 2)", main = "Comparison of expression values for a pair of replicates")

#Add a straight line of slope 1, and intercept 0
abline(a = 0, b = 1)

#Calculate the correlation coefficient and display in a legend
rs = cor(x, y)^2
legend("topleft", paste("R squared = ", round(rs, digits = 3), sep = ""), lwd = 1, col = "black")
dev.off()

#### Plot #3 - Scatter plots with a large number of data points can be misleading ... regenerate this figure as a density scatter plot
pdf(file = "UHR_Rep1_vs_Rep2_SmoothScatter.pdf")
colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
smoothScatter(x = log2(x + min_nonzero), y = log2(y + min_nonzero), xlab = "TPM (UHR, Replicate 1)", ylab = "TPM (UHR, Replicate 2)", main = "Comparison of expression values for a pair of replicates", colramp = colors, nbin = 200)
dev.off()

#### Plot #4 - Scatter plots of all sets of replicates on a single plot
#Create a function that generates an R plot.  This function will take as input the two libraries to be compared and a plot name
plotCor = function(lib1, lib2, name){
	x = gene_expression[, lib1]
	y = gene_expression[, lib2]
	colors = colorRampPalette(c("white", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
	smoothScatter(x = log2(x + min_nonzero), y = log2(y + min_nonzero), xlab = lib1, ylab = lib2, main = name, colramp = colors, nbin = 275)
	abline(a = 0, b = 1)
        zero_count = length(which(x == 0)) + length(which(y == 0))
	rs = cor(x, y, method = "pearson")^2
	legend_text = c(paste("R squared = ", round(rs, digits = 3), sep=""), paste("Zero count = ", zero_count, sep = ""))
	legend("topleft", legend_text, lwd = c(1, NA), col = "black", bg = "white", cex = 0.8)
}

#Now make a call to our custom function created above, once for each library comparison
pdf(file = "UHR_All_Reps_SmoothScatter.pdf")
par(mfrow = c(1, 3))
plotCor("UHR_Rep1", "UHR_Rep2", "UHR_1 vs UHR_2")
plotCor("UHR_Rep2", "UHR_Rep3", "UHR_2 vs UHR_3")
plotCor("UHR_Rep1", "UHR_Rep3", "UHR_1 vs UHR_3")

dev.off()

#### Compare the correlation between all replicates
#Do we see the expected pattern for all six libraries (i.e. replicates most similar, then UHR vs. HBR)?

#Calculate the TPM sum for all 6 libraries
gene_expression[,"sum"] = apply(gene_expression[,data_columns], 1, sum)

#Identify the genes with a grand sum TPM of at least 5 - we will filter out the genes with very low expression across the board
i = which(gene_expression[,"sum"] > 5)

#Calculate the correlation between all pairs of data
r = cor(gene_expression[i,data_columns], use = "pairwise.complete.obs", method = "pearson")

#Print out these correlation values
r


#### Plot #5 - Convert correlation to 'distance', and use 'multi-dimensional scaling' to display the relative differences between libraries
#This step calculates 2-dimensional coordinates to plot points for each library
#Libraries with similar expression patterns (highly correlated to each other) should group together
#What pattern do we expect to see, given the types of libraries we have (technical replicates, biological replicates, tumor/normal)?
pdf(file = "UHR_vs_HBR_MDS.pdf")
d = 1 - r
mds = cmdscale(d, k = 2, eig = TRUE)
par(mfrow = c(1,1))
plot(mds$points, type = "n", xlab = "", ylab = "", main = "MDS distance plot (all non-zero genes)", xlim = c(-0.12, 0.12), ylim = c(-0.12, 0.12))
points(mds$points[, 1], mds$points[, 2], col = "grey", cex = 2, pch = 16)
text(mds$points[, 1], mds$points[, 2], short_names, col = data_colors)
dev.off()


#### Plot #6 - View the distribution of differential expression values as a histogram
#Display only those results that are significant according to DESeq2 (loaded above)
pdf(file = "UHR_vs_HBR_DE_FC_distribution.pdf")
sig = which(results_genes$pvalue < 0.05)
hist(results_genes[sig, "log2FoldChange"], breaks = 50, col = "seagreen", xlab = "log2(Fold change) UHR vs HBR", main = "Distribution of differential expression values")
abline(v = -2, col = "black", lwd = 2, lty = 2)
abline(v = 2, col = "black", lwd = 2, lty = 2)
legend("topleft", "Fold-change > 2", lwd = 2, lty = 2)
dev.off()

#### Plot #7 - Display the mean expression values from UHR and HBR and mark those that are significantly differentially expressed
pdf(file="UHR_vs_HBR_mean_TPM_scatter.pdf")

gene_expression[, "HBR_mean"] = apply(gene_expression[,c(1:3)], 1, mean)
gene_expression[, "UHR_mean"] = apply(gene_expression[,c(4:6)], 1, mean)

x = log2(gene_expression[, "UHR_mean"] + min_nonzero)
y = log2(gene_expression[, "HBR_mean"] + min_nonzero)
#Name each value by its Ensembl gene ID, so we can look genes up by ID below
names(x) = rownames(gene_expression)
names(y) = rownames(gene_expression)
plot(x = x, y = y, pch = 16, cex = 0.25, xlab = "UHR TPM (log2)", ylab = "HBR TPM (log2)", main = "UHR vs HBR TPMs")
abline(a = 0, b = 1)

#The TPM values (StringTie) and the DE results (DESeq2) are two different tables, with different genes in a different order,
#so we match them up by Ensembl gene ID rather than by row number
sig_ids = intersect(results_genes[sig, "ensemblID"], names(x))
points(x = x[sig_ids], y = y[sig_ids], col = "magenta", pch = 16, cex = 0.5)
legend("topleft", "Significant", col = "magenta", pch = 16)

#Get the gene symbols for the top N (according to corrected p-value) and display them on the plot
sig_genes = results_genes[sig, ]
sig_genes = sig_genes[sig_genes$ensemblID %in% names(x), ]
top_genes = head(sig_genes[order(sig_genes$padj), ], 25)
text(x[top_genes$ensemblID], y[top_genes$ensemblID], top_genes$Symbol, col = "black", cex = 0.75, srt = 45)

dev.off()

#### Plot #8 - Create a heatmap to visualize expression differences between the six samples
#Define custom dist and hclust functions for use with heatmaps
mydist = function(c) {dist(c, method = "euclidian")}
myclust = function(c) {hclust(c, method = "average")}

#Create a subset of significant genes with p-value<0.05 and log2 fold-change >= 2
sigpi = which(results_genes[, "pvalue"] < 0.05)
sigp = results_genes[sigpi, ]
sigfc = which(abs(sigp[, "log2FoldChange"]) >= 2)
sigDE = sigp[sigfc, ]

pdf(file = "UHR_vs_HBR_heatmap.pdf")
main_title = "sig DE Genes"
par(cex.main = 0.8)
sigDE_genes = sigDE[, "ensemblID"]
sigDE_genenames = sigDE[, "Symbol"]

data = log2(as.matrix(gene_expression[as.vector(sigDE_genes), data_columns]) + 1)
heatmap.2(data, hclustfun = myclust, distfun = mydist, na.rm = TRUE, scale = "none", dendrogram = "both", margins = c(10,4), Rowv = TRUE, Colv = TRUE, symbreaks = FALSE, key = TRUE, symkey = FALSE, density.info = "none", trace = "none", main = main_title, cexRow = 0.3, cexCol = 1, labRow = sigDE_genenames, col = rev(heat.colors(75)))
dev.off()

#### Plot #9 - Volcano plot

# Set differential expression status for each gene - default all genes to "no change"
results_genes$diffexpressed = "No"

# if log2Foldchange > 2 and padj < 0.05, set as "Up regulated"
results_genes$diffexpressed[results_genes$log2FoldChange >= 2 & results_genes$padj < 0.05] = "Higher in UHR"

# if log2Foldchange < -2 and padj < 0.05, set as "Down regulated"
results_genes$diffexpressed[results_genes$log2FoldChange <= -2 & results_genes$padj < 0.05] = "Higher in HBR"

# write the gene names of those significantly upregulated/downregulated to a new column
results_genes$gene_label = NA
results_genes$gene_label[results_genes$diffexpressed != "No"] = results_genes$Symbol[results_genes$diffexpressed != "No"]

pdf(file = "UHR_vs_HBR_volcano.pdf")

ggplot(data = results_genes[results_genes$diffexpressed != "No",], aes(x = log2FoldChange, y = -log10(padj), label = gene_label, color = diffexpressed)) +
  xlab("log2Foldchange") +
  scale_color_manual(name = "Differentially expressed", values=c("blue", "red")) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  geom_vline(xintercept = c(-2.0, 2.0), col = "black", linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dotted") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  geom_point(data = results_genes[results_genes$diffexpressed == "No",], aes(x = log2FoldChange, y = -log10(padj)), colour = "black")

dev.off()

#To exit R type:
quit(save = "no")
```

To view the plots, open these files in `rnabio-workspace/rnaseq/visualization_advanced/` from your mounted storage:

- `All_samples_TPM_boxplots.pdf`
- `UHR_Rep1_vs_Rep2_scatter.pdf`
- `UHR_Rep1_vs_Rep2_SmoothScatter.pdf`
- `UHR_All_Reps_SmoothScatter.pdf`
- `UHR_vs_HBR_MDS.pdf`
- `UHR_vs_HBR_DE_FC_distribution.pdf`
- `UHR_vs_HBR_mean_TPM_scatter.pdf`
- `UHR_vs_HBR_heatmap.pdf`
- `UHR_vs_HBR_volcano.pdf`

### Visual comparison of example genes from the volcano plot

One can manually explore interesting-looking genes from the volcano plot. In this case our analysis involves comparison of RNA isolated from tissues of different types (HBR -> brain tissue, UHR -> a collection of cancer cell lines). So, in this analysis it might make sense to explore candidates in a tissue expression atlas such as [GTEx](https://www.gtexportal.org/).

Looking at our gene plot, two example genes we could look at are SEPT3 (significantly higher in HBR) and PRAME (significantly higher in UHR).

*Expression of SEPT3 across tissues according to GTEx*

![SEPT3](https://rnabio.org/assets/module_3/SEPT3.png)

Note that this gene appears to be highly expressed in brain tissues.

*Expression of PRAME across tissues according to GTEx*

![PRAME](https://rnabio.org/assets/module_3/PRAME.png)

Note that this gene appears to be almost uniquely expressed in testis tissue. Since one of the cell lines in the UHR sample is a testicular cancer cell line, this makes sense.

### Practical exercise 10 (advanced)

Assignment: Use R to create a volcano plot for the differentially expressed genes you identified with Ballgown in Practical Exercise 9.

- Hint: Follow the example R code above.
- Hint: You could import the Ballgown data object (e.g., `bg.rda`) that you should have saved in Practical Exercise 9 as a source of DE results.

Solution: When you are ready you can check your approach against the [Solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-10---volcano-plot).

---

## 9. For-credit deliverable

For-credit students: please provide a screenshot of a volcano plot representing the differentially expressed genes in the dataset. Set the p-value cutoff to 0.01 and the log fold change threshold to 2. Send it to John, as usual.

---

## 10. Cleaning up

When you're done for the day, type `quit(save = "no")` to leave R (if you're still in it), then `exit` to leave the container and release your interactive job.

**Don't delete anything in `rnabio-workspace`!** You may want these results later in the course.

---

## 11. Cheat sheet

```bash
# --- on the login node: start the RNA-seq toolkit container ---
srun -A compute2-workshop -p workshop -c 2 --mem=8G --pty \
  --container-image=griffithlab/rnaseq-toolkit:1.1.1 \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch/$USER/rnabio-workspace:/workspace" \
  --container-workdir=/workspace \
  /bin/bash

# --- inside the container: load the course variables (every session!) ---
source /workspace/rnabio_env.sh
echo $RNA_HOME                                  # should print /workspace/rnaseq

# --- where this week's results live ---
$RNA_HOME/de/ballgown/ref_only/                 # Ballgown DE results and plots
$RNA_HOME/de/htseq_counts/edgeR/                # edgeR DE results
$RNA_HOME/de/htseq_counts/deseq2/               # DESeq2 DE results, .rds objects, and plots
$RNA_HOME/de/ercc_spikein_analysis/             # ERCC DE comparisons
$RNA_HOME/visualization_advanced/               # "advanced" plots, including the volcano plot

# --- R ---
R                                               # start R (prompt becomes >)
quit(save = "no")                               # leave R
# in R, always use /workspace/rnaseq/... paths, never ~/workspace/...

# --- view PDFs: mount storage on your laptop (Mac) ---
open smb://storage1.ris.wustl.edu/c.a.miller/Active/bfx-workshop-scratch/<washukey>
```


----

This page is adapted from the [RNA-seq Bioinformatics course](https://rnabio.org/course/) by the Griffith Lab (Malachi Griffith, Obi Griffith, and colleagues), released under the [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) license. The commands have been adjusted to run in our container on Compute2.