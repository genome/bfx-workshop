# RNA-seq Part 1: Inputs and Alignment on the Cluster

This week we start working through parts of the [RNA-seq Bioinformatics course](https://rnabio.org/course/) (rnabio.org), developed by the Griffith Lab and used here with permission. That course was written to run on an Amazon (AWS) cloud machine. We'll run it on Compute2 instead, inside a Docker container that has all of the course's tools pre-installed.

By the end of this exercise you should be able to:

- Launch an interactive Slurm job inside the RNA-seq toolkit container, with a disk mounted to hold your analyses
- Set up bash environment variables to store info for reuse
- Download and explore a reference genome, annotation, and RNA-seq data
- Build a HISAT2 index, then align RNA-seq reads with HISAT2
- Run pre- and post-alignment QC with FastQC, MultiQC, Picard, and RSeQC

> [!NOTE]
> As always, we strongly suggest typing the commands out yourself rather than copying and pasting. The exceptions are the very long commands from the course website (like the `hisat2` commands), where a typo can cost you a lot of time. Copy those, but read them carefully first so you know what each piece does.

The first part of this page sets up your workspace and container. The rest of it contains Modules 1 and 2 of the course, with all of the commands already adjusted to work on the cluster. **Follow the commands on this page, rather than the ones on rnabio.org.** Each section heading links to the original lesson if you'd like to see it.

## Contents

**Setup**

1. [The container](#1-the-container)
2. [Set up your workspace](#2-set-up-your-workspace)
3. [Launch an interactive job in the container](#3-launch-an-interactive-job-in-the-container)
4. [Set up the course environment variables](#4-set-up-the-course-environment-variables)
5. [A few notes before you start](#5-a-few-notes-before-you-start)

**Module 1: Inputs** 

- 6.1 [Reference genomes](#61-reference-genomes)
- 6.2 [Annotations](#62-annotations)
- 6.3 [Indexing](#63-indexing)
- 6.4 [RNA-seq data](#64-rna-seq-data)
- 6.5 [Pre-alignment QC](#65-pre-alignment-qc)

**Module 2: Alignment** 

- 7.1 [Adapter trimming (optional)](#71-adapter-trimming-optional)
- 7.2 [Alignment](#72-alignment)
- 7.3 [Indexing BAM files](#73-indexing-bam-files)
- 7.4 [Allele expression](#74-allele-expression)
- 7.5 [Alignment QC](#75-alignment-qc)

**Wrapping up**

8. [Optional challenge: alignment on a new dataset](#8-optional-challenge-alignment-on-a-new-dataset)
9. [Cleaning up](#9-cleaning-up)
10. [Cheat sheet](#10-cheat-sheet)

---

## 1. The container

All of the tools for this course (samtools, HISAT2, StringTie, FastQC, MultiQC, Picard, fastp, RSeQC, R, and many more) have been bundled into a single docker image:

```
griffithlab/rnaseq-toolkit:1.1.1
```

Two things are different from the containers we've used so far:

- **This image is big** (over 2 GB). The first time you launch it, it may take 5–10 minutes to download and unpack. Be patient, and don't cancel it!
- **The course expects all of its work to live under a directory called `/workspace`.** In week 3 we mounted our storage directory at the *same* path inside the container. This time, we'll mount a directory from our storage allocation at a *different* path inside the container, `/workspace`, so that the course's paths line up.

---

## 2. Set up your workspace

Log in to Compute2:

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu
```
(as always, replace \<washukey\> with your username)

Then create a directory for the course inside your workshop folder. We'll be using this same workspace for the rest of the RNA-seq weeks, so we put it at the top level of `~/workshop`, rather than inside a `week08` folder:

```bash
cd ~/workshop
mkdir -p rnabio-workspace/rnaseq
```

Remember that `~/workshop` is a softlink (shortcut) that points to your directory on the storage allocation. If you don't have this link set up, revisit that section of the [week 02 assignment](../week_02/compute2-slurm-intro.md#3-where-is-my-data).

The container won't know about that shortcut, so when we mount the directory we need to give Slurm the *real* path. Let's see what it is:

```bash
readlink -f ~/workshop/rnabio-workspace
```

It should print something like:

```
/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/<washukey>/rnabio-workspace
```

### Question 1

What would happen if we created our workspace in our home directory (`~/rnabio-workspace`) instead of under `~/workshop`?

<details>
<summary>Solution</summary>

It would work for a little while, but home directories on the cluster are very small, and you'd run out of space partway through the course (the alignments alone take up a few GB). Data belongs on a storage allocation, not in your home directory.

</details>

---

## 3. Launch an interactive job in the container

Now we'll request an interactive job running inside the RNA-seq toolkit image, with our workspace mounted at `/workspace`:

```bash
srun -A compute2-workshop -p workshop -c 2 --mem=8G --pty \
  --container-image=griffithlab/rnaseq-toolkit:1.1.1 \
  --container-mounts="/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/$USER/rnabio-workspace:/workspace" \
  --container-workdir=/workspace \
  /bin/bash
```

Most of this should look familiar from weeks 3 and 4:

| Option | Meaning |
|---|---|
| `-A compute2-workshop` | charge the job to the workshop account |
| `-p workshop` | submit to the workshop's interactive partition |
| `-c 2 -N 1 --mem=8G` | two CPUs and 8 GB of RAM |
| `--pty ... /bin/bash` | give me an interactive shell |
| `--container-image=...` | run that shell inside the RNA-seq toolkit image |
| `--container-mounts="<source>:<destination>"` | make `<source>` on the cluster visible at `<destination>` inside the container. This time the two sides are **different**: our `rnabio-workspace` folder shows up inside the container as `/workspace` |
| `--container-workdir=/workspace` | **new!** start the shell in `/workspace`, instead of dropping us at `/` like we saw in week 3 |

`$USER` is an environment variable that holds your username, so that path fills itself in with your washukey. (Try `echo $USER` on the login node to see it.)

> [!NOTE]
> Remember, the first launch of this image can take 5–10 minutes. You won't see much output while it's downloading. When you get a prompt back, you're in!

Once you have a prompt, let's verify a few things:

1) that we're on a compute node, not the head node

```
hostname
```
> c2-node-004

2) that one of the tools we expect to be in the docker image are actually there:

```
hisat2 --version | head -n 1
```

3) That our disk mounting worked correctly:

```bash
ls /workspace
```
> rnaseq

----

> [!IMPORTANT]
> Unless otherwise noted, **all of the commands for the rest of this assignment are run inside this container**. If you get disconnected or `exit`, re-run the `srun` command above and then `source` your environment file (see the next section). Everything you've created will still be there, because it was written to your storage directory, not into the container.

---

## 4. Set up the course environment variables

The course uses a set of environment variables (like `$RNA_HOME` and `$RNA_REFS_DIR`) as shortcuts for paths. Almost every command on rnabio.org refers to them, so if they're not set, commands will fail or put files in strange places.

We'll write the variables into a small file in our workspace, and load it each time we start the container. Go ahead and copy/paste this block inside the container to create this file full of variables:

```bash
cat > /workspace/rnabio_env.sh <<'EOF'
# RNAbio course environment variables, adapted for Compute2
export RNA_HOME=/workspace/rnaseq
export RNA_DATA_DIR=$RNA_HOME/data
export RNA_DATA_TRIM_DIR=$RNA_DATA_DIR/trimmed
export RNA_REFS_DIR=$RNA_HOME/refs
export RNA_REF_INDEX=$RNA_REFS_DIR/chr22_with_ERCC92
export RNA_REF_FASTA=$RNA_REF_INDEX.fa
export RNA_REF_GTF=$RNA_REF_INDEX.gtf
export RNA_ALIGN_DIR=$RNA_HOME/alignments/hisat2
export PICARD=/home/ubuntu/bin/picard.jar
EOF
```

(This is a "heredoc": everything between `<<'EOF'` and the final `EOF` is written into the file)

Now load the variables and check that they're set:

```bash
source /workspace/rnabio_env.sh
env | grep RNA
```

**You'll need to re-run `source /workspace/rnabio_env.sh` every time you start a new container session.** If a command from the course ever complains about a missing directory, or something like `cd: too many arguments`, the first thing to check is `echo $RNA_HOME`.

### Question 2

Why do we put the `export` commands in a file and `source` it, rather than just typing them in once?

<details>
<summary>Solution</summary>

Environment variables only last as long as the shell they were set in. When you `exit` the container (or get disconnected), that shell is gone and the variables go with it. Keeping them in a file means one short command restores all of them, exactly the same way, every time.

</details>

> [!TIP]
> **Optional shortcut:** if you'd rather not type the `source` command every time, add this line to the end of your `~/.bashrc` (using `nano ~/.bashrc`):
>
> ```bash
> [ -f /workspace/rnabio_env.sh ] && source /workspace/rnabio_env.sh
> ```
>
> It says "if `/workspace/rnabio_env.sh` exists, load it." That file only exists inside the container, so this line does nothing on the login node, and loads the course variables automatically whenever you start the toolkit container.

---

## 5. A few notes before you start

### Viewing HTML reports

Several of the steps below produce HTML reports (from FastQC, MultiQC, and fastp) that are meant to be opened in a web browser. There's no browser on the cluster, so we'll look at them from our laptops.

The easiest way is to mount your storage on your laptop, the same way we did in [week 4](../week_04/alignment_assignment.md#5-viewing-the-data-with-igv):

- **Mac:** in a Terminal on your laptop, run `open smb://storage1.ris.wustl.edu/c.a.miller/Active/bfx-workshop-scratch/<washukey>`
- **Windows:** in File Explorer, go to `\\storage1.ris.wustl.edu\c.a.miller\Active\bfx-workshop-scratch\<washukey>` (and log in as `accounts\<washukey>`)

Your course files will be in the `rnabio-workspace/rnaseq/` folder, and you can double-click any HTML file to open it in your browser. (Alternatively, you can `scp` the files to your laptop as we did in [week 3](../week_03/sequence-data-and-containers.md#8-looking-at-the-results).)

### CPUs

We asked Slurm for 2 CPUs, so whenever a tool takes a "threads" option (like `-p` for HISAT2 or `-@` for samtools), we'll tell it to use 2. Asking a program for more threads than Slurm gave you doesn't make it faster, it just makes the threads compete for the same CPUs.

### Running out of memory

8 GB should be plenty for this week. If a command dies with a message like `Killed` or `oom-kill`, exit the container, re-launch it with `--mem=16G`, `source` your environment file, and try again.

### Practical exercises

Each lesson ends with a practical exercise that asks you to repeat what you just did on a second dataset. The solutions are on the [rnabio.org solutions page](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/). Try them yourself first!

---


## 6. Module 1: Inputs

![RNA-seq_Flowchart](https://rnabio.org/assets/module_1/RNA-seq_Flowchart2.png)

---

### 6.1 Reference genomes

#### Obtain a reference genome from Ensembl, iGenomes, NCBI or UCSC

In this example analysis we will use the human GRCh38 version of the genome from Ensembl. We are actually going to perform the analysis using only a single chromosome (chr22) and the ERCC spike-in to make it run faster.

First, make sure your environment variables are loaded (see [Section 4](#4-set-up-the-course-environment-variables)), then create the necessary working directory:

```bash
source /workspace/rnabio_env.sh
cd $RNA_HOME
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
```

The complete data from which these files were obtained can be found at: [ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/](ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/). You could use wget to download the `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz` file, then unzip/untar.

We've prepared a simplified reference that contains chr22 (and ERCC transcript) fasta files in both a single combined file and individual files. Download the reference genome file to the rnaseq working directory:

```bash
cd $RNA_REFS_DIR
wget https://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
ls
```

**Use your favorite LLM to help you answer the following questions:**

1) View the first 10 lines of this file. Why does it look like this?

```bash
head chr22_with_ERCC92.fa
```

2) How many lines and characters are in this file? How long is this chromosome (in bases and Mbp)?

```bash
wc chr22_with_ERCC92.fa
```

3) View 10 lines from approximately the middle of this file. What is the significance of the upper and lower case characters?

```bash
head -n 425000 chr22_with_ERCC92.fa | tail
```

4) What is the count of each base in the entire reference genome file (skipping the header lines for each sequence)?

```bash
cat chr22_with_ERCC92.fa | grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
```

<details>
<summary>Answers</summary>

- 1) After the `>22` header line, the sequence is all `N`s. Regions in GRCh38 that are unassembled are replaced with Ns. These are typically things like centromeres and telomeres.  Newer references from the T2T consortium may have resolution in these regions.
- 2) 51751056 is the count of bytes (characters) from `wc`. There are some newlines and header lines included there, but ~51 Mbp is a reasonable estimate of the chromosome's length.
- 3) Lowercase letters represent "soft-masking", generally of regions that are flagged by RepeatMasker. For things like alignment here, we can just ignore case.
- 4) That's a perl solution to the problem, but solutions in awk, python, or other such tools would be just as valid. On the command line, there is almost always more than one way to do something!

</details>


In your own work, instead of the above, you could also consider getting reference genomes and associated annotations from [UCSC, e.g. the UCSC GRCh38 download](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/).

Wherever you get them from, remember that the names of your reference sequences (chromosomes) must match those in your annotation GTF files (described in the next section).

View a list of all sequences in our reference genome fasta file:

```bash
grep ">" chr22_with_ERCC92.fa
```

---

### 6.2 Annotations

#### Obtain known gene/transcript annotations

In this tutorial we will use annotations obtained from Ensembl ([Homo_sapiens.GRCh38.86.gtf.gz](ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz)) for chromosome 22 only. For time reasons, these have been prepared for you. But you should get familiar with sources of gene annotations for RNA-seq analysis.

Download the gene annotation file to the working directory:

```bash
echo $RNA_REFS_DIR
cd $RNA_REFS_DIR
ls

wget https://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf
```

Take a look at the contents of the GTF file. Press `q` to exit the `less` display.

```bash
echo $RNA_REF_GTF
less start_codon -S $RNA_REF_GTF
```

Note how the `-S` option makes it easier to view this file with `less`. Make the formatting a bit nicer still:

```bash
cat chr22_with_ERCC92.gtf | column -t | less -S
```

How many unique gene IDs are in the GTF file? We can use a perl command-line command to find out:

```bash
perl -ne 'if ($_ =~ /(gene_id\s\"ENSG\w+\")/){print "$1\n"}' $RNA_REF_GTF | sort | uniq | wc -l
```

- Using `perl -ne ''` will execute the code between single quotes, on the GTF file, line-by-line.
- The `$_` variable holds the contents of each line.
- The `'if ($_ =~//)'` is a pattern-matching command which will look for the pattern "gene_id" followed by a space followed by "ENSG" and one or more word characters (indicated by `\w+`) surrounded by double quotes.
- The pattern to be matched is enclosed in parentheses. This allows us to print it out from the special variable `$1`.
- The output of this perl command will be a long list of ENSG IDs.
- By piping to `sort`, then `uniq`, then word count we can count the unique number of genes in the file.

We can also use `grep` to find this same information:

```bash
cat chr22_with_ERCC92.gtf | grep -w gene | wc -l
```

- `grep -w gene` is telling grep to do an exact match for the string 'gene'. This means that it will return lines that are of the feature type `gene`.

Now view the structure of a single transcript in GTF format. Press `q` to exit the `less` display when you are done.

```bash
grep ENST00000342247 $RNA_REF_GTF | less -p "exon\s" -S
```

To learn more, see:

- [https://perldoc.perl.org/perlre.html#Regular-Expressions](https://perldoc.perl.org/perlre.html#Regular-Expressions)
- [https://www.perl.com/pub/2004/08/09/commandline.html](https://www.perl.com/pub/2004/08/09/commandline.html)

#### Definitions

**Reference genome** - The nucleotide sequence of the chromosomes of a species. Genes are the functional units of a reference genome and gene annotations describe the structure of transcripts expressed from those gene loci.

**Gene annotations** - Descriptions of gene/transcript models for a genome. A transcript model consists of the coordinates of the exons of a transcript on a reference genome. Additional information such as the strand the transcript is generated from, gene name, coding portion of the transcript, alternate transcript start sites, and other information may be provided.

**GTF (.gtf) file** - A common file format referred to as Gene Transfer Format used to store gene and transcript annotation information. You can learn more about this format here: [https://genome.ucsc.edu/FAQ/FAQformat#format4](https://genome.ucsc.edu/FAQ/FAQformat#format4)

#### The purpose of gene annotations (.gtf file)

When running the HISAT2/StringTie/Ballgown pipeline, known gene/transcript annotations are used for several purposes:

- During the HISAT2 index creation step, annotations may be provided to create local indexes to represent transcripts as well as a global index for the entire reference genome. This allows for faster mapping and better mapping across exon boundaries and splice sites. If an alignment still can not be found it will attempt to determine if the read corresponds to a novel exon-exon junction. See the [Indexing section](#63-indexing) and the HISAT2 publication for more details.
- During the StringTie step, a GTF file can be used to specify transcript models to guide the assembly process and limit expression estimates to predefined transcripts using the `-G` and `-e` options together. The `-e` option will give you one expression estimate for each of the transcripts in your GTF file, giving you a 'microarray like' expression result.
- During the StringTie step, if the `-G` option is specified without the `-e` option the GTF file is used only to 'guide' the assembly of transcripts. Instead of assuming only the known transcript models are correct, the resulting expression estimates will correspond to both known and novel/predicted transcripts.
- During the StringTie and gffcompare steps, a GTF file is used to determine the transcripts that will be examined for differential expression using Ballgown. These may be known transcripts that you download from a public source, or a GTF of transcripts predicted by StringTie from the read data in an earlier step.

#### Sources for obtaining gene annotation files formatted for HISAT2/StringTie/Ballgown

There are many possible sources of GTF gene/transcript annotation files, for example from Ensembl, UCSC, RefSeq, etc. You don't need to download anything from these for this assignment, but you should know where to find them:

- **Ensembl FTP site:** Based on Ensembl annotations only. Available for many species. [https://useast.ensembl.org/info/data/ftp/index.html](https://useast.ensembl.org/info/data/ftp/index.html)
- **UCSC Table Browser:** Based on UCSC annotations or several other possible annotation sources collected by UCSC. You might choose this option if you want a lot of flexibility in the annotations you obtain (for example, grabbing only the GENCODE transcripts from chromosome 22). From [https://genome.ucsc.edu/](https://genome.ucsc.edu/), choose *Tools > Table Browser*, pick your genome, track, and region, and set the output format to *GTF - gene transfer format*. The same tool can also give you gene or exon BED files, and gene symbols and descriptions. See the [original lesson](https://rnabio.org/module-01-inputs/0001/03/01/Annotations/) for step-by-step instructions.
- **HISAT2 precomputed genome indexes:** HISAT2 has prebuilt reference genome index files for both DNA and RNA alignment. Various versions of the index files include SNPs and/or transcript splice sites. Versions of both the Ensembl and UCSC genomes for human build 38 are linked from the main HISAT2 page: [https://ccb.jhu.edu/software/hisat2/index.shtml](https://ccb.jhu.edu/software/hisat2/index.shtml)

#### Important notes

**On chromosome naming conventions:** In order for your RNA-seq analysis to work, the chromosome names in your GTF file must match those in your reference genome (i.e. your reference genome fasta file). If you get a StringTie result where all transcripts have an expression value of 0, you may have overlooked this. Unfortunately, Ensembl, NCBI, and UCSC can not agree on how to name the chromosomes in many species, so this problem may come up often. You can avoid this by getting a complete reference genome and gene annotation package from the same source (e.g., Ensembl) to maintain consistency.

**On reference genome builds:** Your annotations must correspond to the same reference genome build as your reference genome fasta file, e.g., both correspond to UCSC human build 'hg38', NCBI human build 'GRCh38', etc. Even if both your reference genome and annotations are from UCSC or Ensembl they could still correspond to different versions of that genome. This would cause problems in any RNA-seq pipeline.

A more detailed discussion of commonly used versions of the human reference genome can be found in a companion workshop: [PMBIO Reference Genomes](https://pmbio.org/module-02-inputs/0002/02/01/Reference_Genome/).

---

### 6.3 Indexing

Remember that indexing is the "lookup table" that the aligner uses to find matches very quickly. If you want a detailed refresher on indexing, the RNAbio folks have a [indexing mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2026/mini/RNASeq_MiniLecture_01_02_Indexing.pdf) available.

#### Create a HISAT2 index

Create a HISAT2 index for chr22 and the ERCC spike-in sequences. HISAT2 can incorporate exons and splice sites into the index file for alignment. First create a splice site file, then an exon file. Finally make the aligner FM index.

To learn more about how the HISAT2 indexing strategy is distinct from other next gen aligners refer to the [HISAT publication](https://www.ncbi.nlm.nih.gov/pubmed/25751142).

```bash
cd $RNA_REFS_DIR
hisat2_extract_splice_sites.py $RNA_REF_GTF > $RNA_REFS_DIR/splicesites.tsv
hisat2_extract_exons.py $RNA_REF_GTF > $RNA_REFS_DIR/exons.tsv
hisat2-build -p 2 --ss $RNA_REFS_DIR/splicesites.tsv --exon $RNA_REFS_DIR/exons.tsv $RNA_REF_FASTA $RNA_REF_INDEX
ls
```

(`-p 2` tells `hisat2-build` to use the 2 CPUs we requested. Building the index takes a few minutes.)

Perform a visual survey of the contents of your refs directory. What is the source/purpose of each file?

> [!NOTE]
> To create an index for all chromosomes instead of just chr22, you would do something like the following. **Don't run this!** To index the entire human genome, HISAT2 requires about 160 GB of RAM, and we only asked Slurm for 8 GB.
>
> ```bash
> #hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.86.gtf > splicesites.tsv
> #hisat2_extract_exons.py Homo_sapiens.GRCh38.86.gtf > exons.tsv
> #hisat2-build -p 2 --ss splicesites.tsv --exon exons.tsv Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa Homo_sapiens.GRCh38.dna_sm.primary_assembly
> ```

---

### 6.4 RNA-seq data

#### Obtain RNA-seq test data

The test data consists of two commercially available RNA samples: [Universal Human Reference (UHR)](https://rnabio.org/assets/module_1/UHR.pdf) and [Human Brain Reference (HBR)](https://rnabio.org/assets/module_1/HBR.pdf). The UHR is total RNA isolated from a diverse set of 10 cancer cell lines (breast, liver, cervix, testis, brain, skin, fatty tissue, histocyte, macrophage, T cell, B cell). The HBR is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old.

In addition, a spike-in control was used. Specifically, an aliquot of the [ERCC ExFold RNA Spike-In Control Mixes](https://rnabio.org/assets/module_1/ERCC.pdf) was added to each sample. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies). This range allows us to test the degree to which the RNA-seq assay (including all laboratory and analysis steps) accurately reflects the relative abundance of transcript species within a sample. There are two 'mixes' of these transcripts to allow an assessment of differential expression output between samples if you put one mix in each of your two comparisons. In our case, Mix1 was added to the UHR sample, and Mix2 was added to the HBR sample. We also have 3 complete experimental replicates for each sample. This allows us to assess the technical variability of our overall process of producing RNA-seq data in the lab.

From the ERCC product sheet:

> The ERCC plasmid reference library produces well-characterized transcripts generated largely from random unique sequences. Sequence comparisons have been made to multiple databases available at the time of design including mouse, rat, human, Drosophila, bacteria, mosquito, and other nonhuman species. The control collection contains some sequences with homology to Bacillus subtilis.

For all libraries, low-throughput (Set A) TruSeq Stranded Total RNA Sample Prep Kit libraries were prepared with Ribo-Zero Gold to remove both cytoplasmic and mitochondrial rRNA. Triplicate, indexed libraries were made starting with 100ng Agilent/Strategene Universal Human Reference total RNA and 100ng Ambion Human Brain Reference total RNA. The Universal Human Reference replicates received 2 ul of 1:1000 ERCC Mix 1. The Human Brain Reference replicates received 2 ul of 1:1000 ERCC Mix 2. The libraries were quantified with KAPA Library Quantification qPCR and adjusted to the appropriate concentration for sequencing. The triplicate, indexed libraries were then pooled prior to sequencing. Each pool of three replicate libraries was sequenced across 2 lanes of a HiSeq 2000 using paired-end sequence chemistry with 100 bp read lengths.

So to summarize we have:

- UHR + ERCC Spike-In Mix1, Replicate 1
- UHR + ERCC Spike-In Mix1, Replicate 2
- UHR + ERCC Spike-In Mix1, Replicate 3
- HBR + ERCC Spike-In Mix2, Replicate 1
- HBR + ERCC Spike-In Mix2, Replicate 2
- HBR + ERCC Spike-In Mix2, Replicate 3

Each data set has a corresponding pair of FASTQ files (read 1 and read 2 of paired end reads).

```bash
echo $RNA_DATA_DIR
mkdir -p $RNA_DATA_DIR
cd $RNA_DATA_DIR
wget https://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
```

Unpack the test data using `tar`. You should see 6 sets of paired end fastq files, one for each of our sample replicates above. We have 6 pairs (12 files) because in fastq format, read 1 and read 2 of each read pair (fragment) are stored in separate files.

```bash
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls
```

The reads are paired-end 101-mers generated on an Illumina HiSeq instrument. The test data has been pre-filtered for reads that appear to map to chromosome 22. View the first two read records of a file (in fastq format each read corresponds to 4 lines of data):

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8
```

Identify the following components of each read: read name, read sequence, and quality string.

How many reads are there in the first library? Decompress the file on the fly with `zcat`, pipe into `grep`, search for the read name prefix and pipe into `wc` to do a word count (`-l` gives lines):

```bash
zcat UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | grep -P "^\@HWI" | wc -l
```

#### Determining the strandedness of RNA-seq data (optional)

To determine strandedness, we will be using [check_strandedness](https://github.com/betsig/how_are_we_stranded_here). To use this tool, there are a few steps we need to get our inputs ready, specifically creating a fasta of our GTF file.

To create a transcripts FASTA file from our reference transcriptome in GTF format we need to get our transcript information into BED12 format. To get to that format we will convert GTF -> GenePred -> BED12. The GenePred format is just an intermediate that we create because there isn't a tool that converts straight from GTF to BED12.

Once we have our transcriptome information (the position of exons for each transcript on each chromosome) we can use `bedtools getfasta` to splice the exon sequences from the chromosome sequence together into a full length transcript sequence.

```bash
cd $RNA_HOME/refs/

# Convert our reference transcriptome GTF to genePred format
gtfToGenePred chr22_with_ERCC92.gtf chr22_with_ERCC92.genePred

# Convert the genePred format to bed12 format
genePredToBed chr22_with_ERCC92.genePred chr22_with_ERCC92.bed12

# Use bedtools to create fasta from GTF
bedtools getfasta -fi chr22_with_ERCC92.fa -bed chr22_with_ERCC92.bed12 -s -split -name -fo chr22_ERCC92_transcripts.fa
```

Use `less` to view the file `chr22_ERCC92_transcripts.fa`. Note that this file has messy transcript names. Use the following hairball perl one-liner to tidy up the header line for each fasta sequence:

```bash
cd $RNA_HOME/refs
cat chr22_ERCC92_transcripts.fa | perl -ne 'if($_ =~/^\>\S+\:\:(ERCC\-\d+)\:.*/){print ">$1\n"}elsif ($_ =~/^\>(\S+)\:\:.*/){print ">$1\n"}else{print $_}' > chr22_ERCC92_transcripts.clean.fa
```

View the resulting 'clean' file using `less chr22_ERCC92_transcripts.clean.fa` (use `q` to exit). View the end of this file using `tail chr22_ERCC92_transcripts.clean.fa`. Note that we have one fasta record for each Ensembl transcript on chromosome 22 and we have an additional fasta record for each ERCC spike-in sequence.

We also need to reformat our GTF file slightly. Rows that correspond to genes are missing the "transcript_id" field. We are going to add in this field but leave it blank for these rows using the following command:

```bash
cd $RNA_HOME/refs
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' chr22_with_ERCC92.gtf > chr22_with_ERCC92_tidy.gtf
```

Now that we have created our input files, we can run the check_strandedness tool on some of our instrument data. `check_strandedness` (along with the `kallisto` and RSeQC tools it relies on) is already installed in our toolkit container, so we can run it right here:

```bash
cd $RNA_HOME
check_strandedness --gtf $RNA_REFS_DIR/chr22_with_ERCC92_tidy.gtf --transcripts $RNA_REFS_DIR/chr22_ERCC92_transcripts.clean.fa --reads_1 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads_2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```

`check_strandedness` writes its working files (a `kallisto_index` file and a `stranded_test_...` folder) into whatever directory you run it from, which is why we `cd $RNA_HOME` first.

The output of this command should look like this:

```
Fraction of reads failed to determine: 0.1123
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0155
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8722
Over 75% of reads explained by "1+-,1-+,2++,2--"
Data is likely RF/fr-firststrand
```

The checkstrandedness test asks a single question: in each read pair, which mate aligns in the same orientation as the gene it overlaps? The tool encodes every read as `read#` + `aligned strand` + `gene strand`, then sorts reads into two patterns. For the 'reverse' pattern `1+-,1-+,2++,2--`: read 2 always agrees with the gene's strand (2++, 2--) while read 1 always disagrees (1+-, 1-+). The 'forward' pattern `1++,1--,2+-,2-+` is the exact mirror, with read 1 always agreeing and read 2 always disagreeing. The four codes within each pattern aren't independent observations; they're that one orientation rule applied to genes on the + strand and the - strand, which is why all four always appear together.

In a dUTP/fr-firststrand prep, first-strand cDNA is synthesized antisense to the mRNA, the dUTP-marked second strand is degraded before amplification, and the adapters are arranged so read 1 sequences the surviving antisense strand while read 2 reads its sense complement. So read 2 matches the transcript and read 1 opposes it. The two names describe the same molecule from different angles: "first-strand" refers to the cDNA strand that is retained, and "reverse" refers to read 1 pointing opposite the transcript.

In this case, 87.22% of reads fall in the reverse pattern, only 1.55% in the forward pattern, and 11.23% are undetermined (overlapping genes, multimappers, and ambiguous orientation, all expected and harmless). This strong bias toward `1+-,1-+,2++,2--` suggests an unambiguous RF / fr-firststrand call. A forward-stranded library would have put the bulk in the other bin, and an unstranded library would have split close to 50:50.

Using this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/), we can see if this result is what we expect. Since the UHR and HBR data were generated with the TruSeq Stranded Kit, as mentioned above, the correct strand setting for kallisto is `--rf-stranded`, which is what check_strandedness confirms. Similarly, when we run HISAT2 we will use `--rna-strandness RF`, when we run StringTie we will use `--rf`, and when we run htseq-count we will use `--stranded reverse`.

#### Practical exercise 3

Download an additional dataset and unpack it. This data will be used in future practical exercises.

- Hint: Do this in a separate working directory called `practice` (i.e. `$RNA_HOME/practice`), and create sub-directories for organization (data, alignments, etc).
- In this exercise you will download an archive of publicly available read data from here: [https://genomedata.org/rnaseq-tutorial/practical.tar](https://genomedata.org/rnaseq-tutorial/practical.tar)

The practice dataset includes 3 replicates of data from the HCC1395 breast cancer cell line and 3 replicates of data from the HCC1395BL matched lymphoblastoid line. So, this will be a tumor vs normal (cell line) comparison. The reads are paired-end 151-mers generated on an Illumina HiSeq instrument from TruSeq Stranded Total RNA Libraries. The test data has been pre-filtered for reads that appear to map to chromosome 22.

**Questions:**

- How many data files were contained in the `practical.tar` archive? What commonly used sequence data file format are they?
- In the first read of the hcc1395, normal, replicate 1, read 1 file, what was the physical location of the read on the flow cell (i.e. lane, tile, x, y)?
- In the first read of this same file, how many 'T' bases are there?

When you are ready you can check your approach against the [solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-3---data).

---

### 6.5 Pre-alignment QC

#### FastQC

In our [sequence data](../week_03) exercise we used `FastQC` to look at DNA sequence data.  Here, we show that you can also use it to get a sense of your **RNA** data quality before alignment. 

More info on FastQC:
- [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Video tutorial: [https://www.youtube.com/watch?v=bz93ReOv87Y](https://www.youtube.com/watch?v=bz93ReOv87Y)

Run FastQC on your fastq files:

```bash
cd $RNA_HOME/data
fastqc *.fastq.gz
```

Then, on your laptop, open the mounted storage (see [Viewing HTML reports](#viewing-html-reports)), go to `rnabio-workspace/rnaseq/data/`, and open any of the `*_fastqc.html` files to view the FastQC report.

**Exercise:** Investigate the source/explanation for over-represented sequences. (Hint: try BLASTing them.)

#### Fastp

`Fastp` is a similar alternative tool. QC results for this tool can be produced as follows (this is a good block to copy/paste):

```bash
cd $RNA_HOME/data
mkdir -p fastp

fastp -i $RNA_HOME/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/HBR_Rep1.fastp.json --html fastp/HBR_Rep1.fastp.html 2>fastp/HBR_Rep1.fastp.log

fastp -i $RNA_HOME/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/HBR_Rep2.fastp.json --html fastp/HBR_Rep2.fastp.html 2>fastp/HBR_Rep2.fastp.log

fastp -i $RNA_HOME/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/HBR_Rep3.fastp.json --html fastp/HBR_Rep3.fastp.html 2>fastp/HBR_Rep3.fastp.log

fastp -i $RNA_HOME/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/UHR_Rep1.fastp.json --html fastp/UHR_Rep1.fastp.html 2>fastp/UHR_Rep1.fastp.log

fastp -i $RNA_HOME/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/UHR_Rep2.fastp.json --html fastp/UHR_Rep2.fastp.html 2>fastp/UHR_Rep2.fastp.log

fastp -i $RNA_HOME/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -I $RNA_HOME/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz --json fastp/UHR_Rep3.fastp.json --html fastp/UHR_Rep3.fastp.html 2>fastp/UHR_Rep3.fastp.log
```

The fastp HTML reports will be in `rnabio-workspace/rnaseq/data/fastp/` on your mounted storage.

#### MultiQC

Run MultiQC on your FastQC reports to generate a single summary report across all samples/replicates:

```bash
cd $RNA_HOME/data
multiqc ./
```

Then open `rnabio-workspace/rnaseq/data/multiqc_report.html` from your mounted storage.

#### Clean up

Move all the FastQC files into their own directory:

```bash
cd $RNA_HOME/data
mkdir fastqc
mv *_fastqc* fastqc
```

#### Practical exercise 4

Run FastQC on one of the additional fastq files you downloaded in the previous practical exercise. (Hint: remember that you stored this data in a separate working directory called `practice`.)

Run FastQC on the file `hcc1395_normal_1.fastq.gz` and answer these questions by examining the output:

- How many total sequences are there?
- What is the range (x - y) of read lengths observed?
- What is the most common average sequence quality score?
- What does the Adapter Content warning tell us?

When you are ready you can check your approach against the [solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-4---data-qc).

---

## 7. Module 2: Alignment

![RNA-seq_Flowchart3](https://rnabio.org/assets/module_2/RNA-seq_Flowchart3.png)

### Key concepts

- Splice-aware alignment, SAM/BAM format, SAM flags, CIGAR strings, SAM sorting, etc.

### Learning objectives

- RNA-seq alignment challenges and common questions
- Alignment strategies
- HISAT2
- Introduction to the BAM and BED formats
- Basic manipulation of BAMs
- Alignment QC assessment
- BAM read counting and determination of variant allele expression status


### 7.1 Adapter trimming (optional RNAseq step)

Use fastp to trim sequence adapters from the read FASTQ files and also perform basic data quality cleanup. The output of this step will be trimmed and filtered FASTQ files for each data set. (We'll align the *untrimmed* reads in the next section, so this step is optional, but it's good practice.)

Refer to the fastp project and manual for a more detailed explanation:

- [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)

Fastp basic usage:

```bash
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```

The `-i` and `-I` parameters specify input R1 and R2 data files (raw data). The `-o` and `-O` parameters specify output R1 and R2 data files (trimmed and quality filtered). Extra options specified below:

- `-l 25`: the minimum read length allowed after trimming is 25bp
- `--adapter_fasta`: the path to the adapter FASTA file containing adapter sequences to trim
- `--trim_front1 13`: trim a fixed number (13 in this case) of bases off the left end of read1
- `--trim_front2 13`: trim a fixed number (13 in this case) of bases off the left end of read2
- `--json`: the path to store a log file in JSON file format
- `--html`: the path to store a web report file
- `2>`: use to store the information that would be printed to the screen into a file instead

#### Read trimming with fastp

First, set up some directories for output:

```bash
echo $RNA_DATA_TRIM_DIR
mkdir -p $RNA_DATA_TRIM_DIR
```

Download the necessary Illumina adapter sequence files:

```bash
echo $RNA_REFS_DIR
mkdir -p $RNA_REFS_DIR
cd $RNA_REFS_DIR
wget https://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
```

Use fastp to remove Illumina adapter sequences (if any), trim the first 13 bases of each read, and perform default read quality filtering to remove reads that are too short, have too many low quality bases or have too many N's:

```bash
cd $RNA_HOME

export S1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S1}.read1.fastq.gz -I $RNA_DATA_DIR/${S1}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S1}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S1}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S1}.fastp.json --html $RNA_DATA_TRIM_DIR/${S1}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S1}.fastp.log

export S2=UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S2}.read1.fastq.gz -I $RNA_DATA_DIR/${S2}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S2}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S2}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S2}.fastp.json --html $RNA_DATA_TRIM_DIR/${S2}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S2}.fastp.log

export S3=UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S3}.read1.fastq.gz -I $RNA_DATA_DIR/${S3}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S3}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S3}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S3}.fastp.json --html $RNA_DATA_TRIM_DIR/${S3}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S3}.fastp.log

export S4=HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S4}.read1.fastq.gz -I $RNA_DATA_DIR/${S4}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S4}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S4}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S4}.fastp.json --html $RNA_DATA_TRIM_DIR/${S4}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S4}.fastp.log

export S5=HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S5}.read1.fastq.gz -I $RNA_DATA_DIR/${S5}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S5}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S5}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S5}.fastp.json --html $RNA_DATA_TRIM_DIR/${S5}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S5}.fastp.log

export S6=HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22
fastp -i $RNA_DATA_DIR/${S6}.read1.fastq.gz -I $RNA_DATA_DIR/${S6}.read2.fastq.gz -o $RNA_DATA_TRIM_DIR/${S6}.read1.fastq.gz -O $RNA_DATA_TRIM_DIR/${S6}.read2.fastq.gz -l 25 --adapter_fasta $RNA_REFS_DIR/illumina_multiplex.fa --trim_front1 13 --trim_front2 13 --json $RNA_DATA_TRIM_DIR/${S6}.fastp.json --html $RNA_DATA_TRIM_DIR/${S6}.fastp.html 2>$RNA_DATA_TRIM_DIR/${S6}.fastp.log
```

#### Use FastQC and MultiQC to compare the impact of trimming

Compare the FastQC reports for fastq files before and after trimming. All FastQC reports can be generated on the command line:

```bash
cd $RNA_DATA_TRIM_DIR
fastqc *.fastq.gz

multiqc ./
```

To view the resulting HTML reports, open these folders on your mounted storage and compare the reports:

- `rnabio-workspace/rnaseq/data/fastqc/` (before trimming)
- `rnabio-workspace/rnaseq/data/trimmed/` (after trimming)

#### Clean up

Move the FastQC and fastp results into sub-directories to keep things tidy:

```bash
cd $RNA_DATA_TRIM_DIR
mkdir fastqc
mv *_fastqc* fastqc
mkdir fastp
mv *fastp.* fastp
```

#### Practical exercise 5

Using the approach above, trim the reads for both normal and tumor samples that you downloaded for the previous practical exercise. Try dropping the hard left trim option used above (`--trim_front1 13` and `--trim_front2 13`). Once you have trimmed the reads, compare a pre- and post-trimming FastQ file using the FastQC and MultiQC tools.

- Hint: These files should have been downloaded to `$RNA_HOME/practice/data/`.

Answer these questions by examining the FastQC reports:

- After trimming, what is the range of read lengths observed for hcc1395 normal replicate 1, read 1?
- Which sections of the FastQC report are most informative for observing the effect of trimming?
- In the 'Per base sequence content' section, what pattern do you see? What could explain this pattern?

When you are ready you can check your approach against the [solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-5---trim).

---

### 7.2 [Alignment](https://rnabio.org/module-02-alignment/0002/03/01/Alignment/)

If you would like a refresher on alignment, the course authors have created an [alignment mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2026/mini/RNASeq_MiniLecture_02_01_Alignment.pdf). There are also mini lectures describing the differences between [alignment, assembly, and pseudoalignment](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2026/mini/RNASeq_MiniLecture_02_02_Alignment_vs_Assembly_vs_Kmer.pdf) and describing [SAM, BAM, and BED file formats](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2026/mini/RNASeq_MiniLecture_02_03_SAM_BAM_BED.pdf).

#### HISAT2 alignment

Perform alignments with HISAT2 to the genome and transcriptome.

First, begin by making the appropriate output directory for our alignment results:

```bash
echo $RNA_ALIGN_DIR
mkdir -p $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR
```

HISAT2 uses a graph-based alignment and has succeeded HISAT and TOPHAT2. The output of this step will be a SAM/BAM file for each data set.

Refer to the HISAT2 manual for a more detailed explanation:

- [https://daehwankimlab.github.io/hisat2/manual/](https://daehwankimlab.github.io/hisat2/manual/)

HISAT2 basic usage:

```bash
#hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
```

Extra options specified below:

- `-p 2` tells HISAT2 to use 2 CPUs (the number we asked Slurm for).
- `--rna-strandness RF` specifies the strandedness of the RNA-seq library. We specify RF since the TruSeq strand-specific library was used to make these libraries (and since that's what check_strandedness told us above).
- `--rg-id $ID` specifies a read group ID that is a unique identifier.
- `--rg SM:$SAMPLE_NAME` specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
- `--rg LB:$LIBRARY_NAME` specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
- `--rg PL:ILLUMINA` specifies a read group sequencing platform.
- `--rg PU:$PLATFORM_UNIT` specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
- `--dta` reports alignments tailored for transcript assemblers.
- `-x /path/to/hisat2/index` is the HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier, including splice sites and exons.
- `-1 /path/to/read1.fastq.gz` is the read 1 FASTQ file, optionally gzip (.gz) or bzip2 (.bz2) compressed.
- `-2 /path/to/read2.fastq.gz` is the read 2 FASTQ file, optionally gzip (.gz) or bzip2 (.bz2) compressed.
- `-S /path/to/output.sam` is the output SAM format text file of alignments.

```bash
hisat2 -p 2 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep1.sam
hisat2 -p 2 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep2.sam
hisat2 -p 2 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep3.sam

hisat2 -p 2 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep1.sam
hisat2 -p 2 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep2.sam
hisat2 -p 2 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep3.sam
```

Note: in the above alignments, we are treating each library as an independent data set. If you had multiple lanes of data for a single library, you could align them all together in one HISAT2 command. Similarly you might combine technical replicates into a single alignment run (perhaps after examining them and removing outliers...). To combine multiple lanes, you would provide all the read1 files as a comma separated list for the `-1` input argument, and then all read2 files as a comma separated list for the `-2` input argument (where both lists have the same order). You can also use samtools merge to combine bam files after alignment. This is the approach we will take.

**HISAT2 alignment summary:** HISAT2 generates a summary of the alignments printed to the terminal (on STDERR). Notice the number of total reads, reads aligned, and various metrics regarding how the reads aligned to the reference.

#### SAM to BAM conversion

Convert the HISAT2 SAM files to BAM files and sort by aligned position:

```bash
samtools sort -@ 2 -o UHR_Rep1.bam UHR_Rep1.sam
samtools sort -@ 2 -o UHR_Rep2.bam UHR_Rep2.sam
samtools sort -@ 2 -o UHR_Rep3.bam UHR_Rep3.sam
samtools sort -@ 2 -o HBR_Rep1.bam HBR_Rep1.sam
samtools sort -@ 2 -o HBR_Rep2.bam HBR_Rep2.sam
samtools sort -@ 2 -o HBR_Rep3.bam HBR_Rep3.sam
```

#### Merge HISAT2 BAM files

Make a single BAM file combining all UHR data and another for all HBR data. This could be done in several ways such as `samtools merge`, `bamtools merge`, or using picard-tools (see below). Here we're choosing the third method because it does the best job of merging the BAM header information. (Note: sambamba also retains header info.)

```bash
cd $RNA_HOME/alignments/hisat2
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT UHR.bam -INPUT UHR_Rep1.bam -INPUT UHR_Rep2.bam -INPUT UHR_Rep3.bam
java -Xmx2g -jar $PICARD MergeSamFiles -OUTPUT HBR.bam -INPUT HBR_Rep1.bam -INPUT HBR_Rep2.bam -INPUT HBR_Rep3.bam
```

(`$PICARD` is one of the environment variables we set in [Section 4](#4-set-up-the-course-environment-variables). It points to the copy of Picard that comes pre-installed in the docker image.)

Count the alignment (BAM) files to make sure all were created successfully (you should have 8 total):

```bash
ls -l *.bam | wc -l
ls -l *.bam
```

#### Practical exercise 6

Perform some alignments on additional read data sets. Align the reads using the skills you learned above. Try using the HISAT2 aligner. Also practice converting SAM to BAM files, and merging BAM files.

- Hint: Do this analysis on the additional data, and in the separate working directory called `practice` that you created in Practical Exercise 3.

**Questions:**

- What is the difference between a .sam and .bam file?
- If you sorted the resulting BAM file as we did above, is the result sorted by read name? Or position?
- Which columns of the BAM file can be viewed to determine the style of sorting?
- What command can you use to view only the BAM header?

When you are ready you can check your approach against the [solutions](https://rnabio.org/module-09-appendix/0009/05/01/Practical_Exercise_Solutions/#practical-exercise-6---alignment).

---

### 7.3 [Indexing BAM files](https://rnabio.org/module-02-alignment/0002/05/01/AlignVis_Prep/)

Before we can look at specific regions of our alignments (with samtools, or later with IGV), we need to index our BAM files. We will use `samtools index` for this purpose. For convenience later, index all of the BAM files:

```bash
echo $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR

samtools index -M *.bam
```

The `-M` flag interprets all filename arguments as files to be indexed, allowing multiple files to be indexed at the same time. To index an individual file, use `samtools index input.bam`. Note that we could also have created and run a `samtools index` command for every file ending in `.bam` using the following construct:

```bash
# for i in *.bam;do samtools index $i;done
```
(we also could have used `xargs` - remember, there's more than one way to do it!)

---

### 7.4 [Allele expression](https://rnabio.org/module-02-alignment/0002/05/03/AlignVis_AlleleExpression/)

In this section we will demonstrate how to assess expression of specific variant alleles in the RNA-seq BAM using tools designed to interrogate read alignments and sequence base identities at particular positions.

#### BAM read counting

Using a known variant position, count the number of supporting reference and variant reads. First, use `samtools mpileup` to visualize a region of alignment with a variant.

```bash
cd $RNA_HOME
mkdir bam_readcount
cd bam_readcount
```

Create a faidx indexed reference sequence file for use with mpileup:

```bash
echo $RNA_REF_FASTA
samtools faidx $RNA_REF_FASTA
```

Run `samtools mpileup` on a region of interest:

```bash
samtools mpileup -f $RNA_REF_FASTA -r 22:18918457-18918467 $RNA_ALIGN_DIR/UHR.bam $RNA_ALIGN_DIR/HBR.bam
```

Each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities. At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN` for a mismatch on the forward strand and `acgtn` for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+` indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. See the samtools pileup/mpileup documentation for more explanation of the output:

- [https://samtools.sourceforge.net/pileup.shtml](https://samtools.sourceforge.net/pileup.shtml)
- [https://samtools.sourceforge.net/mpileup.shtml](https://samtools.sourceforge.net/mpileup.shtml)

Now, use `bam-readcount` to count reference and variant bases at a specific position. First, create a bed file with some positions of interest (we will create a file called `snvs.bed` using the `echo` command).

It will contain a single line specifying a variant position on chr22, e.g.: `22:38483683-38483683`

Create the bed file:

```bash
echo "22 38483683 38483683"
echo "22 38483683 38483683" > snvs.bed
```

Run `bam-readcount` on this list for the UHR and HBR merged BAM files:

```bash
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/UHR.bam 2>/dev/null
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/HBR.bam 2>/dev/null
```

Now, run it again, but ignore STDERR and redirect STDOUT to a file:

```bash
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/UHR.bam 2>/dev/null 1>UHR_bam-readcounts.txt
bam-readcount -l snvs.bed -f $RNA_REF_FASTA $RNA_ALIGN_DIR/HBR.bam 2>/dev/null 1>HBR_bam-readcounts.txt
```

From this output, use AWK to pull out the read counts for each base (A, C, G, T) in columns 6-9 and then print a simplified summary. Each `split` command takes the specified column of bam-readcount data (6-9), parses the values separated by ":", and stores them in a `_data` array. The `printf` command prints out the summary line, including the total count for each base (from the second position of the data array for each nucleotide).

```bash
#UHR counts
awk '{
      split($6, A_data, ":")
      split($7, C_data, ":")
      split($8, G_data, ":")
      split($9, T_data, ":")
      printf "UHR Counts\t%s\t%s\tA: %s\tC: %s\tT: %s\tG: %s\n",
      $1, $2, A_data[2], C_data[2], T_data[2], G_data[2]
}' UHR_bam-readcounts.txt

#HBR counts
awk '{
      split($6, A_data, ":")
      split($7, C_data, ":")
      split($8, G_data, ":")
      split($9, T_data, ":")
      printf "HBR Counts\t%s\t%s\tA: %s\tC: %s\tT: %s\tG: %s\n",
      $1, $2, A_data[2], C_data[2], T_data[2], G_data[2]
}' HBR_bam-readcounts.txt
```

If you wish to extend this concept to a more complex analysis, here is a [bam-readcount tutorial](https://github.com/genome/bam-readcount/tree/master/tutorial) that uses python to parse output from bam-readcount to identify an Omicron SARS-CoV-2 variant of concern from raw sequence data.

---

### 7.5 [Alignment QC](https://rnabio.org/module-02-alignment/0002/06/01/Alignment_QC/)

If you would like a refresher on alignment QC, the RNAbio authors have made a [mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cbw/2026/mini/RNASeq_MiniLecture_02_04_alignmentQC.pdf) briefly covering the topic.

#### Use samtools and FastQC to evaluate the alignments

Use `samtools view` to see the format of a SAM/BAM alignment file:

```bash
cd $RNA_ALIGN_DIR
samtools view -H UHR.bam
samtools view UHR.bam | head
samtools view UHR.bam | head | column -t | less -S
```

Try filtering the BAM file to require or exclude certain flags. This can be done with the `samtools view -f -F` options:

- `-f INT` required flag
- `-F INT` filtering flag

"Samtools flags explained":

- [https://broadinstitute.github.io/picard/explain-flags.html](https://broadinstitute.github.io/picard/explain-flags.html)

Try requiring that alignments are 'paired' and 'mapped in a proper pair' (=3). Also filter out alignments that are 'unmapped', the 'mate is unmapped', and 'not primary alignment' (=268):

```bash
samtools view -f 3 -F 268 UHR.bam | head | column -t | less -S
```

Now require that the alignments be only for 'PCR or optical duplicate'. How many reads meet this criteria? Why?

```bash
samtools view -f 1024 UHR.bam | head
```

Use `samtools flagstat` to get a basic summary of an alignment. What percent of reads are mapped? Is this realistic? Why?

```bash
cd $RNA_ALIGN_DIR
mkdir flagstat

samtools flagstat HBR_Rep1.bam > flagstat/HBR_Rep1.bam.flagstat
samtools flagstat HBR_Rep2.bam > flagstat/HBR_Rep2.bam.flagstat
samtools flagstat HBR_Rep3.bam > flagstat/HBR_Rep3.bam.flagstat
samtools flagstat UHR_Rep1.bam > flagstat/UHR_Rep1.bam.flagstat
samtools flagstat UHR_Rep2.bam > flagstat/UHR_Rep2.bam.flagstat
samtools flagstat UHR_Rep3.bam > flagstat/UHR_Rep3.bam.flagstat

# Note that we could have created and run a samtools flagstat command for all files ending in *Rep*.bam using the following construct:
# find *Rep*.bam -exec echo samtools flagstat {} \> flagstat/{}.flagstat \; | sh

# View an example
cat flagstat/UHR_Rep1.bam.flagstat
```

Details of the SAM/BAM format can be found here: [https://samtools.github.io/hts-specs/SAMv1.pdf](https://samtools.github.io/hts-specs/SAMv1.pdf)

#### Create versions of our BAM files with only the chromosome 22 alignments

Note that our alignments contain those for the spiked-in ERCC control sequences. For some QC analyses, these sequences complicate interpretation because they are not like human genes (e.g. no introns or UTRs). To allow us flexibility for the QC analysis, we will create a version of each BAM file that only has the chromosome 22 alignments.

```bash
cd $RNA_ALIGN_DIR
mkdir chr22_only_bams
find *Rep*.bam -exec echo samtools view -b {} 22 -o chr22_only_bams/{} \; | sh
samtools index -M chr22_only_bams/*.bam
```

#### Using FastQC

You can use FastQC to perform basic QC of your BAM file (see [Pre-alignment QC](#65-pre-alignment-qc)). This will give you output very similar to when you ran FastQC on your fastq files.

```bash
cd $RNA_ALIGN_DIR
fastqc UHR_Rep1.bam UHR_Rep2.bam UHR_Rep3.bam HBR_Rep1.bam HBR_Rep2.bam HBR_Rep3.bam
mkdir fastqc
mv *fastqc.html fastqc/
mv *fastqc.zip fastqc/
```

The reports will be in `rnabio-workspace/rnaseq/alignments/hisat2/fastqc/` on your mounted storage.

#### Using Picard

You can use Picard to generate RNA-seq specific quality metrics and figures.

In this section we need to create some additional formats of our reference transcriptome files.

Picard uses a "sequence dictionary" file for many commands (simply a list of reference sequences and their sizes).

We will also filter our transcriptome GTF to one with only ribosomal features, convert it to BED format and then to IntervalList format. This is all done to get the IntervalList format needed for Picard CollectRnaSeqMetrics.

We will also create a version of our whole transcriptome GTF in the [RefFlat](https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat) format needed for Picard CollectRnaSeqMetrics. To get to the RefFlat format we will convert GTF to [GenePredExt](https://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#GenePredExt) format and then simplify this to RefFlat.

```bash
# Generating the necessary input files for picard CollectRnaSeqMetrics
cd $RNA_HOME/refs

# Create a .dict file for our reference
java -jar $PICARD CreateSequenceDictionary -R chr22_with_ERCC92.fa -O chr22_with_ERCC92.dict

# Create a bed file of the location of ribosomal sequences in our reference (first extract them from the GTF then convert to BED format)
# Note that here we pull all the "rrna" transcripts from the GTF. This is a good strategy for the whole transcriptome ...
# ... but on chr22 there is very little "rrna" content, leading to 0 coverage for all samples, so we are also adding a single protein coding ribosomal gene "RRP7A" (normally we would not do this)
# Note the convert2bed command will convert our GTF to BED format
# "<" is used to feed the GTF file into the tool.  "2>/dev/null" is used to throw away a harmless warning. "1>" is use to save our result to a file
grep --color=none -i -P "rrna|rrp7a" chr22_with_ERCC92.gtf > ref_ribosome.gtf
convert2bed --input=gtf --output=bed < ref_ribosome.gtf 2>/dev/null 1>ref_ribosome.bed

# Create interval list file for the location of just the ribosomal sequences in our reference
java -jar $PICARD BedToIntervalList -I ref_ribosome.bed -O ref_ribosome.interval_list -SD chr22_with_ERCC92.dict

# Create a genePred file for our whole reference transcriptome
gtfToGenePred -genePredExt chr22_with_ERCC92.gtf chr22_with_ERCC92.genePredExt

# Reformat this genePred file to first add the Ensembl gene ID column to the beginning of the dataframe using "awk", and then subset it down to the first 11 columns using "cut".
cat chr22_with_ERCC92.genePredExt | awk '{print $12"\t"$0}' | cut -d$'\t' -f1-11 > chr22_with_ERCC92.refFlat.txt

# Use the "find" command to run "picard CollectRnaSeqMetrics" on all 6 BAM files.
# The basic structure of this kind of automation is: find <search pattern> -exec command {} \;
# The "{}" will insert the file found by the "find" command using <search pattern>.  "\;" indicates the end of the command.
cd $RNA_HOME/alignments/hisat2/
mkdir picard
find *Rep*.bam -exec echo java -jar $PICARD CollectRnaSeqMetrics I={} O=picard/{}.RNA_Metrics REF_FLAT=$RNA_HOME/refs/chr22_with_ERCC92.refFlat.txt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=$RNA_HOME/refs/ref_ribosome.interval_list \; | sh
```

#### RSeQC (optional)

**Background:** RSeQC is a tool that can be used to generate QC reports for RNA-seq. For more information, see the [RSeQC tool homepage](https://rseqc.sourceforge.net/).

Files needed:

- Aligned BAM files
- An index file for each BAM file
- A transcript BED file (in BED12 format)

```bash
cd $RNA_HOME/refs/

# Convert GTF to genePred
gtfToGenePred chr22_with_ERCC92.gtf chr22_with_ERCC92.genePred

# Convert genePred to BED12
genePredToBed chr22_with_ERCC92.genePred chr22_with_ERCC92.bed12

# Create a version of the BED12 with the ERCC92 sequences removed (because they confuse the QC analysis that looks at how many reads are coding, UTR, intronic, etc.)
grep -v ERCC chr22_with_ERCC92.bed12 > chr22_without_ERCC92.bed12

cd $RNA_ALIGN_DIR
mkdir rseqc
geneBody_coverage.py -i chr22_only_bams/UHR_Rep1.bam,chr22_only_bams/UHR_Rep2.bam,chr22_only_bams/UHR_Rep3.bam -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 -o rseqc/UHR
geneBody_coverage.py -i chr22_only_bams/HBR_Rep1.bam,chr22_only_bams/HBR_Rep2.bam,chr22_only_bams/HBR_Rep3.bam -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 -o rseqc/HBR

# Calculate the inner distance of RNA-seq fragments.
#              RNA fragment
#  _________________||_________________
# |                                    |
# |                                    |
# ||||||||||------------------||||||||||
#   read_1    inner_distance    read_2
#
# fragment size = read_1 + inner_distance + read_2

find *Rep*.bam -exec echo inner_distance.py -i {} -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 -o rseqc/{} \; | sh

# Annotate exon-exon junctions observed in RNA-seq alignments compared to known exon-exon junctions
cd $RNA_ALIGN_DIR/chr22_only_bams
find *Rep*.bam -exec echo junction_annotation.py -i {} -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 -o ../rseqc/{} \; | sh

# Perform a saturation analysis using only exon-exon junction mapping reads
find *Rep*.bam -exec echo junction_saturation.py -i {} -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 -o ../rseqc/{} \; | sh

# Determine the distribution of reads with respect to the parts of transcripts they align to (e.g. 5' UTR, CDS, 3'UTR, intron, etc.)
find *Rep*.bam -exec echo read_distribution.py  -i {} -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 \> ../rseqc/{}.read_dist.txt \; | sh

# Calculate the RNA fragment sizes and produce statistics for each transcript
cd $RNA_ALIGN_DIR
find *Rep*.bam -exec echo RNA_fragment_size.py -i {} -r $RNA_HOME/refs/chr22_without_ERCC92.bed12 \> rseqc/{}.frag_size.txt \; | sh

# Summarizing mapping statistics of each BAM file
find *Rep*.bam -exec echo bam_stat.py -i {} \> rseqc/{}.bam_stat.txt \; | sh

rm -f log.txt
```

#### MultiQC

We will now use MultiQC to compile a QC report from all the QC tools above:

```bash
cd $RNA_ALIGN_DIR
multiqc ./
```

Open `rnabio-workspace/rnaseq/alignments/hisat2/multiqc_report.html` from your mounted storage. It should look something like this:

![MultiQC](https://rnabio.org/assets/module_2/multiqc.png)

#### View a pre-generated MultiQC report for full BAM files

Our data has been heavily down-sampled to make it run quickly. To see what a MultiQC report looks like for QC reports from non-downsampled BAM files, download this example:

```bash
mkdir $RNA_ALIGN_DIR/example_QC
cd $RNA_ALIGN_DIR/example_QC
wget https://genomedata.org/rnaseq-tutorial/multiqc_report.html
```

Then open `rnabio-workspace/rnaseq/alignments/hisat2/example_QC/multiqc_report.html` from your mounted storage. Below is a brief description of each of the samples included in that report:

| Name     |        Sample type         |
|----------|----------------------------|
| Sample 1 | Brain metastasis           |
| Sample 2 | Melanoma xenograft         |
| Sample 3 | Melanoma cell line         |
| Sample 4 | Melanoma                   |
| Sample 5 | Small Cell Lung Cancer FFPE|
| Sample 6 | Brain metastasis           |

---

## 8. Optional challenge: alignment on a new dataset

*Adapted from the RNAbio [Team Assignment - Alignment](https://rnabio.org/module-02-alignment/0002/07/01/Team_Assignment_Alignment/).*

This section is optional. The goal is to gain hands-on experience by working on recently published RNA-seq data, applying the concepts you have learned up through alignment. To complete it, you'll need to review the commands we ran in earlier sections. There are no step-by-step commands here!

**Background on the dataset:** We will be using subsets of the GSE136366 dataset (Roczniak-Ferguson A, Ferguson SM. Pleiotropic requirements for human TDP-43 in the regulation of cell and organelle homeostasis. Life Sci Alliance 2019 Oct;2(5). PMID: 31527135). This dataset consists of 6 RNA sequencing files of human cells that either express or lack the TDP-43 protein.

**Experimental details:**

- The libraries are prepared as paired end.
- The samples are sequenced on an Illumina HiSeq 2500.
- Each read is 63 bp long.
- The data are RF/fr-firststrand stranded (dUTP).
- The source dataset is located here: [GSE136366](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136366)
- 3 samples are from TDP-43 knockout HeLa cells, and 3 samples are cells in which a wildtype TDP-43 transgene was re-expressed.
- The reads have been subset to a single chromosome. Pick one of these datasets:
  - Dataset A: chr11
  - Dataset B: chr12
  - Dataset C: chr17
  - Dataset D: chr19
  - Dataset E: chr6
- The files are named based on their SRR IDs, and obey the following key:
  - SRR10045016 = KO sample 1
  - SRR10045017 = KO sample 2
  - SRR10045018 = KO sample 3
  - SRR10045019 = Rescue sample 1
  - SRR10045020 = Rescue sample 2
  - SRR10045021 = Rescue sample 3

> [!NOTE]
> Building a HISAT2 index for a whole chromosome takes more memory than our small chr22 example did. If `hisat2-build` gets killed, exit your container and re-launch it with `--mem=16G`.

### Part I - Obtaining the dataset and reference files

**Goals:**

- Obtain the files necessary for data processing
- Review reference and annotation file formats
- Review sequence FASTQ format

Download the dataset you picked:

```bash
cd $RNA_HOME/
mkdir -p team_exercise/untrimmed
cd team_exercise/untrimmed

# Fill in the "XX" below with your dataset's letter (A, B, C...)
wget -c https://genomedata.org/seq-tec-workshop/read_data/rna_alignment-de_exercise/dataset_XX/dataset.tar.gz
tar -xzvf dataset.tar.gz
```

Then create a separate directory and download the corresponding reference files needed for alignment and expression analysis. Don't forget to modify the commands below to **use your dataset's chromosome**:

```bash
mkdir -p $RNA_HOME/team_exercise/references
cd $RNA_HOME/team_exercise/references

## Adapter trimming
wget -c https://genomedata.org/seq-tec-workshop/references/RNA/illumina_multiplex.fa

## Reference fasta corresponding to your assigned chromosome (e.g. chr6)
wget -c https://genomedata.org/seq-tec-workshop/references/RNA/chrXX.fa

## Obtain annotated reference gtf file corresponding to your assigned chromosome (e.g. chr6)
wget -c https://genomedata.org/seq-tec-workshop/references/RNA/chrXX_Homo_sapiens.GRCh38.95.gtf
```

### Part II - Data preprocessing (QC and trimming)

**Goals:**

- Perform adapter trimming on your data, and also pre-trim 5 bases from the **end (right)** of reads
- Perform QC on your data with `fastqc` and `multiqc` before and after trimming

**Q1.** What is the average percentage of reads that are trimmed?

**Q2.** How do you expect the sequence length distribution to look prior to and after trimming? Is your answer confirmed by the MultiQC report results?

**Q3.** Are there any metrics where the sample(s) failed?

### Part III - Alignment

**Goals:**

- Create HISAT2 index files **for your chromosome**
- Review HISAT2 alignment options
- Perform alignments (remember to use `-p 2`!)
- Obtain an alignment summary
- Sort and convert your alignments into compressed BAM format

*A useful option to add to the end of your commands is `2>`, which redirects the STDERR from any command into a specific file. This can be used to redirect your STDERR into a summary file, like this: `my_alignment_command 2> alignment_metrics.txt`. The advantage of this is being able to view the alignment metrics later on.*

**Q4.** What were the percentages of reads that aligned to the reference for each sample?

**Q5.** By compressing your SAM format to BAM, approximately how much space is saved (fold change in size)?

### Part IV - Post-alignment QC

**Goals:**

- Perform post-alignment QC analysis using `fastqc` and `multiqc`
- Merge BAM files (one for each condition)
- Index the BAM files

**Q6.** How does the information from your post-alignment QC report differ from pre-alignment QC?

---

## 9. Cleaning up

When you're done for the day, type `exit` to leave the container and release your interactive job. (Remember, the job holds its CPUs and memory as long as it's open, even if you're not using it.)

**Don't delete anything in `rnabio-workspace`!** We'll keep building on these reference files, indexes, and alignments in the coming weeks.

---

## 10. Cheat sheet

```bash
# --- on the login node: start the RNA-seq toolkit container ---
srun -A compute2-workshop -p workshop -c 2 --mem=8G --pty \
  --container-image=griffithlab/rnaseq-toolkit:1.1.1 \
  --container-mounts="/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/$USER/rnabio-workspace:/workspace" \
  --container-workdir=/workspace \
  /bin/bash

# --- inside the container: load the course variables (every session!) ---
source /workspace/rnabio_env.sh
echo $RNA_HOME                                  # should print /workspace/rnaseq

# --- view HTML reports: mount storage on your laptop (Mac) ---
open smb://storage1.ris.wustl.edu/c.a.miller/Active/bfx-workshop-scratch/<washukey>
```
