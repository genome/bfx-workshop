# Aligning Sequence Data: Indexing, Alignment, and BAM files on the Cluster

Last week we downloaded some sequence data, unpacked it, and ran QC on it using a container. This week we'll take the next step and align reads to a reference genome. By the end you should be able to:

- Index a reference genome with `samtools` and `bwa`
- Align paired-end FASTQ files with `bwa mem`, and attach read group information
- Stream the whole process (FASTQ → alignment → BAM) in a single piped command
- Look inside aligned BAM files, and decode SAM flags
- Merge read groups into a single BAM, and understand what the header is telling you

> [!NOTE]
> As always, we strongly suggest typing the commands out yourself rather than copying and pasting. It really does help it stick.

Throughout this lesson, questions are followed by expandable **Hint** and **Solution** boxes. Click through to get help or check your work, but take a stab at doing it yourself first!

---

## 1. The tools

We'll be using two tools this week, both of which are pre-installed in a commonly used demonstration image described in the O'Reilly book [Genomics in the Cloud](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/):

```
broadinstitute/genomes-in-the-cloud:2.3.1-1512499786
```

(Since there's no registry at the front of that name, it comes from [Docker Hub](https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud).)

- **[Samtools](http://www.htslib.org/)** is a suite of programs for interacting with high-throughput sequencing data. Within the container, it is installed at `/usr/local/bin/samtools`, which is on your `PATH`, so you can just type `samtools`.

- **[BWA](https://github.com/lh3/bwa)** is a software package for mapping DNA sequences against a large reference genome, such as the human genome. Within the container, it is installed at `/usr/gitc/bwa`. That directory is *not* on your `PATH`, so we'll always type out the full path.

---

## 2. Setup

### Download the inputs

We are using a toy data set based on the HCC1395 blood normal cell line. The sequence reads and genome reference are a subset targeting chr6 (the HLA genes region) and chr17 (genes TP53 and BRCA1).

Log in to Compute2, then move into your workshop directory and make a folder for this week, with a few subdirectories to keep things organized:

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu
cd ~/workshop
mkdir week04
cd week04
mkdir -p ref
mkdir -p unaligned/normal
mkdir -p aligned/normal
```

(The `-p` flag tells `mkdir` to create any missing parent directories, and not to complain if the directory already exists.)

Next, download the input files: a reference FASTA, and two pairs of FASTQ files (feel free to copy and paste these commands!) Just like last week, downloading these smallish files is fine to do on the login node:

```bash
curl -O https://storage.googleapis.com/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes.fa
curl -O https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/alignment_exercise/2895499331_1.fastq.gz
curl -O https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/alignment_exercise/2895499331_2.fastq.gz
curl -O https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/alignment_exercise/2895499399_1.fastq.gz
curl -O https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/alignment_exercise/2895499399_2.fastq.gz
```

Then move them into the directories we created:

```bash
mv hla_and_brca_genes.fa ref
mv *.fastq.gz unaligned/normal
```

Use `ls` liberally to check that everything landed where you expected:

```bash
ls -lh ref unaligned/normal
```

The reference FASTA is about 246 MB, and each of the FASTQ files is under 200 KB.

### Launch an interactive job in the container

Now we need our tools. Just like last week, they're not installed on the head node, so we'll ask Slurm for an interactive job running inside a docker container, with our storage directory mounted so that we can see our data:

```bash
srun -A compute2-workshop -p workshop -c 1 --mem=8G --pty \
  --container-image=broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch:/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch" \
  /bin/bash
```

This is the same command we used for FastQC, with two changes: a different `--container-image`, and a little more memory (`--mem=8G`), since BWA is hungrier than FastQC.

> [!NOTE]
> This image is about 2 GB, so the first launch may take a few minutes while it downloads. Be patient!

Once you have a prompt, move back to this week's directory:

```bash
cd ~/workshop/week04
```

> [!IMPORTANT]
> Unless otherwise noted, **all of the remaining commands in this lesson are run inside this container**. If you get disconnected or `exit`, just re-run the `srun` command above and `cd` back to `~/workshop/week04`. Everything you've created will still be there, because it was written to your storage directory.

---

## 3. Indexing

Index files are used all over genomics to provide fast, random access to specific records, locations, or content within much larger files. For example: pulling the sequence of one gene out of an entire genome, finding the reads that overlap a position out of billions of reads, or finding one variant among millions of records.

### Samtools faidx

Using `samtools faidx`, we create an index (`.fai`) of the reference FASTA (`.fa`). Remember, this kind of *file index* is like a table of contents.
This lets samtools (and many other tools) quickly jump to the nucleotides at any position in the genome, without reading through the whole file.

```bash
samtools faidx ref/hla_and_brca_genes.fa
```

We should now see both the FASTA and the FAI files:

```bash
ls ref
```

### Question 1

Take a look at the first twenty lines of the FASTA. What do you see? What do you NOT see?

```bash
head -n 20 ref/hla_and_brca_genes.fa
```

<details>
<summary>Solution</summary>

You'll see one header line (`>chr6`), followed by 19 lines of nothing but `N` characters. You do NOT see any actual A, C, G, or T bases!

`N` means "any/unknown base." Two things are going on here:

- The very beginning of every chromosome is a telomere, which is made of highly repetitive sequence that can't be assembled well, so even the full human reference starts each chromosome with a big run of `N`s.
- This toy reference has also been **masked**: every base outside of our target genes has been replaced with `N`.  This makes alignment way faster for example purposes. 

The `.fai` shows that chr6 and chr17 are their full lengths (170,805,979 and 83,257,441 bp), but only about 247,000 of those ~254 million bases are real sequence. Keeping the full chromosome lengths means that coordinates in this reference match the real human genome (GRCh38), which will come in handy when we view our alignments in IGV.

</details>

### Question 2

Now use `samtools faidx` again, but this time pass it a chromosome name and start and end positions:

```bash
samtools faidx ref/hla_and_brca_genes.fa chr17:43044295-43170245
```

What do you see in the output? Does this look more familiar than the `head` command run earlier?

<details>
<summary>Solution</summary>

This time we get real sequence (A, C, G, T), in FASTA format, with a header line named after the region we asked for (`>chr17:43044295-43045295`).

That region is the beginning of the  **BRCA1** gene, one of the regions that wasn't masked. Notice that samtools returned it instantly, even though it's ~174 million bytes into the file. That's the index making things fast - rather than reading through the big fasta file from the top, samtools used the `.fai` to calculate exactly where to jump.

</details>

### BWA index

BWA needs its own, much more elaborate index of the reference. It actually consists of several files and data structures that all work together (and are all required) for the alignment commands (`mem`, `aln`, `sampe`, etc). These are the same concept a the hash we used in our alignment example in the lecture. We can use `bwa index` to create this index for our genome:

```bash
/usr/gitc/bwa index ref/hla_and_brca_genes.fa
```

> [!NOTE]
> This step takes a few minutes. Even though most of our reference is `N`s, BWA still has to process all ~254 million positions. Indexing the full human genome takes an hour or more, which is why you'll usually want to find a pre-built index rather than making your own!

When it finishes, see what it created:

```bash
ls ref
```

You'll see five new files (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`). When you run `bwa mem`, you only give it the path to the `.fa` file, and it finds the others automatically, so they always need to stay together in the same directory.

---

## 4. Alignment

We have two read groups from the same sample (the HCC1395 Normal). A **read group** is a set of reads that were all sequenced together, typically from one library on one lane of a flowcell. Each read group comes as a pair of FASTQ files, one for read 1 and one for read 2. We'll align them in two slightly different ways:

1. Step by step: align the FASTQ files to make a SAM file, then convert that to BAM
2. All at once: stream the data straight from FASTQ to an aligned BAM, using a pipe

### Looking at the FASTQ files

```bash
ls unaligned/normal
```

### Question 3

The FASTQ files are gzip-compressed. What happens if you look at one directly with `head`?

```bash
head unaligned/normal/2895499331_1.fastq.gz
```

Does this look like a FASTQ file? How could you look at the contents without decompressing the file?

<details>
<summary>Solution</summary>

You'll see a screen full of garbled characters (and it may even mess up your terminal - type `reset` and hit enter if it does). That's because `head` is showing the raw **compressed binary** data, which isn't meant for humans to read.

Last week, we used `gunzip` to decompress the whole file, then used `gzip` to compress it again. For a real data set, that wastes a lot of time and disk space. Instead, use `gunzip -c` (`-c` means "unzip the file without modifying the original") to decompress the data and pipe the result to `head`:

```bash
gunzip -c unaligned/normal/2895499331_1.fastq.gz | head -n 8
```

That looks like a proper FASTQ file: records of four lines each, with a read name (ending in `/1`, since this is the read 1 file), the sequence, a `+` separator, and quality scores. Notice that these reads are 151 bp long.

</details>

### Question 4

How many read pairs are in each read group?

<details>
<summary>Hint</summary>

Remember from last week: every FASTQ record is four lines. You'll need to combine `gunzip -c` with the line-counting command we used last week.

</details>

<details>
<summary>Solution</summary>

```bash
gunzip -c unaligned/normal/2895499331_1.fastq.gz | wc -l
gunzip -c unaligned/normal/2895499399_1.fastq.gz | wc -l
```

That's 7,088 lines ÷ 4 = **1,772 read pairs** for 2895499331, and 6,544 lines ÷ 4 = **1,636 read pairs** for 2895499399. (Check the `_2` files too - they should match!) That's a tiny amount of data, which is why the alignments below will only take a few seconds.

</details>

### Align step by step

Now let's align! Give `bwa mem` the reference, followed by the read 1 and read 2 FASTQ files, and pipe the output to a file. 

```bash
/usr/gitc/bwa mem ref/hla_and_brca_genes.fa unaligned/normal/2895499331_1.fastq.gz unaligned/normal/2895499331_2.fastq.gz >aligned/normal/2895499331.sam
```
`bwa mem` writes its alignments to **STDOUT** in **SAM** format, and writes its progress messages to **STDERR**, which still showed up on the screen. Now let's look at the results:

```bash
less aligned/normal/2895499331.sam
```

> [!Note]
> Reminder - to exit `less`, type `q`, for "quit".


The first lines of our SAM file are a header. `@SQ` lines describe the reference sequences.  (Remember, we're using a toy reference with two chromosomes, so there are only two lines here, instead of a list of all the human chromosomes).  Next is a `@PG` line recording the `bwa mem` command that made the file. In the body of the file , we see one line per read with its name, the chromosome and position it aligned to, a mapping quality, a CIGAR string, and the sequence and qualities from the FASTQ.

It's human-readable, but SAM is plain text and not compressed, so it won't scale well. Compare the size of the SAM file to the two gzipped FASTQ files it came from:

```bash
ls -lh aligned/normal/2895499331.sam unaligned/normal/2895499331_*.fastq.gz
```

The SAM is about 1.5 MB, compared to about 350 KB for the two FASTQs combined, even though it holds the same reads! For a whole genome, that difference is hundreds of gigabytes. That's why we almost always store alignments in the compressed binary **BAM** (or even more compressed **CRAM**) format instead.

Notice one other thing: this SAM header has no `@RG` (read group) line! FASTQ files don't carry any information about which sample, library, or lane the reads came from, so there was nothing for `bwa mem` to put there. We'll fix that next.

### Adding read groups and converting to BAM

Let's try aligning this read group again, this time outputting to compressed BAM format. To do so, we'll pipe the SAM output directly to `samtools view`, which can convert between formats.

We'll also give `bwa mem` a `-R` argument containing the read group information, which it will put in the header and attach to each read. This information usually comes from whoever did the sequencing (your sequencing core, or the metadata in a public repository). Ours came from the sequencing center, and looks like this:

| Read group ID | Platform | Platform unit (flowcell-barcode.lane) | Library | Sample |
|---|---|---|---|---|
| 2895499331 | ILLUMINA | H7HY2CCXX-TGACCACG.3 | H\_NJ-HCC1395-HCC1395_BL-lg21-lib1 | H\_NJ-HCC1395-HCC1395_BL |
| 2895499399 | ILLUMINA | H7HY2CCXX-TGACCACG.4 | H\_NJ-HCC1395-HCC1395_BL-lg21-lib1 | H\_NJ-HCC1395-HCC1395_BL |

```bash
/usr/gitc/bwa mem -R "@RG\tID:2895499331\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.3\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" ref/hla_and_brca_genes.fa unaligned/normal/2895499331_1.fastq.gz unaligned/normal/2895499331_2.fastq.gz | samtools view -o aligned/normal/2895499331.bam -
```

That's a long one, and another one where it's okay to copy/paste, but make sure you look at all the new pieces here:

| Piece | Meaning |
|---|---|
| `-R "@RG\t..."` | the read group header line to add. `\t` stands for a tab character |
| `|` | pipe the SAM output of `bwa mem` into `samtools view` |
| `-o aligned/...bam` | the output file |
| `-` | read the input from STDIN (i.e. from the pipe) rather than from a file |

This file is compressed, like our fastq file was eariler. That means we can't just use `less` or `cat` to look at the file directly, we have to uncompress it first.  We use `samtools` to do this:

Take a look at the header of the BAM file we just created (`-H` means "header"):

```bash
samtools view -H aligned/normal/2895499331.bam
```

Now look at the alignments stored in the body of the BAM file. There are going to be thousands of these, so let's pipe it to `less` again:

```bash
samtools view aligned/normal/2895499331.bam | less
```

### Question 5

The second column of each record is the **SAM flag**, a single number that packs a bunch of yes/no facts about the read. What are the flags of the first two reads, and what do they mean? Use this web app to decode them: https://broadinstitute.github.io/picard/explain-flags.html

<details>
<summary>Hint</summary>

To see just the first two reads, use `head -n 2`. To see just the first few columns, try `cut -f 1-4`.

```bash
samtools view aligned/normal/2895499331.bam | head -n 2 | cut -f 1-4
```

</details>

<details>
<summary>Solution</summary>

The flags are **83** and **163**:

| Flag | Decoded |
|---|---|
| 83 = 1 + 2 + 16 + 64 | read paired, read mapped in proper pair, read on the **reverse** strand, **first** in pair |
| 163 = 1 + 2 + 32 + 128 | read paired, read mapped in proper pair, **mate** on the reverse strand, **second** in pair |

So these are the two mates of one read pair (notice they have the same name in the first column). Both aligned to chr17, about 50 bp apart (columns 3 and 4), facing each other: one on the forward strand and one on the reverse strand. That's exactly what we expect from a normal paired-end fragment, so `bwa mem` marked them as a "proper pair."

Try decoding a few more. Can you find reads that aren't in a proper pair, or that didn't map at all?

</details>

### Align in one step

Now let's align our second read group, 2895499399. This time we'll go straight from FASTQ to BAM with a single piped command, and skip writing the SAM file.

### Question 6

Write the command to align the second read group, 2895499399, directly to a BAM file named `aligned/normal/2895499399.bam`. What needs to change from the command we used for the first read group?

<details>
<summary>Hint</summary>

Use the up arrow to find the `bwa mem -R ... | samtools view ...` command above, then edit it to use this new data instead. Look carefully at the read group table: which values are different for the second read group? And don't forget the input and output file names. 

</details>

<details>
<summary>Solution</summary>

```bash
/usr/gitc/bwa mem -R "@RG\tID:2895499399\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.4\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" ref/hla_and_brca_genes.fa unaligned/normal/2895499399_1.fastq.gz unaligned/normal/2895499399_2.fastq.gz | samtools view -o aligned/normal/2895499399.bam -
```

A few things change:

- the read group `ID`
- the lane at the end of the `PU` (`.3` becomes `.4`)
- the two input FASTQ files
- the output BAM file

The library (`LB`) and sample (`SM`) stay the same, because these reads came from the same library of the same sample, just sequenced on a different lane. Getting these right matters: downstream tools use `SM` to decide which reads belong to the same sample, and `LB` to find PCR duplicates.

</details>

### Question 7

Which files do you now have for read group 2895499399? How does that compare to read group 2895499331?

<details>
<summary>Solution</summary>

```bash
ls -lh aligned/normal
```

For **2895499331**, we have a SAM file and a BAM file. For **2895499399**, we have only the BAM file.

The piped approach streamed the alignments from `bwa mem` straight into `samtools view` in memory, so no intermediate SAM file was ever written to disk. For a whole genome, that SAM file could be hundreds of gigabytes, and writing it out and then reading it back in takes a lot of time. It also means there's nothing left over to clean up afterwards.

</details>


Since we don't need it anymore, go ahead and delete the SAM file:

```bash
rm aligned/normal/2895499331.sam
```

### Merge and sort alignments

Now that we have two read group BAMs, we want to combine them into a single BAM with all of the alignments for the Normal HCC1395 sample. At the same time, we'll **sort** the alignments by position. Right now, the reads are in the order they came off the sequencer, which means reads from the same part of the genome are scattered all through the file. Most downstream tools (including IGV, which we'll use for homework) need them ordered by chromosome and position instead.

We can do both in one step, by piping `samtools merge` into `samtools sort`:

```bash
samtools merge - aligned/normal/2895499331.bam aligned/normal/2895499399.bam | samtools sort -o aligned/normal.bam -
```

Just like with `samtools view` earlier, the `-` tells `samtools merge` to write its output to STDOUT (into the pipe), and tells `samtools sort` to read its input from STDIN (from the pipe).  

> [!NOTE]
> You'll see a warning that says `No @HD tag found`. That's harmless: it's just telling us that our input BAMs don't declare how they're sorted (because they aren't!). `samtools sort` takes care of that.

As always, use `ls` to make sure the file showed up as expected.

Let's take a closer look at this BAM header - how is it different from the previous BAM headers we have looked at?

```bash
samtools view -H aligned/normal.bam
```

- **`@RG`** - We now have two read group tags. Every read in the body of the file has an `RG:Z:` tag pointing back to one of these IDs, so we know where each read came from.
- **`@PG`** - we now also have several records documenting which program (and which version, and the exact command line) was used to produce the data. It's a built-in record of how the file was made!

There's also a new **`@HD`** line at the top with `SO:coordinate`. That was added by `samtools sort`, and it records that the reads are now ordered by chromosome and position, rather than in the order they came off the sequencer. That will be important for the homework!

</details>

### Index and view the data

Bams are most useful when they have an index file, which is that Table of Contents that lets us jump straight to various portions of the genome.  Create that index for the `normal.bam` file:

```
samtools index aligned/normal.bam
ls aligned
```

This creates `aligned/normal.bam.bai`.

Note that a BAM file must be **sorted by coordinate** before it can be indexed. Our merged BAM is, because we ran it through `samtools sort`. If you try to index one of the unsorted read group BAMs, like `aligned/normal/2895499331.bam`, samtools will give you an error.


## 5 Viewing the data with IGV

Let's get to the payoff and see what we've created. To do this, we'll view the indexed `normal.bam` file in [IGV](https://igv.org/doc/desktop/), on your laptop.

There are multiple ways to do this - we could use `scp` to copy the data from the cluster to our laptop, but real bam files are huge, and we don't want to fill up our hard drive.  Instead, we're going to "mount" the cluster's network storage drive to our laptop, so you can browse it just like files on your own machine.  

This procedure is slightly different, depending on whether you're on Mac or Windows, but the end result is the same.  Use the RIS documentation to see how you can connect:

#### MacOS
- [Connecting to Storage from MacOS](https://washu.atlassian.net/wiki/spaces/RUD/pages/1795784747/Connecting+to+Storage+from+MacOS)

If you prefer to use Finder like they describe, that's fine, but I think it's sometime easier to use a Terminal and just type:

```
open smb://storage1.ris.wustl.edu/c.a.miller/Active/bfx-workshop-scratch/<wustl_key>/
```
(as always, replacing \<wustl_key> with your actual username)

Login with your wustl credentials when prompted. 

The resulting shortcut folder will show up under /Volumes/<wustl_key> in Finder

#### Windows

- [Connecting to Storage from+Windows](https://washu.atlassian.net/wiki/spaces/RUD/pages/1795588135/Connecting+to+Storage+from+Windows)

You'll want to use 
`\\storage1.ris.wustl.edu\c.a.miller\Active\bfx-workshop-scratch\<wustl_key>`
(as always, replacing \<wustl_key> with your actual username)


### Fire up IGV

1. If you don't have it installed, grab it from [https://igv.org/doc/desktop/](https://igv.org/doc/desktop/)

2. Choose the **Human (hg38)** genome from the dropdown at the top left. Our toy reference uses full-length chr6 and chr17 with GRCh38 coordinates, so the alignments line up with the standard genome.
3. Choose **File > Load from File...** and select `normal.bam`.
4. Type `BRCA1` into the search box and hit Go.

Zoom in until you see the reads.  Notice how they cluster over exons - this was clearly exome sequencing.  

We'll dive much more deeply into IGV in future sessions.

### Cleaning up

When you're done, type `exit` in the container to release your interactive job.

---

## 5. Cheat sheet

```bash
# --- interactive shell in the alignment container ---
srun -A compute2-workshop -p workshop -c 1 --mem=8G --pty \
  --container-image=broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch:/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch" \
  /bin/bash

# --- indexing ---
samtools faidx ref.fa                          # index a FASTA (.fai)
samtools faidx ref.fa chr17:43044295-43170245  # pull out one region
/usr/gitc/bwa index ref.fa                     # build the BWA index (slow!)

# --- looking at FASTQ/SAM/BAM files ---
gunzip -c reads.fastq.gz | head -n 8           # peek at a gzipped FASTQ
samtools view -H file.bam                      # header only
samtools view file.bam | head                  # first few records
samtools view -o out.bam -                     # convert SAM on STDIN to BAM

# --- alignment ---
/usr/gitc/bwa mem -R "@RG\tID:..." ref.fa reads_1.fastq.gz reads_2.fastq.gz | samtools view -o out.bam -

# --- merge and sort ---
samtools merge - a.bam b.bam | samtools sort -o merged.bam -
samtools index merged.bam                      # index a sorted BAM (.bai)

# --- piping safely ---
set -o pipefail
```

---
