# Tool Introduction

## Docker

As discussed in Week 1, we will be using Docker throughout this workshop.
Hopefully everyone has Docker installed on their local environment, if not please see [Week1](../week_01).
We are pulling a commonly used demonstration image described in the O'Reilly book [Genomics in the Cloud](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/).
"Pulling" the image means that Docker is downloading the binary image that includes all of the necessary software tools pre-installed.

```
docker pull broadinstitute/genomes-in-the-cloud:2.3.1-1512499786
```

## Samtools

[Samtools](http://www.htslib.org/) is a suite of programs for interacting with high-throughput sequencing data. Within the Docker container, it is installed at `/usr/local/bin/samtools`.

## BWA

[BWA] is a software package for mapping DNA sequences against a large reference genome, such as the human genome. Within the Docker container, it is installed at `/usr/gitc/bwa`.

## Picard

[Picard](https://broadinstitute.github.io/picard/) is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF. It is run via Java, which is also available in the provided image; it can be launched using the following command: `java -Xms2G -jar /usr/gitc/picard.jar`

# Setup

## Launch Docker
Since the Docker image has all of the tools we need, we will do most of our work within a container. Launch an interactive shell using this command: 
```
docker run -it -v $PWD:/data broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 bash
```

On some devices, namely newer Macs, you may see a warning message such as
> WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested

This is just a warning, not an error. However, the message can be annoying, so you may want to add `--platform linux/amd64` to your `docker run` command: 
```
docker run -it -v $PWD:/data --platform linux/amd64 broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 bash
```

Note: the `-v $PWD:/data` portion of the command tells Docker to mount your current working directory within the container at the path `/data`. This allows us to make files within the container at this path and still access them after we close the container. By default, the filesystem within the container is completely separate from the filesystem on your device, and any files made within the container will be lost after exiting. So, make sure anything you want to save is in the `/data` directory before you `exit` the container!

If you docker command is successful, you'll notice your shell prompt (the portion to the left of your cursor) will change to something like `root@b42034319df5:/usr/gitc#`. This is telling you that within the image, you're running as a user named "root", with a machine name "b42034319df5", located at the path `/usr/gitc`. 

Let's switch to the working directory we mounted from our device filesystem:
```
cd /data
```

Now, let's set up a few directories for our data:

```
mkdir -p ref
mkdir -p unaligned/normal
mkdir -p aligned/normal
```

## Download Inputs

We are using a toy example data set based on the HCC1395 blood normal cell line. The sequence reads and genome reference are a subset targeting chr6, genes HLA-A and HLA-B-C, and chr17, genes TP53 and BRCA1. Download the following files:

- [FASTA](https://storage.googleapis.com/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes.fa)
- [Normal reads lane 3](https://storage.googleapis.com/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499331.bam)
- [Normal reads lane 4](https://storage.googleapis.com/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499399.bam)

Example download command:
```
wget <url>
```

Move the files to the directories we created above:

```
mv hla_and_brca_genes.fa ref
mv 2895499331.bam unaligned/normal
mv 2895499399.bam unaligned/normal
```

# Indexing

Index files of various formats and data structures are used commonly in genomics to provide random access to specific records, locations, or content within much larger domains and coordinate spaces. Ex. entire genome nucleotide sequences, billions of sequence reads, millions of variant records, etc.

## Samtools

Using `samtools faidx`, we create an index (.fai) of the reference FASTA (.fa) which provides random access to the nucleotides at specific positions within the complete genome reference. The FAIDX index is used both by the samtools as well as other toolkits, algorithms, and libraries that require random access to specific coordinates in a timely manner.

```
/usr/local/bin/samtools faidx ref/hla_and_brca_genes.fa
```

We should now see the FASTA and FAIDX format files.

```
ls ref
```

Take a look at the first twenty lines of this FASTA. What do you see? What do you NOT see?

```
head -n 20 ref/hla_and_brca_genes.fa
```

Using `samtools faidx` once again, let's pass the name of the chromosome and some example start and end positions.

```
/usr/local/bin/samtools faidx ref/hla_and_brca_genes.fa chr17:43044295-43170245
```

Now, what do you see in the command-line output? Does this look more familiar than the head command run earlier?

## BWA

BWA index will apply the correct index algorithm depending on the size of the reference provided. The index actually consists of several files and data structures that all work together (and are required) to perform alignment with subsequent alignment commands, ex. mem, aln, sampe, etc.

```
/usr/gitc/bwa index ref/hla_and_brca_genes.fa

ls ref
```

# Alignment

Let us demonstrate the alignment of two Normal HCC1395 read groups (from the same sample) using both FASTQ format files (created from our unaligned BAMs) as well as a performance improvement to align sequence reads directly from unaligned BAM to an aligned BAM in one streaming process (using pipes).

## Align from FASTQ

Let's walk through each step of the alignment process including the conversion from unaligned BAM to compressed FASTQ.

```
ls unaligned/normal
```

The `-H` flag is provided to the `samtools view` subcommand so we can interogate the unaligned BAM header.

```
/usr/local/bin/samtools view -H unaligned/normal/2895499331.bam
```

Without the `-H` flag, the default behavior for the `samtools view` subcommand is to print each SAM record.

```
/usr/local/bin/samtools view unaligned/normal/2895499331.bam
```

Let's look at the SAM flag in this web app to decode what values are set for the first two reads: https://broadinstitute.github.io/picard/explain-flags.html

Now, we are going to use the `SamToFastq` command from the Picard toolkit to convert our unaligned BAM file to a pair of FASTQ files (the default input to BWA). The OUTPUT_PER_RG option is used to output one pair of FASTQ (both read1 and read2 for paired-end) to the output directory of our choosing. The output compressed FASTQ files will be named based on the read group ID.

```
java -Xms2G -jar /usr/gitc/picard.jar SamToFastq INPUT=unaligned/normal/2895499331.bam OUTPUT_PER_RG=true COMPRESS_OUTPUTS_PER_RG=true RG_TAG=ID OUTPUT_DIR=unaligned/normal

ls unaligned/normal
```

Since the FASTQ files are compressed binary files, we can not view them as-is. As an example, we will decompress one of the read group files with gunzip. However, for an actual data set, keep the compressed binary format as-is to save space and rely on tools such as `zcat`, `zgrep`, `zless`, `zdiff`, etc. to perform basic manipulation (when necessary).

```
head unaligned/normal/2895499331_1.fastq.gz
```

What are the characters above? Does this look like a FASTQ format file?

```
gunzip -f unaligned/normal/2895499331_1.fastq.gz

ls unaligned/normal

head unaligned/normal/2895499331_1.fastq
```

How does the decompressed FASTQ look now? This format is much easier for humans to read, but less space efficient. Let's binary compress the file again.

```
gzip unaligned/normal/2895499331_1.fastq
```

The following steps are a series of BWA mem and Samtools view commands that:
- align the input FASTQ files to the FASTA reference
- redirect the output to a SAM format file
- convert SAM to BAM (compressed binary) file formats

```
ls unaligned/normal

/usr/gitc/bwa mem ref/hla_and_brca_genes.fa unaligned/normal/2895499331_1.fastq.gz unaligned/normal/2895499331_2.fastq.gz
```

What does the output of the above command look like? Is it a familiar file format we have discussed? What will we do with this output format? Would this scale for a whole genome?

Maybe we should re-run this command and output to a proper file. Let's take a look at `samtools view`, a basic command often used to read and write alignment formats.

```
/usr/gitc/bwa mem ref/hla_and_brca_genes.fa unaligned/normal/2895499331_1.fastq.gz unaligned/normal/2895499331_2.fastq.gz > aligned/normal/2895499331.sam

head aligned/normal/2895499331.sam
```

Does the above output look familiar? What format is it? Will this scale?

Let's look at an example alignment with a binary compressed format output called BAM. To do so, we will pipe the SAM format directly to the samtools view command mentioned earlier. Note the longer command in one set of quotes.

```
/usr/gitc/bwa mem -R "@RG\tID:2895499331\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.3\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" ref/hla_and_brca_genes.fa unaligned/normal/2895499331_1.fastq.gz unaligned/normal/2895499331_2.fastq.gz | /usr/local/bin/samtools view -1 -o aligned/normal/2895499331.bam -
```

Take a look at the header of the BAM format file we just created.

```
/usr/local/bin/samtools view -H aligned/normal/2895499331.bam
```

Now look at the alignments stored in the body of the BAM file.

```
/usr/local/bin/samtools view aligned/normal/2895499331.bam | head
```

## Align from Unaligned BAM

Let's walk through the alignment process again for another read group. This time we will pipe multiple commands together in one single streaming process from unaligned BAM to aligned BAM.

```
/usr/local/bin/samtools view -H unaligned/normal/2895499399.bam
```

When running commands piped together, best practice is to set these bash arguments to generate appropriate error messages if something goes wrong such as an intermediate step is interupted or the streaming data becomes corrupt/incomplete.

```
set -o pipefail
set -o errexit

java -Xms2G -jar /usr/gitc/picard.jar SamToFastq INPUT=unaligned/normal/2895499399.bam FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true | /usr/gitc/bwa mem -R "@RG\tID:2895499399\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.4\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" -p ref/hla_and_brca_genes.fa /dev/stdin | /usr/local/bin/samtools view -1 -o aligned/normal/2895499399.bam -

ls aligned/normal
```

Which file formats do you now see for read group 2895499399? As compared to 2895499331 the read group aligned earlier?

Take a closer look at the BAM output from the single, piped command:

```
/usr/local/bin/samtools view -H aligned/normal/2895499399.bam 

/usr/local/bin/samtools view aligned/normal/2895499399.bam | head 
```

## Merge Alignments

Now that we have two read group BAMs, we want to `MergeSamFiles` (which works on BAM inputs) to produce a single BAM with all of our alignments for the Normal HCC1395 sample.

```
java -Xms2G -jar /usr/gitc/picard.jar MergeSamFiles OUTPUT=aligned/normal.bam INPUT=aligned/normal/2895499331.bam INPUT=aligned/normal/2895499399.bam

ls aligned

/usr/local/bin/samtools view -H aligned/normal.bam 
```

Do we really need to see another BAM header? Actually, take a closer look. How is this BAM header different from previous BAM headers we have looked at?

@RG and @PG tags - What do those stand for? How many tag lines do you see?

# Homework

- Index the normal.bam file. HINT: samtools index OR igvtools
- View the indexed normal.bam file with IGV (outside of the docker container) HINT: Search for BRCA1.
- Make a list of questions and/or observations about the alignments to discuss next week.
- Are there other post-alignment processing steps we've missed? Bring suggestions for next week.

