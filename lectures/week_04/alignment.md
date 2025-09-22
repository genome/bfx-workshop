# Setup

Let's setup a few directories and explore a few software tools that we will use for this week's tutorial.

```
echo $PWD
mkdir -p $PWD/ref
mkdir -p $PWD/unaligned/normal
mkdir -p $PWD/aligned/normal
```

## Docker

As discussed in Week 1, we will be using Docker throughout this workshop.
Hopefully everyone has Docker installed on their local environment, if not please see [Week1](../week_01).
We are pulling a commonly used demonstration image described in the O'Reilly book [Genomics in the Cloud](https://www.oreilly.com/library/view/genomics-in-the/9781491975183/).
"Pulling" the image means that Docker is downloading the binary image that includes all of the necessary software tools pre-installed.

```
docker pull broadinstitute/genomes-in-the-cloud:2.3.1-1512499786
```

## Samtools

[Samtools](http://www.htslib.org/) is a suite of programs for interacting with high-throughput sequencing data.

If you see a warning message like
> WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested

you may want to add `--platform linux/amd64` to each `docker run` command below

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools
```

## BWA

[BWA] is a software package for mapping DNA sequences against a large reference genome, such as the human genome.

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/gitc/bwa
```

## Picard

[Picard](https://broadinstitute.github.io/picard/) is a set of command line tools for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 java -Xms2G -jar /usr/gitc/picard.jar
```

## Inputs

We are using a toy example data set based on the HCC1395 blood normal cell line. The sequence reads and genome reference are a subset targeting chr6, genes HLA-A and HLA-B-C, and chr17, genes TP53 and BRCA1. Download the following files (`wget` recommended):

- [FASTA](https://storage.googleapis.com/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes.fa)
- [Normal reads lane 3](https://storage.googleapis.com/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499331.bam)
- [Normal reads lane 4](https://storage.googleapis.com/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499399.bam)

In this example, each file was downloaded to ~/Downloads. If you saved the downloaded files in another folder or location, the following paths will need to be updated to account for those differences.

```
mv ~/Downloads/hla_and_brca_genes.fa $PWD/ref
mv ~/Downloads/2895499331.bam $PWD/unaligned/normal
mv ~/Downloads/2895499399.bam $PWD/unaligned/normal
```

# Indexing

Index files of various formats and data structures are used commonly in genomics to provide random access to specific records, locations, or content within much larger domains and coordinate spaces. Ex. entire genome nucleotide sequences, billions of sequence reads, millions of variant records, etc.

## Samtools

Using `samtools faidx`, we create an index (.fai) of the reference FASTA (.fa) which provides random access to the nucleotides at specific positions within the complete genome reference. The FAIDX index is used both by the samtools faidx command as well as other toolkits, algorithms, and libraries that require random access to specific coordinates in a timely manner.

```
ls $PWD/ref
```

We only have the FASTA sequence, let's make a FAIDX format file using Samtools by running the command with no additional arguments (other than the FASTA file):

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools faidx /data/ref/hla_and_brca_genes.fa
```

We should now see the FASTA and FAIDX format files.

```
ls $PWD/ref
```

Take a look at the first twenty lines of this FASTA. What do you see? What do you NOT see?

```
head -n 20 $PWD/ref/hla_and_brca_genes.fa
```

Using `samtools faidx` once again, let's pass the name of the chromosome and some example start and end positions.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools faidx /data/ref/hla_and_brca_genes.fa chr17:43044295-43170245
```

Now, what do you see in the command-line output? Does this look more familiar than the head command run earlier?

## BWA

BWA index will apply the correct index algorithm depending on the size of the reference provided. The index actually consists of several files and data structures that all work together (and are required) to perform alignment with subsequent alignment commands, ex. mem, aln, sampe, etc.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/gitc/bwa index /data/ref/hla_and_brca_genes.fa

ls $PWD/ref
```

# Alignment

Let us demonstrate the alignment of two Normal HCC1395 read groups (from the same sample) using both FASTQ format files (created from our unaligned BAMs) as well as a performance improvement to align sequence reads directly from unaligned BAM to an aligned BAM in one streaming process (using pipes).

## Align from FASTQ

Let's walk through each step of the alignment process including the conversion from unaligned BAM to compressed FASTQ.

```
ls $PWD/unaligned/normal
```

The `-H` flag is provided to the `samtools view` subcommand so we can interogate the unaligned BAM header.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools view -H /data/unaligned/normal/2895499331.bam
```

Without the `-H` flag, the default behavior for the `samtools view` subcommand is to print each SAM record.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools view /data/unaligned/normal/2895499331.bam
```

Let's look at the SAM flag in this web app to decode what values are set for the first two reads: https://broadinstitute.github.io/picard/explain-flags.html

Now, we are going to use the `SamToFastq` command from the Picard toolkit to convert our unaligned BAM file to a pair of FASTQ files (the default input to BWA).

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 java -Xms2G -jar /usr/gitc/picard.jar SamToFastq 
```

The OUTPUT_PER_RG option is used to output one pair of FASTQ (both read1 and read2 for paired-end) to the output directory of our choosing. The output compressed FASTQ files will be named based on the read group ID.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 java -Xms2G -jar /usr/gitc/picard.jar SamToFastq \
        INPUT=/data/unaligned/normal/2895499331.bam OUTPUT_PER_RG=true COMPRESS_OUTPUTS_PER_RG=true RG_TAG=ID OUTPUT_DIR=/data/unaligned/normal

ls $PWD/unaligned/normal
```

Since the FASTQ files are compressed binary files, we can not view them as-is. As an example, we will decompress one of the read group files with gunzip. However, for an actual data set, keep the compressed binary format as-is to save space and rely on tools such as zcat, zgrep, zless, zdiff, etc. to perform basic manipulation (when necessary).

```
head $PWD/unaligned/normal/2895499331_1.fastq.gz
```

What are the characters above? Does this look like a FASTQ format file?

```
gunzip -f $PWD/unaligned/normal/2895499331_1.fastq.gz

ls $PWD/unaligned/normal

head $PWD/unaligned/normal/2895499331_1.fastq
```

How does the decompressed FASTQ look now? This format is much easier for humans to read, but no computers. Let's binary compress the file again to save space.

```
gzip $PWD/unaligned/normal/2895499331_1.fastq
```

The following steps are a series of BWA mem and Samtools view commands that:
- align the input FASTQ files to the FASTA reference
- redirect the output to a SAM format file
- convert SAM to BAM (compressed binary) file formats

```
ls $PWD/unaligned/normal

docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/gitc/bwa mem

docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/gitc/bwa mem /data/ref/hla_and_brca_genes.fa /data/unaligned/normal/2895499331_1.fastq.gz /data/unaligned/normal/2895499331_2.fastq.gz
```

What does the output of the above command look like? Is it a familiar file format we have discussed? What will we do with this output format? Would this scale for a whole genome?

Maybe we should re-run this command and output to a proper file. Let's take a look at `samtools view`, a basic command often used to read and write alignment formats.

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools view
```

When working with special characters, ex. the stdout redirect >, best practice is to run the command when quoted to avoid interpreting those characters in the original command, ie. docker run.

NOTE: The command is wrapped as a call out to the shell running within our container. This happens to be BASH, i.e., `/bin/bash -c`. The `-c` argument tells bash to run the command provided.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/gitc/bwa mem /data/ref/hla_and_brca_genes.fa /data/unaligned/normal/2895499331_1.fastq.gz /data/unaligned/normal/2895499331_2.fastq.gz > /data/aligned/normal/2895499331.sam"
```

We also need to clarify that the command in double quotes should be passed to the bash interpreter as the command to execute. We do so by specificly executing bash (`/bin/bash`) with the `-c` argument. `-c` is documented as: "Read and execute commands from the first non-option argument command_string, then exit."

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/gitc/bwa mem /data/ref/hla_and_brca_genes.fa /data/unaligned/normal/2895499331_1.fastq.gz /data/unaligned/normal/2895499331_2.fastq.gz > /data/aligned/normal/2895499331.sam" 

head $PWD/aligned/normal/2895499331.sam
```

Does the above output look familiar? What format is it? Will this scale?

Let's look at an example alignment with a binary compressed format output called BAM. To do so, we will pipe the SAM format directly to the samtools view command mentioned earlier. Note the longer command in one set of quotes.

```
!docker run --platform linux/amd64 -v $PWD:/data -v $PWD/aligned/normal:/data/aligned/normal -v $PWD/ref:/data/ref -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c '/usr/gitc/bwa mem -R "@RG\tID:2895499331\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.3\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" /data/ref/hla_and_brca_genes.fa /data/unaligned/normal/2895499331_1.fastq.gz /data/unaligned/normal/2895499331_2.fastq.gz | /usr/local/bin/samtools view -1 -o /data/aligned/normal/2895499331.bam -' 
```

Take a look at the header of the BAM format file we just created.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/local/bin/samtools view -H /data/aligned/normal/2895499331.bam" 
```

Now look at the alignments stored in the body of the BAM file.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/local/bin/samtools view /data/aligned/normal/2895499331.bam | head" 
```

## Align from Unaligned BAM

Let's walk through the alignment process again for another read group. This time we will "pipe" multiple commands together in one single streaming process from unaligned BAM to aligned BAM.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /usr/local/bin/samtools view -H /data/unaligned/normal/2895499399.bam
```

When running commands piped together, best practice is to set these bash arguments to generate appropriate error messages if something goes wrong such as an intermediate step is interupted or the streaming data becomes corrupt/incomplete.

```
set -o pipefail
set -o errexit

docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c \
    'java -Xms2G -jar /usr/gitc/picard.jar SamToFastq INPUT=/data/unaligned/normal/2895499399.bam FASTQ=/dev/stdout INTERLEAVE=true NON_PF=true | /usr/gitc/bwa mem -R "@RG\tID:2895499399\tPL:ILLUMINA\tPU:H7HY2CCXX-TGACCACG.4\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tSM:H_NJ-HCC1395-HCC1395_BL\tCN:MGI" -p /data/ref/hla_and_brca_genes.fa /dev/stdin | /usr/local/bin/samtools view -1 -o /data/aligned/normal/2895499399.bam -' 

ls $PWD/aligned/normal
```

Which file formats do you now see for read group 2895499399? As compared to 2895499331 the read group aligned earlier?

Take a closer look at the BAM output from the single, piped command:

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/local/bin/samtools view -H /data/aligned/normal/2895499399.bam" 

docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/local/bin/samtools view /data/aligned/normal/2895499399.bam | head" 
```

## Merge Alignments

Now that we have two read group BAMs, we want to `MergeSamFiles` (which works on BAM inputs) to produce a single BAM with all of our alignments for the Normal HCC1395 sample.

Next, we will look at the command-line usage for the Picard utility. Don't worry about ERROR message for a missing INPUT. We'll fix that soon.

```
docker run --platform linux/amd64 -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 java -Xms2G -jar /usr/gitc/picard.jar MergeSamFiles
```

Now, we will run the 'normal' merge using both of the newly minted read-group BAM files from earlier steps.

```
docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 java -Xms2G -jar /usr/gitc/picard.jar MergeSamFiles OUTPUT=/data/aligned/normal.bam INPUT=/data/aligned/normal/2895499331.bam INPUT=/data/aligned/normal/2895499399.bam

ls $PWD/aligned

docker run --platform linux/amd64 -v $PWD:/data -it broadinstitute/genomes-in-the-cloud:2.3.1-1512499786 /bin/bash -c "/usr/local/bin/samtools view -H /data/aligned/normal.bam" 
```

Do we really need to see another BAM header? Actually, take a closer look. How is this BAM header different from previous BAM headers we have looked at.

@RG and @PG tags - What do those stand for? How many tag lines do you see?

# Homework

- Index the normal.bam file. HINT: samtools index OR igvtools
- View the indexed normal.bam file with IGV HINT: Search for BRCA1.
- Make a list of questions and/or observations about the alignments to discuss next week.
- Are there other post-alignment processing steps we've missed? Bring suggestions for next week.

