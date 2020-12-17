## Week 12 - RNAseq part 1

- [Lecture Slides](https://github.com/genome/bfx-workshop/blob/master/week_12/RNASeq_IntrotoRNA.pdf)
- [Lecture Recording](https://wustl.box.com/s/6bguhi1d6su0ezv27glkd0cq0ea17d7v)

##Notes

### Obtain RNA-seq test data.
The test data consists of two commercially available RNA samples: [Universal Human Reference (UHR)](/assets/module_1/UHR.pdf) and [Human Brain Reference (HBR)](/assets/module_1/HBR.pdf). The UHR is total RNA isolated from a diverse set of 10 cancer cell lines. The HBR is total RNA isolated from the brains of 23 Caucasians, male and female, of varying age but mostly 60-80 years old.

In addition, a spike-in control was used. Specifically we added an aliquot of the [ERCC ExFold RNA Spike-In Control Mixes](/assets/module_1/ERCC.pdf) to each sample. The spike-in consists of 92 transcripts that are present in known concentrations across a wide abundance range (from very few copies to many copies). This range allows us to test the degree to which the RNA-seq assay (including all laboratory and analysis steps) accurately reflects the relative abundance of transcript species within a sample. There are two 'mixes' of these transcripts to allow an assessment of differential expression output between samples if you put one mix in each of your two comparisons. In our case, Mix1 was added to the UHR sample, and Mix2 was added to the HBR sample. We also have 3 complete experimental replicates for each sample. This allows us to assess the technical variability of our overall process of producing RNA-seq data in the lab.

For all libraries we prepared low-throughput (Set A) TruSeq Stranded Total RNA Sample Prep Kit libraries with Ribo-Zero Gold to remove both cytoplasmic and mitochondrial rRNA. Triplicate, indexed libraries were made starting with 100ng Agilent/Strategene Universal Human Reference total RNA and 100ng Ambion Human Brain Reference total RNA. The Universal Human Reference replicates received 2 ul of 1:1000 ERCC Mix 1. The Human Brain Reference replicates received 1:1000 ERCC Mix 2. The libraries were quantified with KAPA Library Quantification qPCR and adjusted to the appropriate concentration for sequencing. The triplicate, indexed libraries were then pooled prior to sequencing. Each pool of three replicate libraries were sequenced across 2 lanes of a HiSeq 2000 using paired-end sequence chemistry with 100bp read lengths.

So to summarize we have:

* UHR + ERCC Spike-In Mix1, Replicate 1
* UHR + ERCC Spike-In Mix1, Replicate 2
* UHR + ERCC Spike-In Mix1, Replicate 3
* HBR + ERCC Spike-In Mix2, Replicate 1
* HBR + ERCC Spike-In Mix2, Replicate 2
* HBR + ERCC Spike-In Mix2, Replicate 3

Each data set has a corresponding pair of FASTQ files (read 1 and read 2 of paired end reads).

```bash
cd ~
mkdir -p rnaseq/data
cd rnaseq/data
wget http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar

```

Unpack the test data using tar. You should see 6 sets of paired end fastq files. One for each of our sample replicates above. We have 6 pairs (12 files) because in fastq format, read 1 and read 2 of a each read pair (fragment) are stored in separate files.

```bash
tar -xvf HBR_UHR_ERCC_ds_5pc.tar
ls

```

Enter the data directory and view the first two read records of a file (in fastq format each read corresponds to 4 lines of data)

The reads are paired-end 101-mers generated on an Illumina HiSeq instrument. The test data has been pre-filtered for reads that appear to map to chromosome 22. Lets copy the raw input data to our tutorial working directory.
```bash
zcat < UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz | head -n 8

```

### Determining the strandedness of RNA-seq data 

In order to determine strandedness, we will be using [check_strandedness](https://github.com/betsig/how_are_we_stranded_here)([docker image](https://hub.docker.com/r/chrisamiller/how_stranded)). In order use this tool, there are a few steps we need to get our inputs ready, specifically creating a fasta of our GTF file.

Now that we have created our input files, we can now run the check_strandedness tool on some of our instrument data. Note: we are using a docker image for this tool.

```bash
cd ~/rnaseq
mkdir refs
cd refs
wget http://genomedata.org/rnaseq-tutorial/results/cshl2020/rnaseq/refs/chr22_with_ERCC92_tidy.gtf
wget http://genomedata.org/rnaseq-tutorial/results/cshl2020/rnaseq/refs/chr22_ERCC92_transcripts.clean.fa

docker run -v ~/rnaseq:/docker_workspace chrisamiller/how_stranded:latest check_strandedness --gtf /docker_workspace/refs/chr22_with_ERCC92_tidy.gtf --transcripts /docker_workspace/refs/chr22_ERCC92_transcripts.clean.fa --reads_1 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads_2 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```
`docker run` is how you initialize a docker container to run a command

`-v` is the parameter used to mount your workspace so that the docker container can see the files that you're working with. In the example above, `/home/ubuntu/workspace/rnaseq` from the EC2 instance has been mounted as `/docker_workspace` within the docker container. 

`chrisamiller/how_stranded` is the docker container name. The `:latest` refers to the specific tag and release of the docker container.


The output of this command should look like so:
```bash
Fraction of reads failed to determine: 0.1123
Fraction of reads explained by "1++,1--,2+-,2-+": 0.0155
Fraction of reads explained by "1+-,1-+,2++,2--": 0.8722
Over 75% of reads explained by "1+-,1-+,2++,2--"
Data is likely RF/fr-firststrand
```

Using this [table](https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/), we can see if this is what we expect. Note that since the UHR and HBR data were generated with the TruSeq Stranded Kit, as mentioned above, the correct strand setting for kallisto is `--rf-stranded`, which is what check_strandedness confirms. Similarly when we run HISAT we will use `--rna-strandness RF`, when we run StringTie we will use `--rf`, and when we run htseq-count we will use `--stranded reverse`.

### Get references and create HISAT2 index
Create a HISAT2 index for chr22 and the ERCC spike-in sequences. HISAT2 can incorporate exons and splice sites into the index file for alignment. First create a splice site file, then an exon file. Finally make the aligner FM index.

To learn more about how the HISAT2 indexing strategy is distinct from other next gen aligners refer to the HISAT publication.

```bash
cd ~/rnaseq/refs/
wget http://genomedata.org/rnaseq-tutorial/fasta/GRCh38/chr22_with_ERCC92.fa
wget http://genomedata.org/rnaseq-tutorial/annotations/GRCh38/chr22_with_ERCC92.gtf


docker run -v ~/rnaseq:/docker_workspace zlskidmore/hisat2:latest hisat2_extract_splice_sites.py /docker_workspace/refs/chr22_with_ERCC92.gtf > /docker_workspace/refs/splicesites.tsv
docker run -v ~/rnaseq:/docker_workspace zlskidmore/hisat2:latest hisat2_extract_exons.py /docker_workspace/refs/chr22_with_ERCC92.gtf > /docker_workspace/refs/exons.tsv
docker run -v ~/rnaseq:/docker_workspace zlskidmore/hisat2:latest hisat2-build -p 4 --ss /docker_workspace/refs/splicesites.tsv --exon /docker_workspace/refs/exons.tsv /docker_workspace/refs/chr22_with_ERCC92.fa /docker_workspace/refs/chr22_with_ERCC92
ls
```

### HISAT2 alignment
Perform alignments with HISAT2 to the genome and transcriptome.

First, begin by making the appropriate output directory for our alignment results.

```bash
cd ~/rnaseq
mkdir alignment
cd alignment
```

HISAT2 uses a graph-based alignment and has succeeded HISAT and TOPHAT2. The output of this step will be a SAM/BAM file for each data set.

Refer to HISAT2 manual for a more detailed explanation:

* [https://ccb.jhu.edu/software/hisat2/manual.shtml](https://ccb.jhu.edu/software/hisat2/manual.shtml)

HISAT2 basic usage:

```bash
#hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r> | --sra-acc <SRA accession number>} [-S <sam>]
```

Extra options specified below:

* '-p 8' tells HISAT2 to use eight CPUs for bowtie alignments.
* '--rna-strandness RF' specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries. See here for options.
* '--rg-id $ID' specifies a read group ID that is a unique identifier.
* '--rg SM:$SAMPLE_NAME' specifies a read group sample name. This together with rg-id will allow you to determine which reads came from which sample in the merged bam later on.
* '--rg LB:$LIBRARY_NAME' specifies a read group library name. This together with rg-id will allow you to determine which reads came from which library in the merged bam later on.
* '--rg PL:ILLUMINA' specifies a read group sequencing platform.
* '--rg PU:$PLATFORM_UNIT' specifies a read group sequencing platform unit. Typically this consists of FLOWCELL-BARCODE.LANE
* '--dta' Reports alignments tailored for transcript assemblers.
* '-x /path/to/hisat2/index' The HISAT2 index filename prefix (minus the trailing .X.ht2) built earlier including splice sites and exons.
* '-1 /path/to/read1.fastq.gz' The read 1 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
* '-2 /path/to/read2.fastq.gz' The read 2 FASTQ file, optionally gzip(.gz) or bzip2(.bz2) compressed.
* '-S /path/to/output.sam' The output SAM format text file of alignments.

```bash
hisat2 -p 8 --rg-id=UHR_Rep1 --rg SM:UHR --rg LB:UHR_Rep1_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep1.sam
hisat2 -p 8 --rg-id=UHR_Rep2 --rg SM:UHR --rg LB:UHR_Rep2_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep2.sam
hisat2 -p 8 --rg-id=UHR_Rep3 --rg SM:UHR --rg LB:UHR_Rep3_ERCC-Mix1 --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./UHR_Rep3.sam

hisat2 -p 8 --rg-id=HBR_Rep1 --rg SM:HBR --rg LB:HBR_Rep1_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep1.sam
hisat2 -p 8 --rg-id=HBR_Rep2 --rg SM:HBR --rg LB:HBR_Rep2_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep2.sam
hisat2 -p 8 --rg-id=HBR_Rep3 --rg SM:HBR --rg LB:HBR_Rep3_ERCC-Mix2 --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz -2 $RNA_DATA_DIR/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz -S ./HBR_Rep3.sam

```
