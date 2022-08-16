## Week 12 - RNAseq part 1

- [Lecture Slides](RNASeq_IntrotoRNA.pdf)
- [Lecture Recording](https://wustl.box.com/s/6bguhi1d6su0ezv27glkd0cq0ea17d7v)

## Assignment

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


