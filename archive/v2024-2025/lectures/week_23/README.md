# Week 23 - Long Read Alignment

[Lecture slides](long_read_sequencing.pdf)

[Lecture video](https://wustl.box.com/s/pamv6vzs8sy6laexddyep6xfpo4u31ga) 


## Assignment

Let's look at some outputs of long read sequencing from the Oxford Nanopore (ONT) platform. These are sequences from the K562 cell line, prepared with the ONT cDNA sequencing kit (poly-A selected).  Off the machines, the data will consist of a FAST5 or POD5 file, which are a compressed representation of the raw signal. These are subsequently run through a basecalling algorithm (such as [Dorado](https://github.com/nanoporetech/dorado)) to generate FASTQ files.

The choice of basecalling algorithm and parameters goes pretty deep, so we'll assume that reasonable choices have been made. For simplicity, we've also subset the data to include just small portions of the genome, including a few genes of interest. 

Go ahead and pull down this fastq file:

```
wget https://storage.googleapis.com/bfx_workshop_tmp/k562_ont_raw.fastq.gz
```


To start, let's go ahead and trim the data with the [pychopper](https://github.com/epi2me-labs/pychopper) tool from ONT. It's used to identify, orient, trim, and un-fuse reads.  We don't have a local install of this tool, so let's use docker to run it:

```
docker run -it -v $(pwd -P):/data ontresearch/wf-isoforms
cd /data
pychopper -t 4 -m edlib -k PCS111 -r pychopper_report.pdf -S pychopper_stats.txt k562_ont_raw.fastq.gz | gzip > k562_ont_trimmed.fastq.gz
```

In this case, we've told pychopper to use 4 threads (`-t`), its standard model for deconvoluting reads/primers (`-k PCS111`) and to dump out a report and some stats.

After that finishes, exit the docker container, and look at the output reports, in the forms of some statistics and some plots about them.  In this case, the data has been pre-selected, and so the vast majority of the reads are usable, but in some real-life runs, an even larger fraction of them might not be.

### Sequence QC
Next, do some basic QC on this fastq using a tool called [NanoPlot](https://github.com/wdecoster/NanoPlot).  

```
docker run -v $(pwd -P):/data -it nanozoo/nanoplot
cd /data
mkdir nanoplot
NanoPlot --fastq k562_ont_trimmed.fastq.gz --prefix nanoplot/
```
In this case, we've told nanoplot to prepend `nanoplot/` to all of the output file names, meaning that they'll all land in that folder we just created.

Exit the container and open the `NanoPlot-report.html` file. Note that our reads here are quite long compared to short read data - over 1000 bp long.  ONT typically advertises reads reaching tens of thousands of base pairs, though - why isn't this the case here?

<details><summary>Answer</summary>
<p>
This data is created from RNA, which means that the lengths of the molecules are dependent on the lengths of the transcripts, which are not typically tens of thousands of bases long
</p>
</details>

Scrolling down, you can see this length, along with other stats, represented graphically with embedded dynamic plots, or you can retrieve simple png images from the output folder.  Also note the read mean read quality, which is substantially lower than what we'd see with Illumina short-read sequencing.  

Looking at the read lengths plots, they look a little multimodal.  In a real experiment, that would be weird, but nothing to worry about in this case - it's an artifact that appears because we're only using a very small amount of data.


### Aligning the reads
Next up, let's align this data, using a tool called minimap2. It's generally the go-to aligner for long-read data, whether from RNA or DNA. 

Just like with short-read DNA data, we'll need a reference genome to align against. This fasta contains a small portion of the human build 38 reference. Download and untar it:

```
wget https://storage.googleapis.com/bfx_workshop_tmp/ref_genome_sm.tar
tar -xvf ref_genome_sm.tar
```

Your command for running the alignment should look like this:

```
docker run -it -v $(pwd -P):/data quay.io/biocontainers/minimap2:2.17--hed695b0_3
cd /data
minimap2 -ax splice -uf ref_genome_sm.fa.gz k562_ont_trimmed.fastq.gz >k562_ont.sam
```

This will take a while to align, so plan on letting it go in the background for a bit. While that's running, a few things to notice here: 

1. We don't need to create an index ahead of time for minimap. It creates its index of the reference genome very quickly and efficiently compared to short-read aligners like bwa.

2. We're using the `-x splice` parameter, which indicates that it should be doing spliced alignment. This handles RNAseq data that has to account for introns interrupting the reads. We'll talk in more detail about spliced alignment that during the RNAseq portion of this course.  We're also using the `-uf` flags which help it find exon boundaries more effectively, by looking at the sequence context.

Now, we've got an aligned SAM file, and you can look at it with `less` and see that it follows the familiar format of header, followed by alignments.  As before, we're going to want to sort and index this bam in order to make it useful for downstream steps.

```
docker run -it -v $(pwd -P):/data biocontainers/samtools:v1.9-4-deb_cv1
cd /data
samtools view -Sb k562_ont.sam | samtools sort >k562_ont.bam
samtools index k562_ont.bam
```

Exit the container, then open IGV. It's time to visualize this data and get a feel for what you've got going on.  Let's grab a short-read Illumina bam, generated from the same cell line, for comparison:

```
wget https://storage.googleapis.com/bfx_workshop_tmp/k562_illumina.bam https://storage.googleapis.com/bfx_workshop_tmp/k562_illumina.bam https://storage.googleapis.com/bfx_workshop_tmp/k562_illumina.bam.bai

```

Open both the Illumina and ONT bams in IGV, using the Human (hg38) reference genome.  Navigate to the `U2AF1` gene and answer some questions:


1) How many exons are spanned by a typical short read?  How about a typical long read?  Which has more even coverage over the exons of these genes?

2) Right click on the Gene track and choose "Expanded"  Can you match up individual long reads with full-length transcripts?  What about short reads? How might this affect your ability to understand alternative splicing?

3) Zoom into an exon, until you can see basepair level changes.  Which data has a higher error rate?  Are the classes of errors you see between the two platforms different?  What kinds of analyses would be limited by this error rate?  

4) Now navigate to the DNMT3A gene and zoom in on the 3' (left) end.  What do you notice about the reads and the coverage?  Do these look like full-length transcripts to you?

5) Go to the ensembl website and pull look at the gold transcripts for DNMT3A - the most stable, well characterized.  How long is the the dominant isoform? Does that information help explain what you see in IGV?

Compare this to the short read data. Which is more useful in this situation? Think about what kinds of biases this can introduce when estimating gene expression or transcript abundance.
