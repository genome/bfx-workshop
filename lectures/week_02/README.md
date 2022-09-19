## Week 2: DNA Sequence generation and QC

This week, we cover some basics of how sequence is generated, delve into the details of data and formats, and talk through basic QC of your sequence.

- [Lecture Recording](https://wustl.box.com/s/yev38ztocy7uqn1tvwu7woskkswmlluo)

- [Slides](bfx_workshop_02_sequencing.pdf)

## Homework Assignment

Let's play around with some real sequence data and QC it.  See the lecture for some commands that may help, and ask for help in #bfx_workshop if you get stuck!

1) Get an interactive job in a docker container. The image `chrisamiller/docker-genomic-analysis` has a lot of common genomics tools installed that may be useful for the first few steps

2) We're going to work with data from a human cell line posted here: [https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar](https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar) Make a directory called "week2", and download the tarball to the cluster using the command line (`wget` or `curl`).

3) Use `tar -xvf` to extract the directory, then cd into the directory and look around.  We're not going to use all of this data in this week's homework. Let's focus on the contents of `Exome_Tumor.tar`. Untar it, then unzip the fastq files.

4) Look at the first three records of each fastq file. Take a close look at the read names and how they match up across files. (Use the unix commands that we discussed in class to answer this and the following three questions)

5) How many paired end sequences do these files contain? 

6) What is the read length? Is the read length consistent for every record?

7) How many total nucleotides of sequence are contained in these two files?

8) Exit and get into an interactive container that contains the fastqc tool: `quay.io/biocontainers/fastqc:0.11.9--0`

9) Run fastqc on these data:  `fastqc *.fastq.gz`  What is the asterisk doing here?  Note the files it produces - an HTML file, with a user-friendly summary, and a zip file, which you can dig into if you need more details, or wanted to parse the files by hand.

10) Exit the docker image, browse to the html files, and open them up.  Do you see any potential issues with the sequence data?
