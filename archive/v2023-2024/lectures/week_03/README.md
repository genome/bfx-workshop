## Week 3: DNA Sequence generation and QC

This week, we cover some basics of how sequence is generated, delve into the details of data and formats, and talk through basic QC of your sequence.

- [Lecture Recording](https://wustl.box.com/s/4toay23fvzanyg09k3blbtlf3wp2qurs)

- [Slides](bfx_workshop_03_sequencing.pdf)

## Homework Assignment

Let's play around with some real sequence data and QC it. The expectation is that you'll do this on your laptop/desktop, but it should also be possible to do all of this on the cluster or cloud, provided that docker is installed. 

This is not a step-by-step tutorial, but all of the commands you need to complete these steps are introduced in the lecture slides. Remember to use the `-h` or `--help` flags to get usage, or use `man <command>` to see additional help on particular commands. Ask for help in #bfx_workshop if you get stuck!

1) Get an interactive job in a docker container. The image `chrisamiller/docker-genomic-analysis` has a lot of common genomics tools installed that may be useful for the first few steps. Be sure to mount your working directory so that you have access inside the container!

2) We're going to work with data from a human cell line posted here: [https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar](https://storage.googleapis.com/bfx_workshop_tmp/Exome_Tumor.tar) Make a directory called "week2", and download the tarball to your computer using the command line (`wget` or `curl`).

3) Use `tar -xvf` to extract the directory from the tar file, then cd into the directory and look around.  We're not going to use all of this data in this week's homework. Let's focus on the contents of `Exome_Tumor.tar`. Untar it, then unzip the fastq files.

4) Look at the first three records (not first three lines!) of each fastq file. Take a close look at the read names and how they match up across files. 

5) How many paired end sequences do these files contain?

6) What is the read length? Is the read length consistent for every record?

7) How many total nucleotides of sequence are contained in these two files?

8) Use `gzip` to recompress these two fastq files to save space

9) Exit that docker container (type `exit`) and launch a new docker session using a container that has the fastqc tool: `quay.io/biocontainers/fastqc:0.11.9--0`

10) Run fastqc on these data:  `fastqc *.fastq.gz`  What is the asterisk doing here?  Note the files it produces - an HTML file, with a user-friendly summary, and a zip file, which you can dig into if you need more details, or wanted to parse the files by hand.

11) Exit the docker image, browse to the html files, and open them up.  Do you see any potential issues with the sequence data?
