# Week 7
## RNAseq Part 1 - Introduction and Alignment

- Lecture Recording

- [Slides](RNASeq_part1_bfx.pdf)

- [RNA-seq Bioinformatics home](https://rnabio.org/course) (rnabio.org)

## Assignments:
1. Complete the short "RNAseq Course Setup" section below to get your computer ready to do RNAseq analysis.

2. Finish Module 1 of the RNAbio course: [Inputs](https://rnabio.org/module-01-inputs/0001/01/01/Intro_to_Inputs/)
    - Be sure to note the "Known issues" section at the bottom of this page, which has you modify a key check-strandedness command in #1. 
3. Finish Module 2 of the RNAbio course: [Alignment](https://rnabio.org/module-02-alignment/0002/01/01/Intro_to_Alignment/)
    - Skip the "IGV" section, as we've covered most of that in class previously. 
    - The "Team Assignment - Alignment" section is optional. Feel free to complete it if you want a challenge!
    
----

## RNAseq Course Setup
For the BFX Workshop, we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org. 

1) Pull the image to your local Docker client from the griffithlab repository:

```
    docker pull griffithlab/rnaseq-toolkit:latest
```

2) Set up a local workspace directory for the RNAseq course. If you change the path or command used here in Step 2, make sure to update the path to the workspace directory accordingly in Step 3.

```
    mkdir -p ~/workshop/rnabio-workspace/rnaseq
```

3) Initialize a Docker container using the image we pulled above

```
    docker run -it -p 8080:8080 -v ~/workshop/rnabio-workspace:/workspace griffithlab/rnaseq-toolkit:latest
```

Notes about this command:
- in our previous docker exercises, we've often mounted things to `/data/`, but here, we're following the convention of this course, and mounting them to `/workspace`.  This is the base directory for nearly all commands and steps in RNAbio.
- we're running ths docker image interactively and opening port 8080 so that we can access files in specific ways.


--- 
## Known Issues/ Discrepancies from RNAbio website
1. When running the check strandedness tool in the Module 1, RNAseq Data section, the docker run command cannot be run from within your docker session. To run it, we suggest that you open a new terminal window, `cd` into the `rnaseq` directory you created at the beginning of this assignment, and use the following command instead:
```
cd ~/workshop/rnabio-workspace/rnaseq
docker run -v $PWD/:/docker_workspace mgibio/checkstrandedness:latest check_strandedness --gtf /docker_workspace/refs/chr22_with_ERCC92_tidy.gtf --transcripts /docker_workspace/refs/chr22_ERCC92_transcripts.clean.fa --reads_1 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads_2 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```
This is the same command as what is mentioned in the course webpage, except that instead of mounting (`-v` flag) `/home/ubuntu/workspace/rnaseq` to the docker image- which is where the data was stored for students running through the course on an AMI, you will instead mount whatever your current directory is. Also, this is different from an interactive session where we are able to enter the docker and run commands within it. Instead we are executing our command directly all in that one line of code.

2. In various parts of RNAbio, in order to view HTML files, plots, etc., the tutorial suggests going to a public IPV4 address link in your browser window. That is only needed for the AMI. Since you'll be running everything locally, you can either find the files in your Finder window or File Explorer and open them directly; or even better, use `open [your_file.html]` on Mac and `explorer.exe [your_file.html]` on Windows/WSL2 to open the file in your default browser!
