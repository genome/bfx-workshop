# Week 12

## RNAseq Part 3

- [Lecture Recording](https://wustl.box.com/s/1frgkv9qa1vc9xcnlqs7w8g6pjbrg3e4)

- [Differential Expression Slides](bfx_RNASeq_Module3_2024.pdf)

- [RNA-seq Bioinformatics](https://rnabio.org/course)

## RNAseq Course Setup

For the BFX Workshop, we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org.

### Docker Setup

Instructions for setting up docker can be found [here](https://github.com/genome/bfx-workshop/tree/master/lectures/week_07)

### User Setup

**Reminder:** Now that we are running a Docker container, Docker, by default, will log you in as the "root" user. We need to run as the ubuntu user to match the RNAseq course tutorials.

1. Switch User `su` to the unbutu user:

```bash
su ubuntu
```

2. Source the pre-installed `.bashrc` file to configure your environment to match the RNAseq course:

```bash
source ~/.bashrc
```
3. Set the environment variable
```bash
export RNA_HOME=~/workspace/rnaseq
``` 

NOTE: Using Docker and the persistent "workspace" volume we attached will allow you to start/stop as you wish. EVERYTIME YOU LOGIN TO THE DOCKER CONTAINER, YOU MUST LOGIN AS THE `ubuntu` USER *AND* `source ~/.bashrc` UPON EACH LOGIN.

## Homework Assignments

Start with [Differential Expression Ballgown](https://rnabio.org/module-03-expression/0003/03/01/Differential_Expression-Ballgown/), and continue through the pages, finishing with "DE Visualization Advanced".

### For-credit students
Please provide a screenshot of a volcano plot representing the differentially expressed genes in the dataset. Set the p-value cutoff to 0.01 and the log fold change threshold to 2.  Send it to Jenny, as usual.




