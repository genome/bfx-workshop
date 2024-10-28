# Week 8

## RNAseq Part 2

- Lecture Recording

- [Abundance Estimation Slides](RNASeq_Module3_AbundanceEstimation_bfxworkshop_2024.pdf)

- [RNA-seq Bioinformatics](https://rnabio.org/course)

## RNAseq Course Setup

For the BFX Workshop, we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org.

### Docker Setup

Instructions for setting up docker can be found [in last week lecture](https://github.com/genome/bfx-workshop/tree/master/lectures/week_07)

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

NOTE: Using Docker and the persistent "workspace" volume we attached will allow you to start/stop as you wish. EVERYTIME YOU LOGIN TO THE DOCKER CONTAINER, YOU MUST LOGIN AS THE `ubuntu` USER *AND* `source ~/.bashrc` UPON EACH LOGIN.

## Homework Assignments

Complete Expression Analysis portion of Module 3 - [Expression](https://rnabio.org/module-03-expression/0003/02/01/Expression/)

For-credit students: please find the TPM value of gene CECR7 in sample HBR_Rep1 and send that to Jenny. 
