# Week 10
## RNAseq Part 1

- [Lecture Recording] (Link coming soon)

- [Slides](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2023/full/RNASeq_Module1_IntrotoRNA.pdf)

- [RNA-seq Bioinformatics](https://rnabio.org/course)

## RNAseq Course Setup
For the BFX Workshop, we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org.

### Docker Setup

A Docker image is available through the DockerHub repository-
[https://hub.docker.com/layers/griffithlab/rnabio/0.0.1/images/sha256-b13f5e9048941c8be3e83555295c0f4ed21645d5fd9bae4226e6bc4f30f54b52?context=explore](https://hub.docker.com/layers/griffithlab/rnabio/0.0.1/images/sha256-b13f5e9048941c8be3e83555295c0f4ed21645d5fd9bae4226e6bc4f30f54b52?context=explore)

1. Ensure that Docker Desktop is running. 

2. This command will pull the image `rnabio` to your local Docker client with the tag `0.0.1` from the `griffithlab` DockerHub repository:

```bash
docker pull griffithlab/rnabio:0.0.1
```

3. Setup a local workspace directory for the RNAseq course. If you change the path or command used here in Step 3, please update the path to the workspace directory accordingly in Step 4. Also, make a file `test_my_docker_mount` that we will look for later.

```bash
mkdir -p bfx-workshop/rnabio-workspace
echo 'this file helps me test my docker mount' >> bfx-workshop/rnabio-workspace/test_my_docker_mount
```

4. Enter the directory where you created the `rnabio-workspace` folder, and initialize a Docker container using the image we pulled above. `-v` tells Docker to mount our workspace directory within the Docker container as `/workspace` with read-write priveleges. You'll see in the RNAseq course `/workspace` is the base directory for nearly all commands and steps.

```bash
cd bfx-workshop/rnabio-workspace
docker run -v $PWD/:/workspace:rw -it griffithlab/rnabio:0.0.1 /bin/bash
```

5. Use `ls` to see what's in this file, enter the `workspace` folder and then use `ls` again to see what is in the `workspace` folder.

```bash
ls
cd workspace
ls
```

### User Setup

Now that we are running a Docker container, Docker, by default, will log you in as the "root" user. We need to run as the ubuntu user to match the RNAseq course tutorials.

1. Switch User `su` to the unbutu user:

```bash
su ubuntu
```

2. Source the pre-installed `.bashrc` file to configure your environment to match the RNAseq course:

```bash
source ~/.bashrc
```

NOTE: Using Docker and the persistent "workspace" volume we attached will allow you to start/stop as you wish. EVERYTIME YOU LOGIN TO THE DOCKER CONTAINER, YOU MUST LOGIN AS THE `ubuntu` USER *AND* `source ~/.bashrc` UPON EACH LOGIN.

### Environment Setup

Create a working directory and set the ‘RNA_HOME’ environment variable
```
mkdir -p ~/workspace/rnaseq/

export RNA_HOME=~/workspace/rnaseq
```

Make sure whatever the working dir is, that it is set and is valid
```
echo $RNA_HOME
```

Since all the environment variables we set up for the RNA-seq workshop start with ‘RNA’ we can easily view them all by combined use of the env and grep commands as shown below. The env command shows all environment variables currently defined and the grep command identifies string matches.
```
env | grep RNA
```

In order to view the contents of this file, you can type:
```
less ~/.bashrc
```
To exit the file, type `q`.

### Known Issues
1. When running the check strandedness tool in the Module 1, RNAseq Data section, the docker run command cannot be run from within your rnabio docker session. To run it, we suggest entering the `rnaseq` directory you created at the beginning of the course, and using the following command instead-
```
docker run -v $PWD/:/docker_workspace mgibio/checkstrandedness:latest check_strandedness --gtf /docker_workspace/refs/chr22_with_ERCC92_tidy.gtf --transcripts /docker_workspace/refs/chr22_ERCC92_transcripts.clean.fa --reads_1 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz --reads_2 /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```
 
2. `geneBody_coverage.py` in the optional RSeQC section is not correctly in the `PATH`. Use the full path to the python script `/home/ubuntu/.local/bin/geneBody_coverage.py`


## Homework Assignments
1. Finish Module 1 - [Inputs](https://rnabio.org/module-01-inputs/0001/01/01/Intro_to_Inputs/)
2. Complete Module 2 - [Alignment](https://rnabio.org/module-02-alignment/0002/01/01/Intro_to_Alignment/)
