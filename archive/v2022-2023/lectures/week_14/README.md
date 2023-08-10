# Week 14
## RNAseq Part 1

- [Lecture Recording](https://wustl.box.com/s/g0b3tz7xcb8y0hegy3npxaz0wlkgxd45)

- [Slides](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2022/full/RNASeq_Module1_IntrotoRNA.pdf)

- [RNA-seq Bioinformatics](https://rnabio.org/course)

## RNAseq Course Setup
For the BFX Workshop, we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org.

### Docker Setup

A Docker image is available in the ICTS Precision Health Artifact Registry on Google Cloud. For more details, see the console:
https://console.cloud.google.com/artifacts/docker/icts-precision-health/us-central1/bfx-workshop-repo/ami-09b613ae9751a96b1?project=icts-precision-health

1. Use these commands with the Docker client to pull the image. To use these commands, your Docker client must be configured to authenticate with `us-central1-docker.pkg.dev`. If this is the first time that you are pulling an image from `us-central1-docker.pkg.dev` with your Docker client, run the following command on the machine where Docker is installed. You may have to install the [gcloud cli](https://cloud.google.com/sdk/docs/install) if you have not previously.

```bash
gcloud auth configure-docker us-central1-docker.pkg.dev
```

2. This command will pull the image `ami-09b613ae9751a96b1` to your local Docker client with the tag `0.0.1` from the `bfx-workshop-repo` repository in the `icts-precision-health` project:

```bash
docker pull \
    us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/ami-09b613ae9751a96b1:0.0.1
```

3. Setup a local workspace directory for the RNAseq course. If you change the path or command used here in Step 3, please update the path to the workspace directory accordingly in Step 4.

```bash
mkdir -p ~/rnabio-workspace
```

4. Initialize a Docker container using the image we pulled above. `-v` tells Docker to mount our workspace directory within the Docker container as `/workspace` with read-write priveleges. You'll see in the RNAseq course `/workspace` is the base directory for nearly all commands and steps.

```bash
docker run -v ~/rnabio-workspace:/workspace:rw -it us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/ami-09b613ae9751a96b1:0.0.1 /bin/bash
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
1. `geneBody_coverage.py` in the optional RSeQC section is not correctly in the `PATH`. Use the full path to the python script `/home/ubuntu/.local/bin/geneBody_coverage.py`


## Homework Assignments
1. Finish Module 1 - [Inputs](https://rnabio.org/module-01-inputs/0001/01/01/Intro_to_Inputs/)
2. Complete Module 2 - [Alignment](https://rnabio.org/module-02-alignment/0002/01/01/Intro_to_Alignment/)
