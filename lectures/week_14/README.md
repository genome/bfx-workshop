# Week 14
## RNAseq Part 1
TODO: Overview

- [Lecture Recording]()

- [Slides]()

- [RNA-seq Bioinformatics](https://rnabio.org/course()

## RNAseq Course Setup
For the BFX Workshop we will not be using AWS Cloud. Instead, we will use a Docker image created from the AWS AMI used in rnabio.org.

### Docker Setup

A Docker image is available in the ICTS Precision Health Artifact Registry on Google Cloud. For more details, see the console:
https://console.cloud.google.com/artifacts/docker/icts-precision-health/us-central1/bfx-workshop-repo/ami-09b613ae9751a96b1?project=icts-precision-health

1. Use these commands with the Docker client to pull the image. To use these commands, your Docker client must be configured to authenticate with `us-central1-docker.pkg.dev`. If this is the first time that you are pulling an image from `us-central1-docker.pkg.dev` with your Docker client, run the following command on the machine where Docker is installed.
```
gcloud auth configure-docker us-central1-docker.pkg.dev
```
2. This command will pull the image `ami-09b613ae9751a96b1` to your local Docker client with the tag `0.0.1` from the `bfx-workshop-repo` repository in the `icts-precision-health` project:
```
docker pull \
    us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/ami-09b613ae9751a96b1:0.0.1
```
3. Setup a local workspace directory for the RNAseq course. If you change the path or command used here in Step 3, please update the path to the workspace directory accordingly in Step 4.
```
mkdir -p ~/rnabio-workspace
```

4. Initialize a Docker container using the image we pulled above. `-v` tells Docker to mount our workspace directory within the Docker container as `/workspace` with read-write priveleges. You'll see in the RNAseq course `/workspace` is the base directory for nearly all commands and steps.
```
docker run -v ~/rnabio-workspace:/workspace:rw -it us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/ami-09b613ae9751a96b1:0.0.1 /bin/bash
```

### User Setup

Now that we are running a Docker container, Docker, by default, will log you in as the "root" user. We need to run as the ubuntu user to match the RNAseq course tutorials.

1. Switch User `su` to the unbutu user:
```
su ubuntu
```
2. Source the pre-installed `.bashrc` file to configure your environment to match the RNAseq course:
```
source ~/.bashrc
```
NOTE: Using Docker and the persistent "workspace" volume we attached will allow you to start/stop as you wish. EVERYTIME YOU LOGIN TO THE DOCKER CONTAINER, YOU MUST LOGIN AS THE `ubuntu` USER.

## Homework Assignments
1. TODO
2. TODO
3. TODO
4. TODO
