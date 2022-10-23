# Cloud Build Docker Tutorial

Let's build a docker container, but this time, we'll build it and host it in our cloud space, so that it's accessible when we try to run things there. 

1. Cloud Shell Workspace - First, let's set up our workspace:

Log in to a Cloud Shell Terminal: https://shell.cloud.google.com/?show=ide%2Cterminal

Authenticate with WUSTL Key account.
```
gcloud auth login
```

Set the Google Project we are using.
```
gcloud config set project icts-precision-health
```

Show the value for the environment variable `$USER`
```
echo $USER
```

Parts of GCP don't allow characters like underscores that are common in wustl user names. Set a new "clean" shell variable to remove invalid characters:
```
USER_CLEAN=`echo $USER | sed 's/_/-/g'` && export USER_CLEAN
```

Create a directory and setup a cloudshell workspace using the new directory:
```
mkdir gatk-depth-filter-docker
```

```
cloudshell workspace gatk-depth-filter-docker
```
With the buttons at the top right, you can toggle between the Terminal (using `Open/Close Terminal`) and the Editor (using `Open/Close Editor`). You can also leave them both open at the same time.

In a Cloud Shell Terminal, navigate into the newly created directory:
```
cd gatk-depth-filter-docker
``` 

2. Python Script - we need to pull in the python script that we want to include next to GATK in our container.

Using the Cloud Shell Editor, Save a file named `depth_filter.py` in the `gatk-depth-filter-docker` workspace. 

Copy/Paste the contents of `depth_filter.py` from GitHub:
https://raw.githubusercontent.com/genome/docker-depth-filter/master/depth_filter.py

In a Cloud Shell Terminal, add executable permissions to the Python script:
```
chmod +x depth_filter.py
```

(What other terminal commands could you have used to get this file instead of a copy/paste?)

3. Dockerfile - the script that docker will use to create the image

Using the Cloud Shell Editor, Save a file named `Dockerfile` in the `gatk-depth-filter-docker` workspace containing:
```
FROM broadinstitute/gatk:4.3.0.0
RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all
RUN pip install vcfpy pysam
COPY depth_filter.py /usr/bin/depth_filter.py
```

Note that this essentially the same thing we used in our local install, except that we added a library needed for the cloud (libnss-sss).  Generally speaking, a Dockerfile should be platform independent, even though the steps around building them might differ slightly.

4. Cloud Build - create the docker image

On our local machines, we could just say "docker build" to build an image, then "docker push" to shoot it up to dockerhub (assuming you've created an account there and are logged in).  Here on GCP, we need to use our own repository/registry, and set up a config file to tell the build process how to push it there for later use.

Using the Cloud Shell Editor, Save a file named `cloudbuild.yaml` in the `gatk-depth-filter-docker` workspace containing:
NOTE: REPLACE `$USER_CLEAN` WITH THE VALUE FROM THE ENVIRONMENT VARIABLE. This should be your WUSTL Key with `-` instead of `_` characters.
```
steps:
- name: 'gcr.io/cloud-builders/docker'
  args: [ 'build', '-t', 'us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:0.1', '.' ]
images:
- 'us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:0.1'
```
Note that we've given it a unique name ($USER_CLEAN-gatk-depth-filter-image) and a version tag (0.1)

In a Cloud Shell Terminal, submit the build:
```
gcloud builds submit --region=us-central1 --config cloudbuild.yaml
```

If you hit errors like `ERROR: (gcloud.builds.submit) parsing cloudbuild.yaml: while parsing a block collection`, check that your yaml file's spacing and indentation are exactly the same as the example above. 

You can follow your build's progress on the command line, or by looking at the "History" tab in Google Cloud Build: https://console.cloud.google.com/cloud-build/dashboard;region=us-central1?project=icts-precision-health

5. Docker Run - actually use your container

First, from a Cloud Shell Terminal, we need to configure the Artifact Registry credentials for the region we intend to pull the Docker image from.
```
gcloud auth configure-docker us-central1-docker.pkg.dev
```

Now, we can use the Docker `run` command to pulll the image from the Artifact Registry and jump into the container once it's running with Docker.
```
docker run -it us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:tag1 /bin/bash
```

# Germline Variant Detection
We will repeat the Germline Variant Detection steps from the local tutorial. 

However, this time we will perform all of the steps in the cloud including the execution of the depth_filter.py script and upload of the VCF back to a Cloud Bucket for use with IGV.

Please be sure you are using the Docker container terminal (in Cloud Shell) from Step 5 of the Cloud Build tutorial. If not, load the container using Step 5 in Cloud Shell now.

## GATK

```
gatk --list
```

```
gatk HaplotypeCaller --help
```
## Input
Use the cromwell-server API to retrieve the path to the HCC1395 Normal BAM file (ex. `alignment_exome.cwl.bam`).

Replace $BAM with the path returned from Cromwell. 

Example: `gs://icts-precision-health-cromwell-wf-exec/alignment_exome.cwl/d153b0da-ef5e-43ba-94dc-69e24311c83f/call-alignment/sequence_to_bqsr.cwl/a1da9f01-7014-411b-baef-e592bcf34cb6/call-index_bam/H_NJ-HCC1395-HCC1395_BL.bam`

## Output
We are now in a Docker container which does not inherit the environment from our Cloud Shell Terminal.
NOTE: Replace `$USER_CLEAN` in the path below
NOTE: Replace `$BAM` with the path returned by cromwell-server in the [DNA Alignment Workflow Tutorial](../week_06/bfx_workshop_06_alignment.md). 
```
gatk HaplotypeCaller --input $BAM --output gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/H_NJ-HCC1395-HCC1395_BL.vcf --reference gs://analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes.fa
```
# Depth Filter

## Input

Our Python script does not accept `gs://` paths. We must stage the file to run locally in our Cloud Shell Terminal.
NOTE: Again, replace `$USER_CLEAN` with the value.
```
gsutil cp gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/H_NJ-HCC1395-HCC1395_BL.vcf .
```

```
depth_filter.py --minimum_depth=30 H_NJ-HCC1395-HCC1395_BL.vcf H_NJ-HCC1395-HCC1395_BL H_NJ-HCC1395-HCC1395_BL.depth_filter.vcf
```

## Output

The file exists within our Docker container in Cloud Shell. Save it to the BFX Workshop scratch bucket.
NOTE: One more time, replace `$USER_CLEAN` with the actual value.
```
gsutil cp H_NJ-HCC1395-HCC1395_BL.depth_filter.vcf gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/
```

When finished with GATK and depth_filter.py, exit the Docker container in Cloud Shell Terminal returning to the original prompt.
```
exit
```

# IGV

View the BAM and depth filtered VCF in IGV using the BFX Workshop Scratch Cloud Bucket.

# Cleanup
From a Cloud Shell Terminal, remove BFX Workshop scratch space:
NOTE: `$USER_CLEAN` should be set in our Cloud Shell Terminal environment.
```
gsutil rm -r gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/
```
