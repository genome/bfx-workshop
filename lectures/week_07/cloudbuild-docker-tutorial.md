1. Cloud Shell Workspace

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

Set a new "clean" shell variable to remove invalid characters:
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
NOTE: Use the Cloud Shell feature to switch between the Terminal (using `Open Terminal`) and the Dditor (using `Open Editor`) where labeled throughtout the remainder of this tutorial.

In a Cloud Shell Terminal, navigate into the newly created directory:
```
cd gatk-depth-filter-docker
``` 

2. Python Script

Using the Cloud Shell Editor, Save a file named `depth_filter.py` in the `gatk-depth-filter-docker` workspace. 

Copy/Paste the contents of `depth_filter.py` from GitHub:
https://raw.githubusercontent.com/genome/docker-depth-filter/master/depth_filter.py

In a Cloud Shell Terminal, add executable permissions to the Python script:
```
chmod +x depth_filter.py
```

3. Dockerfile

Using the Cloud Shell Editor, Save a file named `Dockerfile` in the `gatk-depth-filter-docker` workspace containing:
```
FROM broadinstitute/gatk:4.3.0.0
RUN apt-get update && apt-get install -y libnss-sss && apt-get clean all
RUN pip install vcfpy pysam
COPY depth_filter.py /usr/bin/depth_filter.py
```

4. Cloud Build

Using the Cloud Shell Editor, Save a file named `cloudbuild.yaml` in the `gatk-depth-filter-docker` workspace containing:
NOTE: REPLACE `$USER_CLEAN` WITH THE VALUE FROM THE ENVIRONMENT VARIALBE, ie. your WUSTL Key with `-` instead of `_` characters.
```
steps:
- name: 'gcr.io/cloud-builders/docker'
  args: [ 'build', '-t', 'us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:tag1', '.' ]
images:
- 'us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:tag1'
```

In a Cloud Shell Terminal, submit the build:
```
gcloud builds submit --region=us-central1 --config cloudbuild.yaml
```

Use Google Cloud Build to view each build: https://console.cloud.google.com/cloud-build/dashboard;region=us-central1?project=icts-precision-health

5. Docker Run

First, from a Cloud Shell Terminal, we need to configure the Artifact Registry credentials for the region Artifacy Registry we intend to pull the Docker image from.
```
gcloud auth configure-docker us-central1-docker.pkg.dev
```

Now, we can use the Docker `run` command to pulll the image from the Artifact Registry and jump into the container once it's running with Docker.
```
docker run -it us-central1-docker.pkg.dev/icts-precision-health/bfx-workshop-repo/$USER_CLEAN-gatk-depth-filter-image:tag1 /bin/bash
```

# GATK
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

Replace `$USER_CLEAN` in the path below

Replace `$BAM` with the path returned by cromwell-server. 
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
NOTE: One more time, replace `$USER_CLEAN` with the value.
```
gsutil cp H_NJ-HCC1395-HCC1395_BL.depth_filter.vcf gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/
```

When finished with GATK and depth_filter.py, exit the Docker container in Cloud Shell Terminal returning to the original prompt.
```
exit
```


6. IGV

View the BAM and depth filtered VCF in IGV.

7. Cleanup

Remove BFX Workshop scratch space:
NOTE: Replace `$USER_CLEAN` with the value.
```
gsutil rm -r gs://icts-precision-health-bfx-workshop-scratch/$USER_CLEAN/
```
