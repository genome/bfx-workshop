# Requirements

To run a workflow we need three things:

1. Tools: The definitions of what to run.
2. Data: The input files to run on.
3. Executor: A program that coordinates the transfer of data, schedules the execution of steps, and collects the results

# Tools

There are several languages available for defining workflows, including [Common Workflow Language](https://www.commonwl.org/) (CWL), [Workflow Description Language](https://openwdl.org) (WDL), [SnakeMake](https://snakemake.readthedocs.io/en/stable/), and [Nextflow](https://www.nextflow.io/) (NF).  

We have a collection of workflows in CWL available at https://github.com/genome/analysis-workflows/

We also have a collection of workflows in WDL available at https://github.com/griffithlab/analysis-wdls

Today we will focus on the following CWL workflow:
https://github.com/genome/analysis-workflows/blob/2022-bfx-workshop/definitions/pipelines/alignment_exome.cwl

Let's take sometime to review the structure of the CWL workflow including the actual alignment commands called.

# Data

The input files are defined in a YAML for each sample.

HCC1395 Normal
https://github.com/genome/analysis-workflows/blob/2022-bfx-workshop/example_data/normal_alignment_exome_gcloud.yaml

HCC1395 Tumor
https://github.com/genome/analysis-workflows/blob/2022-bfx-workshop/example_data/tumor_alignment_exome_gcloud.yaml

Let's take a closer look at each input file.

# Executor

[Cromwell](https://cromwell.readthedocs.io/en/stable/) is the Workflow Mangement System we will use today. Specfically we will use Cromwell in Server mode: https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/

Cromwell is an implementaiton of the [Global Alliance for Genomics and Health](https://www.ga4gh.org/) (GA4GH) [Workflow Execution Service](https://ga4gh.github.io/workflow-execution-service-schemas/docs/) (WES) Application Programing Interface (API) standard.

# Tutorial

## [OPTIONAL] Install gcloud
If you choose to download gcloud, you can use your local terminal rather than the Google Cloud Shell used below.
https://cloud.google.com/sdk/docs/install

## Cloud Shell

Please login to a Google Cloud Shell using your WUSTL Key for authentication: https://shell.cloud.google.com/?show=ide%2Cterminal

## Cromwell Server

### Communications Tunnel
 
In order to access the Cromwell Server with your authenticated WUSTL Key credentials, we need to run the following gcloud command:
```
gcloud compute start-iap-tunnel cromwell-server 8000 --local-host-port=localhost:8080 --zone=us-central1-a --project icts-precision-health
```

The above command opens a communication (IAP for TCP forwarding) channel between the Cromwell Server and your terminal session.

### Cromwell API

Using Cloud Shell, open a Web Preview using port 8080. If you are using a local terminal session with gcloud, simply load http://localhost:8080 in a browser.

You should see a Swagger User Interface (UI) to the Cromwell API commands available.

## Submit Workflow

Select the POST command to submit a workflow to the Cromwell Server.

For `workflowUrl`, please use:
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/definitions/pipelines/alignment_exome.cwl
```

For `workflowInputs`, you will need to download or copy/paste the following YAML file and then 'Choose File' to upload the contents of that file.
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/normal_alignment_exome_gcloud.yaml
```

For `workflowType` select 'CWL' and for `workflowTypeVersion` sekect 'v1.0'

`Execute` the workflow and record the ID returned for use in later steps.

Check the workflow Status using the ID from the submission.

TODO - Add info on checking status, outputs, etc.

Repeat for the Tumor inputs:
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/tumor_alignment_exome_gcloud.yaml
```

