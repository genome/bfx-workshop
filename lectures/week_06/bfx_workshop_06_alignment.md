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

## Install gcloud

Please download and install gcloud. It is not required to execute the workflow. However, it will be used later to access workflow outputs, logs, etc. Using gcloud, you can use your local terminal rather than the Google Cloud Shell example below. For Google installation instructions:
https://cloud.google.com/sdk/docs/install

### glcoud Install HINTS:
1. After installation, gcloud may not automatically appear in your $PATH. If so, `source /Users/<username>/google-cloud-sdk/path.bash.inc`
2. Similarly, auto-compltion may not work immediately. If so, try `source /Users/<username>/google-cloud-sdk/path.bash.inc`
3. If you are authenticating with gcloud for the first time, try `gcloud auth login` which will open a browser window prompting you to login.

## Cloud Shell

Login to a Google Cloud Shell using your WUSTL Key: https://shell.cloud.google.com/?show=ide%2Cterminal

## Cromwell Server

### Communications Tunnel
 
In order to access the Cromwell Server with your authenticated WUSTL Key credentials, we need to run the following gcloud command:
```
gcloud compute start-iap-tunnel cromwell-server 8000 --local-host-port=localhost:8080 --zone=us-central1-a --project icts-precision-health
```

The above command opens a communication (IAP for TCP forwarding) channel between the Cromwell Server and your Cloud Shell or local (if using gcloud locally) terminal session.

### Cromwell API

Using Cloud Shell, open a Web Preview using port 8080. Web Preview is a small icon near the top right of the Cloud Shell interface.

If you are using a local terminal session with gcloud, simply load this URL in a web browser of your choice: http://localhost:8080

Once a Web Preview (via Cloud Shell) or web browser (via local browser) loads, a Swagger User Interface (UI) to the Cromwell API should be available.

## Submit Workflow

Select the POST command to submit a workflow to the Cromwell Server. Please use the following setting when submitting the workflow. All other settings (if not listed) can be left as the default setting or without a value.

For `workflowUrl`, please use:
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/definitions/pipelines/alignment_exome.cwl
```

For `workflowInputs`, you will need to download or copy/paste the following YAML file and then 'Choose File' to upload the contents of that file.
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/normal_alignment_exome_gcloud.yaml
```

For `workflowType` select 'CWL' and for `workflowTypeVersion` sekect 'v1.0'. 
NOTE: Do NOT use `1.0`! There is a difference and `1.0` would be used for `WDL` workflows.

`Execute` the workflow and record the ID returned for use in later steps.
NOTE: Please grab the Workflow ID that is unique to your execution of the workflow. The Sample ID is a different value located near the Workflow ID. Be sure to grab the Workflow ID.

## Workflow API Calls

### Workflow Status

`GET` the workflow status using the Workflow ID from the submission above. Submit the API call and compare the returned value to the list of exmample values in the UI.

### Workflow Metadata

`GET` workflow metadata via an API call. This call to the Cromwell Server API  returns workfow and call-level metadata including logs from the Cromwell Server itself. 
NOTE: The workflow metadata log is useful when the status API returns an ambiguous `501`. The workflow could still be processing or it failed to validate the inputs. The metadata log should provide you additional info on the status.

### Workflow Logs

`GET` the workflow logs via the API. This call returns a Google Cloud path to the *workflow* log that can be examined further with `gsutil`. Once the workflow is indeed running, attempt listing the path with `gsutil ls`. If the log exists, the `gsutil` command will return without error. To see it's contents, try `gsutil cat`.

### Workflow Outputs 

`GET` the workflow outputs via the API. A JSON object containing the output paths from each workflow step configured to output the results as part of the workflow. Again, you may use `gsutil` to interogate these files. For files that you would like to load locally, ex. binary files such as BAM, PDF, PNG, etc., try `gsutil cp` providing a path to the local storage location, ex. Downloads (CWD is the default). 

## Repeat
If the Normal HCC1395 sample succeeds this workflow (~30 minutes or less. You can repeat the tutorial using the Tumor HCC1395 sample inputs:
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/tumor_alignment_exome_gcloud.yaml
```
