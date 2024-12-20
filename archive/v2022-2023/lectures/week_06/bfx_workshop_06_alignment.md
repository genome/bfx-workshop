## Week 06 Alignment Exercise

Today, we're going to improve upon last week's alignment workflow in a few ways.

1. We're going to add a few quality control steps.
2. We're going to add automation, via a CWL pipeline of tools that is run using a service manager called Cromwell.
3. We're going to run it in the cloud instead of on your laptop.
 
# Requirements

To run this workflow we need three things:

1. Tools: The definitions of what to run.
2. Data: The input files to run on.
3. Executor: A program that coordinates the transfer of data, schedules the execution of steps, and collects the results

## Tools

There are several languages available for defining workflows, including [Common Workflow Language](https://www.commonwl.org/) (CWL), [Workflow Description Language](https://openwdl.org) (WDL), [SnakeMake](https://snakemake.readthedocs.io/en/stable/), and [Nextflow](https://www.nextflow.io/) (NF).  

We have a collection of workflows in CWL available at https://github.com/genome/analysis-workflows/

We also have a collection of workflows in WDL available at https://github.com/griffithlab/analysis-wdls

Today we will focus on the following CWL workflow for doing alignment of DNA sequences:
https://github.com/genome/analysis-workflows/blob/2022-bfx-workshop/definitions/pipelines/alignment_exome.cwl

Take some time to review the structure of the CWL workflow, including the actual alignment commands called. Think about how many commands this is producing, and how long it would take you to run them by hand (and how that might not scale to hundreds or thousands of samples!)

## Data

The input files are defined in a YAML for each sample.

```
HCC1395 Normal
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/normal_alignment_exome_gcloud.yaml

HCC1395 Tumor
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/tumor_alignment_exome_gcloud.yaml
```

Download an input file, open it up, and take a closer look.  Notice that all the file paths begin with `gs://` which means they're already stored in a Google Cloud bucket. If you were using your own data, you'd have to upload it first.

## Executor

[Cromwell](https://cromwell.readthedocs.io/en/stable/) is the Workflow Mangement System we will use today. It takes in a description of a workflow, and handles all the behind-the-scenes actions that request appropriate machine sizes, spin up docker containers, execute the code, and keep track of hte outputs.  Specfically we will use Cromwell in Server mode: https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/

Cromwell is an implementaiton of the [Global Alliance for Genomics and Health](https://www.ga4gh.org/) (GA4GH) [Workflow Execution Service](https://ga4gh.github.io/workflow-execution-service-schemas/docs/) (WES) Application Programing Interface (API) standard.

# Tutorial

## Install gcloud

Please download and install gcloud on your local machine. We will not use it to execute the workflow. However, it will be used later to access workflow outputs, logs, etc. Using gcloud, you can use your local terminal for many tasks. For Google installation instructions:
https://cloud.google.com/sdk/docs/install

### gcloud Install Hints:
1. After installation, gcloud may not automatically appear in your $PATH. If so, `source /Users/<username>/google-cloud-sdk/path.bash.inc`
2. Similarly, auto-completion may not work immediately. If so, try `source /Users/<username>/google-cloud-sdk/path.bash.inc`
3. If you are authenticating with gcloud for the first time, try `gcloud auth login` which will open a browser window prompting you to login. Use your WUSTL key to get access to the course workspace and resources.

## Cloud Shell

Login to a Google Cloud Shell using your WUSTL Key: https://shell.cloud.google.com/?show=ide%2Cterminal  This is just like the terminal on your computer, but it's running on one of Google's servers in the cloud.  Be sure that you're logged in using your WUSTL key and not your personal account. A separate browser (firefox/chrome/safari) may be helpful for this, or an Incognito/private window.

## Cromwell Server

### Communications Tunnel
 
We have a cromwell server already running in the cloud, and you need connect to it, using your authenticated WUSTL Key credentials. To do so, run the following gcloud command:
```
gcloud compute start-iap-tunnel cromwell-server 8000 --local-host-port=localhost:8080 --zone=us-central1-a --project icts-precision-health
```

This command opens a communication (IAP for TCP forwarding) channel between the Cromwell Server and your Cloud Shell or local (if using gcloud locally) terminal session.  In other words, it allows the shell you're sitting at to talk to the cromwell server which will run your workflow.

### Cromwell API

Using Cloud Shell, open a Web Preview using port 8080. Web Preview is a small icon near the top right of the Cloud Shell interface.  If you later want to try this using a local terminal session with gcloud (instead of the cloud shell), you would load this URL in a web browser of your choice: http://localhost:8080

Once a Web Preview (via Cloud Shell) or web browser (via local browser) loads, a Swagger User Interface (UI) to the Cromwell API should be available.

## Submit Workflow

Before cromwell can run our workflow, we need to use the POST command to send data to the server. The first line should say `Submit a workflow for execution`.  Click to open it and then click `Try it out`.  Then we need to provide a few key pieces of information: 

1) Where are the files describing the workflow you want to run? (here, a DNA alignment pipeline)

For `workflowUrl`, use
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/definitions/pipelines/alignment_exome.cwl
```

2) Where are the input files for this workflow? (what data are we aligning?)

For `workflowInputs`, you will need to download the following YAML file and then 'Choose File' to upload it to the server.
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/normal_alignment_exome_gcloud.yaml
```

3) What kind of workflow language files are you feeding it?  

For `workflowType` select 'CWL' and for `workflowTypeVersion` sekect 'v1.0'. 
NOTE: Do NOT use `1.0`! There is a difference and `1.0` would be used for `WDL` workflows.

All other settings (if not listed) can be left as the default setting or without a value.

`Execute` the workflow and record the ID returned for use in later steps.

#### Understanding the API output

The output of these executions can be confusing. The first two lines: `Curl` and `Request URL` show you information about how the API call was run. This is generally not going to be that useful to you.  What you're interested in is the `Server Response` section, which contains the actual information that the server sent back to you. Speciically, the `Response Body` will generally have the info you need - in this case, a workflow id like `b15cbfc3-f697-4d4a-8b18-8c4fbcad357d`.  Below that, in the `Responses` block, is information on the possible responses that could have been sent back. These can also generally be ignored.

## Workflow API Calls

Now we have a workflow that is running in the background.  How can we see what's going on with it, or access the data when it's done?  We can use additional API calls to the server to answer those questions.

### Workflow Status

Again, find this option in the list of API calls, click `Try it out`, and then paste in your workflow ID (from the above submission) to `GET` the workflow status.  You should see your workflow as `Running`, and then after about 20-30m, if you check back, it will say `Succeeded`.  If it shows anything else, check out "Troubleshooting" below.

### Workflow Outputs 

`GET` the workflow outputs via the API. This returns a JSON object containing the output paths from each workflow step configured to output the results as part of the workflow. Again, you may use `gsutil` to view or download these files.  

1. Find the CollectAlignmentSummaryMetrics file and look it over.  What percentage of reads were aligned?  What percentage were aligned in pairs?  What could low alignment rates mean for your experiment?
2. Find the PDF containing the Insert Size Histogram, copy it to your filesystem and open it.  What could a very small insert size tell you about your DNA or library preparation?
3. Find the path to the output bam file. Open up IGV, and use `File > Load From URL` to view the bam directly through the gs:// path.  Browse to the BRCA1 gene region and zoom in to an expon.  Notice that IGV is capable of "streaming" bam files.  This works because it doesn't need to load the whole thing, but can instead use the bam index to just grab the chunk of sequence data that is needed for display.  This is really useful when operating on whole genome bams which can be over 100Gb each!

### Viewing the outputs in your browser 
By going to https://console.cloud.google.com/storage/browser/icts-precision-health-cromwell-wf-exec/alignment_exome.cwl you can see a GUI of the bucket that we're all writing to. Click on the folder that matches your workflow ID.  

Cromwell creates a complicated directory structure that follows the nested path of the CWL. Browse through the workflow's alignment portion until you find the mark_duplicates_and_sort output.  In that folder should be a `mark_duplicates_and_sort.log` file. Click on it and you'll see options for accessing it, including an `https://` link that will open in your browser, or a `gs://` path that will allow you to use `gsutil` to cat or copy it. 

Look at the log and compare it to what the outputs look like when you ran those tools locally.  There is extra info about docker starting up, etc, but the same STDERR output from the tool being run is there.  This can be extremely useful for troubleshooting why a certain step in a workflow failed.


## Repeat
If the Normal HCC1395 sample succeeds, repeat the tutorial using the Tumor HCC1395 sample inputs:
```
https://raw.githubusercontent.com/genome/analysis-workflows/2022-bfx-workshop/example_data/tumor_alignment_exome_gcloud.yaml
```
Next week, we'll start using these results to do variant calling.

## Troubleshooting

- My cloud shell timed out and gave me a connection error!

Just open a new one and run the same commands to get to your swagger API page. The server that's actually doing the work should still be chugging around happily in the background 

- My API call returns `404 Error:  Not Found` and `content-length 0 . . .`

This likely means that your cloud shell has timed out. Go back to that page, reconnect, run the connection command and reopen your swagger API window.

- My Execute command returned `{"status": "fail",  "message": . . .`

Read the message, see if you can make sense of what happened (bad input files?) and if not, post in slack to get help!

- When I checked the status of my workflow, I got a `501 Unsupported HTTP method` error.  

If you just executed your workflow, wait a minute or two and try again. It takes a short time for the server to spin up the job.

