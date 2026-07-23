# Week 23: Genomic Workflows/Cloud Computing

- [Lecture Recording](https://wustl.box.com/s/qwmq2kh3xdfbaez1aunb594uzl29e3gz)
- [Lecture Slides](bfx_genomic_workflows.pdf)


## Exercise

For this week's exercise, you'll be running a small section of a larger analysis workflow within Google's cloud platform (GCP). You have free access to a small cloud machine using your wustl key login.

1. Open [google cloud shell](https://shell.cloud.google.com/).
2. Use your @wustl.edu email as the login, which will take you to the wustl SSO login page.
3. Afterwards, you will be taken to the cloud shell main page. 

The rest of the exercise will be conducted in the terminal at the bottom of the screen, so you may want to resize it. 

### Setup

1. Make a working directory with `mkdir workflow`, then `cd` into it.

2. Download [Cromwell](https://github.com/broadinstitute/cromwell) using this command: `wget https://github.com/broadinstitute/cromwell/releases/download/79/cromwell-79.jar`. Note that this is not the latest version of Cromwell, but the last release that supported CWL, which is the language our workflow is written in.

3. Download the input data for this workflow using `docker run -v $PWD:/staging mgibio/data_downloader:0.1.0 gsutil -m cp -r gs://analysis-workflows-example-data /staging/`. This may take a few minutes. 

4. After it completes, there should be a new directory called `analysis-workflows-example-data` within your working directory. This command downloads the new directory as the root user, which may cause some permissions issues later, so let's change ownership to your user: `sudo chown -R $USER:$USER analysis-workflows-example-data/`

5. Look at what you've downloaded using `cd`, `ls` , and `less`.  There are a few types of files:  
	a. YAML files (like `rnaseq.yaml`), which list inputs to the pipeline
	b. sample data files to use in these workflows, like `unaligned_subset_bams/tumor/2895499237.bam` containing raw sequence reads, or `somatic_inputs/hla_and_brca_genes.fa`, containing a tiny reference genome we can align against.

6. Move back to the main directory: `cd ~/workflow`, then download the github repository containing many genomic workflows written in CWL: `git clone https://github.com/genome/analysis-workflows.git`. Again, this make take a few minutes, and after it completes, you should have a new directory `analysis-workflows` within your working directory. 

7. Download a file mapping paths to the inputs downloaded in step 3 to the names expected by the workflow we will be running below: `wget https://raw.githubusercontent.com/genome/bfx-workshop/refs/heads/master/archive/v2024-2025/lectures/week_20/alignment_exome.yaml`.  Look at the inputs in the folder, some of which will look familiar from our alignment exercise way back in week 04. 

### Running a workflow
Now that we've downloaded everything we need, we can launch a small workflow (`alignment_exome.cwl`), which will perform the dozens of steps needed to align reads from a tumor sample to a human reference genome and then generate basic qc metrics.

This is part of a larger workflow (`somatic_exome.cwl`) that also aligns a normal sample, runs variant detection workflows on both, and then reports  germline variants, somatic variants, and many kinds of additional files (annotations, qc, etc). 

The somatic workflow, in turn, is part of an even larger workflow (`immuno.cwl`) for generating personalized cancer vaccines targeting somatic variants. All this is to say, workflows can be nested, and can by used to reproducibly produce some very complex outputs.

Due to time and hardware limitations, we will focus only on the small alignment and QC workflow as a demonstration, but the data is available to launch those larger workflows if you wish. 

The workflow is written in CWL, and we will use Cromwell to launch the workflow and manage inputs and outputs. Cromwell writes output to `stdout`, and while this is useful to monitor the progress of the workflow, we'd also like to save this to a file to make it easier to pull information from later. To do this, we will use the unix `tee` command, which will write to both a specified file and to stdout. 

From your `workflow` directory, run: 

```
java -jar cromwell-79.jar run -t cwl -i alignment_exome.yaml analysis-workflows/definitions/pipelines/alignment_exome.cwl | tee cromwell.log`. 
```

This will take about 10 minutes to run, and you'll see a lot of output in your terminal during this time.  It is _very_ verbose, but all this detail can be useful when debugging issues with the pipeline.

### Examining the outputs
Open the file `cromwell.log` with `less`. Near the end of the file is a large chunk of structured code in JSON format that looks like this:

```
[2026-04-13 04:25:22,15] [info] SingleWorkflowRunnerActor workflow finished with status 'Succeeded'.
{
  "outputs": {
    "alignment_exome.cwl.per_base_hs_metrics": [{
      "format": null,
      "location": "/home/c_a_miller/workflow/cromwell-executions/alignment_exome.cwl/a928e4cb-7f5e-442a-9828-f6ad5670e7f4/call-qc/qc_exome.cwl/6822ac90-20c9-41cb-b1ce-3c3790f9427e/call-collect_detailed_hs_metrics/hs_metrics.cwl/2fbc6568-5354-4632-836c-ab799d0e9f5e/call-collect_per_base_hs_metrics/shard-0/execution/H_NJ-HCC1395-HCC1395.base-clinvar-HsMetrics.txt",
      "contents": null,
      "checksum": null,
      "class": "File",
      "size": 5600,
      "secondaryFiles": []
    }],
    "alignment_exome.cwl.hs_metrics": {
      "format": null,
      "location": "/home/c_a_miller/workflow/cromwell-executions/alignment_exome.cwl/a928e4cb-7f5e-442a-9828-f6ad5670e7f4/call-qc/qc_exome.cwl/6822ac90-20c9-41cb-b1ce-3c3790f9427e/call-collect_roi_hs_metrics/execution/H_NJ-HCC1395-HCC1395.roi-HsMetrics.txt",
      "contents": null,
      "checksum": null,
      "class": "File",
      "size": 5031,
      "secondaryFiles": []
    },
    
    
. . .     
```

This specifies the names and paths to all outputs defined in the workflow CWL file we ran. In this run, most of the files are reporting QC metrics about the alignment. By tracing the inputs and outputs in the CWL file to the other files it calls, you can trace each output to the tool and the command template used to generate it. The commands are also present in the log file, although they can be hard to pull out of the surrounding text.


## Homework
As proof of completion, find the output file `alignment_exome.cwl.insert_size_histogram` in the Cromwell log. Copy the file path, then click the 3 dots at the top right of your window, click download, and paste the path into the window that pops up. Note that the Cromwell log lists the full path, and the download window automatically fills in a prefix corresponding to your home directory, so either delete that and then paste the full path from Cromwell, or trim the prefix from the Cromwell log path before pasting. Email this pdf to John.
