# Week 20: Genomic Workflows/Cloud Computing

- [Lecture Recording](https://wustl.box.com/s/hcgh9shk4a9dddhqxcbcv1i0qi3mgjc7)
- [Lecture Slides](GenomicWorkflows_BFX-Workshop_v2024-2025.pdf)

## Exercise

For this week's exercise, you'll be running a small section of a larger analysis workflow within Google's cloud platform (GCP). You have free access to a small cloud machine using your wustl key login.

Open [google cloud shell](https://shell.cloud.google.com/). Use your @wustl.edu email as the login, which will take you to the wustl SSO login page. Afterwards, you will be taken to the cloud shell main page. The rest of the exercise will be conducted in the terminal at the bottom of the screen, so you may want to resize it. 

### Setup

1. Make a working directory with `mkdir NAMEHERE`, then `cd` into it.
2. Download [Cromwell](https://github.com/broadinstitute/cromwell) using this command: `wget https://github.com/broadinstitute/cromwell/releases/download/79/cromwell-79.jar`. Note that this is not the latest version of Cromwell, but the last release that supported CWL, which is the language our workflow is written in.
3. Download the input data for this workflow using `docker run -v $PWD:/staging mgibio/data_downloader:0.1.0 gsutil -m cp -r gs://analysis-workflows-example-data /staging`. This may take a few minutes. After it completes, there should be a new directory called `analysis-workflows-example-data` within your working directory. This command downloads the new directory as the root user, which may cause some permissions issues later, so let's change ownership to your user: `sudo chown -R $USER:$USER analysis-workflows-example-data/`
4. Download the github repository containing the workflow we will be running: `git clone https://github.com/genome/analysis-workflows.git`. Again, this make take a few minutes, and after it completes, you should have a new directory `analysis-workflows` within your working directory. 
5. Download a file mapping paths to the inputs downloaded in step 3 to the names expected by the workflow we will be using from step 4: `wget https://raw.githubusercontent.com/genome/bfx-workshop/refs/heads/master/lectures/week_20/alignment_exome.yaml`

### Running a workflow
Now that we've downloaded everything we need, we can launch a small workflow (alignment_exome.cwl), which will align reads from a tumor sample to a human reference genome and then perform some basic qc metrics. This is part of a larger workflow (somatic_exome.cwl) that also aligns a normal sample, runs variant detection workflows on both, and then reports variants only found in the normal sample, likely germline variants, and those only found in the tumor sample, likely somatic variants, and in this case, possibly responsible for cancer. This in turn is part of an even larger workflow (immuno.cwl) for generating personalized cancer vaccines targeting tumor-only variants. Due to time and hardware limitations, we will focus on the alignment and QC workflow as a demonstration, but the data is available to launch those larger workflows if you wish. 

The workflow is written in CWL, and we will use Cromwell to launch the workflow and manage inputs and outputs. Cromwell writes output to `stdout`, and while this is useful to monitor the progress of the workflow, we'd also like to save this to a file to make it easier to pull information from later. To do this, we will use the unix `tee` command, which will write to both a specified file and to stdout. From your working directory: `java -jar cromwell-79.jar run -t cwl -i alignment_exome.yaml analysis-workflows/definitions/pipelines/alignment_exome.cwl | tee cromwell.log`. This will take about 10 minutes to run, and you'll see a lot of output in your terminal during this time.

### Examining the outputs
Open the file `cromwell.log` with an editor such as vim**. Near the end of the file is a large chunk of JSON specifying the names and paths to all outputs. These names match the names in the CWL file we ran. Most of these files are reporting QC metrics about the alignment. By tracing the inputs and outputs in the CWL file to the other files it calls, you can trace each output to the tool and the command template used to generate it. The commands are also present in the log file, although they can be hard to pull out of the surrounding text.

** If you are not familiar with vim, it can be quite confusing. Please consult a [guide like this one](https://www.geeksforgeeks.org/basic-vim-commands/) to learn the basics, or use a more user friendly editor like nano. Alternatively, if you'd prefer to examine the log locally, see the Homework section below for instructions on downloading files from cloud shell to your local machine.

## Homework
As proof of completion, find the file `alignment_exome.cwl.insert_size_histogram` in the Cromwell log. Copy the file path, then click the 3 dots at the top right of your window, click download, and paste the path into the window that pops up. Note that the Cromwell log lists the full path, and the download window automatically fills in a prefix corresponding to your home directory, so either delete that and then paste the full path from Cromwell, or trim the prefix from the Cromwell log path before pasting. Email this pdf to John.
