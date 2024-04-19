# Homework: SingleSampleQc

For this week's homework you will execute a modified and downsampled version of a Workflow Definition Language (WDL, pronounced Widdle) from the NHGRI Genomic Data Science Analysis, Visualization, and Informatics Lab-space or [AnVIL](https://anvilproject.org/), a cloud-based genomic data sharing and analysis platform.

## SingleSampleQC WDL

The original AnVIL SingleSampleQC WDL for use on [Terra](https://anvil.terra.bio/):
[genome/qc-analysis-pipeline](https://github.com/genome/qc-analysis-pipeline)

## Prerequisites

* This walkthrough uses miniwdl.  There are a few different ways to install it [listed on the miniwdl docs](https://miniwdl.readthedocs.io/en/latest/getting_started.html#install-miniwdl).  The simplest one is `pip install miniwdl`.

* Running the workflow also requires the use of Docker containers. One option for running these is [Docker Desktop](https://www.docker.com/products/docker-desktop/).

## Getting Started

If you're not already reading this in your own local copy of the bfx-workshop repository, start by cloning it.

```
git clone https://github.com/genome/bfx-workshop.git
```

Then we can switch into our new checkout including this week's lecture and homework. If you already have a checkout of the repository, the following path may change. Please adjust appropriately to navigate to the correct directory containing the SingleSampleQc.wdl file.

```
cd bfx-workshop/lectures/week_23/SingleSampleQc
```

Before we launch the workflow, let's check the miniwdl configuration.  Run:

```
miniwdl configure
```

This will ask a series of questions.  It's recommended to enable the call and download cache.  Optionally set the download cache to `~/.cache/miniwdl_downloads/` instead of the default location.

This configuration only needs to happen once unless you want to change the values.

## Running the Workflow

In miniwdl there are two ways to specify the inputs.  They can be listed directly on the command-line (like `--input=value`) or a JSON specifying the inputs can be provided.  

The repo for the QC workflow includes an example input JSON that can be run.  To do this:

```
miniwdl run -o out.json -i SingleSampleQc.downsample.json SingleSampleQc.wdl
```

The input JSON refers to files hosted in Google Cloud Storage, so `miniwdl` will start by downloading them to your local machine (to the download cache we configured earlier).  This may take awhile depending on the download speeds.  Then it will start running the steps of the workflow.

The terminal will be updated with the progress of the workflow as it goes along.

By default miniwdl creates a directory in your current directory named with a timestamp, so once the workflow finsihes we can look at what it created:

## Examining the Outputs

```
ls *_*_SingleSampleQc/
```

There are detailed docs about the contents of this directory in [miniwdl's runner reference](https://miniwdl.readthedocs.io/en/latest/runner_reference.html#i-o-and-run-directory-structure).  If the workflow succeeded, the important thing inside will be the `out/` directory, which has directories for each output of the workflow that link directly to the output files from the workflow.  These will be the same links that are listed in the `out.json` file.  The `out.json` will also include the non-file outputs.

If something went wrong, the `out.json` file should instead contain an error message and the other files in this workflow's directory may be useful for debugging.

## MultiQC

For usage and install instructions, please see [MultiQC](https://github.com/MultiQC/MultiQC).

Generate an HTML report for the workflow output and submit as homework.


