## Requirements

To run a workflow we need three things:

1. Tools: The definitions of what to run.
2. Data: The input files to run on.
3. Executor: A program that coordinates sending data to the tools and collecting the results.

### Tools

There are several languages available for defining processes, including CWL, WDL, and Nextflow.  We have a collection of workflows in CWL available at https://github.com/genome/analysis-workflows/, so we can pull it down with:

```bash
git clone https://github.com/genome/analysis-workflows.git
```

### Data

For a quick demonstration, the git repository we cloned also includes some example data.  For instance, there is a file `example_data/split_interval_list.10.t.yaml`.  This in turn points to the `detect_variants/chr22.interval_list` file.

### Executor

For CWL, the simplest executor is the official `cwltool`.  It can be installed with pip.

```bash
#create a new virtualenv if desired
mkdir ~/venvs/
virtualenv ~/venvs/cwl
source ~/venvs/cwl/bin/activate

#install cwltool
pip install cwltool
```

## Running a small example.

We can then put all the pieces together into an example:

```bash
cd analysis-workflows #where we cloned the repo of CWL

cwltool definitions/tools/split_interval_list.cwl example_data/split_interval_list.10.t.yaml
```

This will pull down the docker image for the tool that splits interval lists and run it on our input interval list.  By default the output will be in our current directory once it completes.

```bash 
ls *.interval_list
```
The result should show 10 files (matching the scatter count in our input YAML file).
```
1.interval_list  10.interval_list 2.interval_list  3.interval_list  4.interval_list  5.interval_list  6.interval_list  7.interval_list  8.interval_list  9.interval_list
```

## A larger example.

### Another executor.

cwltool is good for running locally, but scaling up requires something capable of more parallelization.  The compute1 cluster also disallows running `docker run` directly, which precludes using cwltool there.  Other executors include Cromwell and toil.  We already have a version of cromwell in a docker image that can be used on the compute1 cluster.  To walk through an example of using that executor, see the [alignment walkthrough](./cromwell_alignment_walkthrough.md).

