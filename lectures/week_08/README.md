# Week 08 - Somatic Variant Calling

This week, we cover somatic variant calling and some of the characteristics of tumor data that make variant calling harder. The assignment walks through making calls from a paired set of tumor/normal bams.

- [Lecture recording]()
- [Slides]()


## Assignment

Start by gathering some data. Here is a set of aligned bam files (human build38) similar to those you generated in last week's assignment:

```
#normal bam 
gs://icts-precision-health-bfx-workshop-scratch/inputs/normal.bam

#tumor bam
gs://icts-precision-health-bfx-workshop-scratch/inputs/tumor.bam

```

These contain reads from a few regions of chromosome 17 of a human cell line tumor/normal pair.  

### Mutation calling with Mutect
Today you're going to run two different somatic variant calling algorithms - Mutect and Varscan.  Let's try running the core callers by hand first.

Open up a cloud shell, and make a directory called `week_08` in your home directory. Move into that directory, 

Fire up a docker image containing mutect and being sure to mount the current directory so you have access to your data: `docker run -it -v $(pwd -P):/week_08 broadinstitute/gatk:4.1.8.1`. When it launches, run `pwd -P` to list your current directory and note that we're not in our "week_08" directory any more!  Use `cd` to get to the location where you mounted it.

Once there, start by setting some variables so you don't have to type the full gs:// path later:

```
NORMAL_BAM=gs://icts-precision-health-bfx-workshop-scratch/inputs/normal.bam
```
make sure you've stored it correctly by typing `echo $NORMAL_BAM`, Then do the same for the tumor bam using the path above.

Mutect, as part of the GATK suite, is capable of operating directly on files at these gs:// paths.  So now you can run Mutect like this:

```
/gatk/gatk Mutect2 --java-options "-Xmx1g" -O mutect.vcf.gz -R gs://icts-precision-health-bfx-workshop-scratch/inputs/chr17.fa -I $TUMOR_BAM -tumor Exome_Tumor -I $NORMAL_BAM -normal Exome_Normal
```

Mutect will display a progress monitor - observe that this takes about 15 minutes to run, even though we're only operating on a very small portion of the genome. That's because of some of the complications that we talked about in this week's lecture - it's more complex that germline calling.

Confirm that your output VCF file exists, then exit the docker container (type `exit`).  Back outside of the docker container, try running `echo $NORMAL_BAM`. What happens? Why?

### Mutation calling with Varscan

Varscan, like many other tools, isn't capable of operating directly on google bucket paths. You'll have to pull the files to your local disk before doing calling.  Copy the normal and tumor bams to your current directory using `gsutil cp <PATH> .` replacing <PATH> with the above bucket links. 

Do the same for the index files which are in the same bucket as the bams (e.g. `tumor.bam.bai`)  We're also going to need the reference fasta and index, so grab `chr17.fa` and the accompanying `chr17.fa.fai` and `chr17.dict` file. 

At this point, running `ls -1` should return this:

```$ ls -1
chr17.fa
chr17.fa.fai
mutect.vcf.gz
mutect.vcf.gz.stats
mutect.vcf.gz.tbi
normal.bam
normal.bam.bai
tumor.bam
tumor.bam.bai
```

To actually run Varscan, drop into a container with VarScan installed: `mgibio/varscan_helper-cwl:1.0.0`. Run it with the below command Note the piping on the command line with `<()`. It takes the input of `samtools mpileup` and feeds it to varscan as if it were a file.

```
java -jar /opt/varscan/VarScan.jar somatic <(/opt/samtools/bin/samtools mpileup --no-baq -f chr17.fa normal.bam tumor.bam) varscan.vcf --strand-filter 0 --min-coverage 8 --min-var-freq 0.1 --p-value 0.99 --mpileup 1 --output-vcf
```

Which caller was faster? Go ahead and exit the varscan container.

Both callers typically require post-processing to give only high-quality variant calls, but let's take a quick look at the raw outputs. Mutect outputs a gzipped file, so you'll have to either run `gunzip -c mutect.vcf.gz | less`, or for a shortcut, just use `zless mutect.vcf.gz`

- Which caller produced more variants in the VCF? (how can you count variants without the header?)
- How could we use bash set operations to do a quick and dirty look at number of identical variants called?
<details><summary>hint</summary>
<p>
Example:

```
#extract just the variant information from the VCFs
zgrep -v "^#" mutect.vcf.gz | cut -f 1,2,4,5 >tmp.mut
zgrep -v "^#" varscan.vcf.snp | cut -f 1,2,4,5 >tmp.var
zgrep -v "^#" varscan.vcf.indel | cut -f 1,2,4,5 >>tmp.var
#use bash utilities to see how many appear twice
cat tmp.mut tmp.var | sort | uniq -d | wc -l
#cleanup
rm tmp.mut tmp.var
```

</p>
</details>

### Using a Pipeline

Varscan is a little different because a) it produces separate snp and indel VCFs, which require later merging and b) is very permissive (sensitive), and then relies on a lot of post-caller filtering. Let's use a CWL to streamline some of these steps, specifically the workflow here: [`varscan_pre_and_post_processing.cwl`](https://github.com/genome/analysis-workflows/blob/master/definitions/subworkflows/varscan_pre_and_post_processing.cwl)

- Click around and trace a path through the steps of this workflow to get a quick understand what it's doing under the hood.

Now, let's put together an inputs yaml file for this pipeline. Open a text editor on your local computer (not in the cloud shell) and paste in the following:


```
reference: 
  class: File
  path: gs://icts-precision-health-bfx-workshop-scratch/inputs/chr17.fa
tumor_bam:
  class: File
  path: /path/to/tumor.bam
  secondaryFiles:
    - class: File
      path: /path/to/tumor.bam.bai
    - class: File
      path: /path/to/tumor.bai
normal_bam:
  class: File
  path: /path/to/normal.bam
  secondaryFiles:
    - class: File
      path: /path/to/normal.bam.bai
    - class: File
      path: /path/to/normal.bai
normal_sample_name: Exome_Normal
tumor_sample_name: Exome_Tumor
scatter_count: 2
```

Lots of things to talk about here: 

1) We need to drop in the actual `gs://` paths to the bam files and bam indices from above (replacing `/path/to/...`). Note that we're passing in the bam's index file twice here: once as `file.bam.bai` and once as `file.bai`. Tools tend to be very opinionated about which one they'll accept, and there's no consensus in the field, sadly. In our pipelines, we often just pass the bam index both ways so we don't have to worry about it.

2) Looking at the inputs that CWL expects, we're missing the `interval_list` parameter.  This is nothing complicated - it's just a list of regions of the genome in which to call variants. In exome data, for example, this would be the genic regions that were targeted. For the purposes of this exercise, we'll use this file of just the exome intervals on chr17 - add it as one of your inputs in this file: `gs://icts-precision-health-bfx-workshop-scratch/inputs/chr17_exome.interval_list`

3) These inputs have defaults that will work for our purposes today.  It's important to be aware of them, though!  For example, `min_var_freq` is set to `0.1`. If you're doing deep sequencing and looking for very rare variants, or have a very impure tumor, this might not be an appropriate setting!

4) Another useful parameter here is the "scatter_count".  This pipeline is set up such that it'll split up your data into chunks and run them in parallel on different machines, to speed things up.  In this case, we're just going to use 2, since our data is very small. Our pipelines often use as many as 50 for full Exome or WGS data.

Save the file to your local hard drive as `inputs_week_08.yaml`.

Now we're all set up to run the workflow, just like we did for the alignment exercise.

### Running the workflow

We have a cromwell server already running in the cloud, and you need to connect to it, using your authenticated WUSTL Key credentials. To do so, run the following gcloud command:

```
gcloud compute start-iap-tunnel cromwell-server 8000 --local-host-port=localhost:8080 --zone=us-central1-a --project icts-precision-health
```
Then, open a Web Preview using port 8080 using the button at the top of your cloud shell window. 

Expand the first section `Submit a workflow for execution`, click `Try it out`, then fill in the following values:

1) `workflowUrl`, set to `https://raw.githubusercontent.com/genome/analysis-workflows/master/definitions/subworkflows/varscan_pre_and_post_processing.cwl`

2) `workflowInputs` upload the yaml file you created and saved on your hard drive.

3) `workflowType` select 'CWL'

4) `workflowTypeVersion` select 'v1.0'.

Execute the workflow and **record the ID** returned for use in later steps.

That's going to run for a bit, shooting out jobs to different machines in the cloud and waiting for them to return results.  

### Filtering mutect results
While you're waiting, let's go back to the cloud shell and get on with filtering our mutect results.  Hit "Ctrl-C" in the terminal to stop our connection to the cromwell server - we'll fire it back up and check the status later. 

There are multiple layers of filtering in our somatic pipeline, but here we'll just run the first-pass, built-in filter, using the same GATK container as above (`broadinstitute/gatk:4.1.8.1`).  Hint: You can use the up arrow to scroll back through your commands and reuse the previous one.

```
/gatk/gatk FilterMutectCalls -R gs://icts-precision-health-bfx-workshop-scratch/inputs/chr17.fa -V mutect.vcf.gz -O mutect.filtered.vcf.gz
```

You'll find some detailed stats in `mutect.filtered.vcf.gz.filteringStats.tsv` afterwards.  Now, compare the number of mutations in these VCFs before and after filtering:

```
zgrep -v "^#" mutect.vcf.gz | wc -l
768

zgrep -v "^#" mutect.filtered.vcf.gz | wc -l
768
```
Huh? I thought we were removing false positives here!  The key to understanding this is to look into the VCF at the FILTER column - if you use `zless` to look at your files, you'll note that in the filtered VCF, there are values filled in: `weak_evidence`, `PASS`, etc.  The key to understanding the FILTER field is that any value besides `.` or `PASS` means that the call is being removed/filtered.  It can be useful to keep those around, though, in case we want to "recover" low VAF variants, or understand why a specific site was removed.

Nothing is marked as filtered in the original VCF, so there are 768 variant calls.  Now let's try a more complicated command to figure out how many variants are in the post-filtering version:

```
zgrep -v "^#" mutect.filtered.vcf.gz | cut -f 7 | egrep "PASS|\." | wc -l
```

That's more like it!

### Getting the Varscan results
To finish this section up, let's check back in on that Varscan workflow that was running. Make sure you're out of the gatk docker container, then re-connect to cromwell:

```
gcloud compute start-iap-tunnel cromwell-server 8000 --local-host-port=localhost:8080 --zone=us-central1-a --project icts-precision-health
```
Then, open a Web Preview using port 8080 using the button at the top of your cloud shell window.

Find the `Status` API call and check on your workflow using the workflow ID you copied down above. If it's Running, check back later. If it's "Succeeded", then you're good to go.  Move down to the section that returns outputs, and execute it on your workflow ID.  

It will return a big block of json that contains the outputs from the pipeline.  Copy down the path for the `filtered_vcf`.  Now switch back to your cloud shell window, use `CTRL-C` to kill the cromwell connection, and then use `gsutil cp` to grab that Varscan file to your working directory (`week_08`).  Also grab the `.tbi` file next to it.  This file contains an index of the VCF that lets tools jump quickly to any point in the file (similar to the index next to our bam files).  

- How many variants are contained in this file, after several filtering steps?


### Merging Variants
To complete this exercise, let's merge these two VCFs into one output file. We're going to accomplish this using the CombineVariants tool bundled with GATK.  Drop into this container:  `mgibio/gatk-cwl:3.6.0`

Then run the command:

```
/usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan -o combined.vcf.gz -R chr17.fa --variant:mutect mutect.filtered.vcf.gz --variant:varscan varscan.vcf.snp
```

Notice that we have to tell the tool what to do with overlaps. If both VCFs contain the same variant call, which line do we keep?  Different tools output different types of auxiliary information, and so we have to choose one. It also adds the "set" attribute, which tells which item came from which VCF (or contains a value like "FilteredInAll" which tells us that it didn't pass anywhere.

- Take a look through this file - how many variants (if any) were called by both callers?

The differences in variants called reflect different underlying statistical models and assumptions. In this case, it also reflects that we didn't run through the whole Mutect filtering process as implemented in our pipelines. Still, it's important to realize that different callers have different strengths and weaknesses!
 
### Wrapup
If all of this seemed a little complicated, well, that's because it is!  Somatic variant calling, filtering, and prioritization is a difficult problem, given the challenges of purity, ploidy, and heterogeneity that we talked about in the lecture.  Thankfully, we've put together a pretty solid pipeline in the [analysis-workflows](https://github.com/genome/analysis-workflows) repository that handles a lot of the complexity. You still need to be aware of some of the underlying issues and relevant parameters, as we discussed above, but using a pre-built pipeline can save you a lot of time.
