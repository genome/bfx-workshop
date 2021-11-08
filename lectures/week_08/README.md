# Week 08 - Somatic Variant Calling

This week, we cover somatic variant calling and some of the characteristics of tumor data that make variant calling harder. The assignment walks through making calls from a paired set of tumor/normal bams.

- [Lecture recording](https://wustl.box.com/s/yizndcjs9kgutm3yg9pj7su2dmz0vv1d)
- [Slides](https://docs.google.com/presentation/d/1kn6Z_j-2a3WW4Au6FVYPUxG7YiT0UEcIHujgWUZCRCQ/edit?usp=sharing)


## Assignment

Start by gathering some data. Pull down a set of aligned bam files (human build38) like so:

```
wget https://www.dropbox.com/s/jhru5qu4bxcmfwq/normal.bam https://www.dropbox.com/s/4gimazw35rcxr8x/normal.bam.bai https://www.dropbox.com/s/uxk7duij016qbcq/normal.bai https://www.dropbox.com/s/0fwwcm92fry66gu/tumor.bam https://www.dropbox.com/s/ihxeqy9cnbiwdn3/tumor.bam.bai https://www.dropbox.com/s/v5h222yxpos169w/tumor.bai
```
These contain reads from chromosome 17 of human cell line tumor/normal pair.  You might wonder why we have both `tumor.bai` and `tumor.bam.bai` files here. The short answer is that they contain the same information. Some tools strongly prefer one format, and some the other, so to avoid complications in our pipelines, we just keep both around.

### Mutation Calling
We're going to run two different somatic variant calling algorithms today - Mutect and Varscan.  Let's try running the core callers by hand first.

Fire up image `broadinstitute/gatk:4.1.8.1` and give it 32G of RAM.  Once inside, run Mutect like this:

```
/gatk/gatk Mutect2 --java-options "-Xmx28g" -O mutect.vcf.gz -R /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa -I tumor.bam -tumor "Exome_Tumor" -I normal.bam -normal "Exome_Normal"
```
For Varscan, drop into a container with VarScan installed: `mgibio/varscan_helper-cwl:1.0.0`  16G of RAM should be fine. Run it like this:

```
java -jar /opt/varscan/VarScan.jar somatic <(/opt/samtools/bin/samtools mpileup --no-baq -f /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa normal.bam tumor.bam) varscan.vcf --strand-filter 0 --min-coverage 8 --min-var-freq 0.1 --p-value 0.99 --mpileup 1 --output-vcf
```
Note the piping on the command line with `<()` that takes the input of `samtools mpileup` and feeds it to varscan.

- Which caller was faster?

Both callers typically require post-processing to give only high-quality variant calls, but let's take a quick look at the raw outputs. Mutect outputs a gzipped files, so you'll have to either run `gunzip -c mutect.vcf.gz | less`, or for a shortcut, just use `zless mutect.vcf.gz`

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

- Trace a path through the steps of this workflow to understand what it's doing under the hood.

To run it, first check out the analysis-workflows repository, if you don't already have it

``` 
git clone git@github.com:genome/analysis-workflows.git 
```

Now, let's put together an inputs file for our CWL.  Lots of these inputs have defaults that will work for our purposes today.  It's important to be aware of them, though!  For example, `min_var_freq` is set to `0.1`. If you're doing deep sequencing and looking for very rare variants, or have a very impure tumor, this would not be an appropriate setting!

Create a file called `inputs.yaml` and fill it with the following information (substituting paths where necessary):

```
reference: /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa
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

Uh oh - looking at the inputs that CWL expects, we're missing the `interval_list` parameter.  This is nothing complicated - it's just a list of regions of the genome in which to call variants. In exome data, for example, this would be the genic regions that were targeted. For the purposes of this exercise, we'll use this exome interval list - add it to your inputs.yaml: `/storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/interval-list/c83be453937b44e7aa0d7e36f6806710/94b23d66f667472690a7ee165e2037b6.interval_list`

Another useful parameter here is the "scatter_count".  This pipeline is set up such that it'll split up your data into chunks and run them in parallel, to speed things up.  In this case, we're just going to use 2, since our data is very small. Our pipelines often use as many as 50 for WGS data.

One last thing we need for our workflow is the cromwell configuration that tells it about our environment and where to stick the results.  Make a copy of the one that you created in the [week 05 homework](https://github.com/genome/bfx-workshop/blob/master/archive/lectures/week_05/cromwell_alignment_walkthrough.md), changing the paths as need. 

Okay, now we're set up to actually run.  Fire up a container that contains cromwell (our workflow manager), giving it 4 gigs of memory: 
``` 
gsub -m 4
```

Then launch the job, pointing to the key three inputs: the cromwell config, the cwl workflow, and the inputs file that we just created:

```
/usr/bin/java -Dconfig.file=cromwell.config -jar /opt/cromwell.jar run -t cwl -i inputs.yaml analysis-workflows/definitions/subworkflows/varscan_pre_and_post_processing.cwl
```

That's going to run for a bit, shooting out jobs to the cluster and waiting for them to come back.  

While you're waiting, let's open a new window and get on with filtering our mutect results. There are multiple layers of filtering in our somatic pipeline, but here we'll just run the first-pass, built-in filter, using the same GATK container as above:

```
/gatk/gatk FilterMutectCalls -R /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa -V mutect.vcf.gz -O mutect.filtered.vcf.gz
```

You'll find some detailed stats in `mutect.filtered.vcf.gz.filteringStats.tsv` afterwards.  Now, compare the number of mutations in these VCFs before and after filtering:

```
zgrep -v "^#" mutect.vcf.gz | wc -l
768

zgrep -v "^#" mutect.filtered.vcf.gz | wc -l
768
```
Huh? I thought we were removing false positives here.  The key to understanding this is to look into the VCF at the FILTER column - if you use `zless` to look at your files, you'll note that in the filtered VCF, there are values filled in - `weak_evidence`, `PASS`, etc.  The key to understanding the FILTER field is that any value besides `.` or `PASS` means that the call is being removed/filtered.  It can be useful to keep those around, though, in case we want to "recover" low VAF variants, or understand why a specific site was removed.

Nothing is marked as filtered in the original VCF, so there are 768 variant calls.  Now let's try a more complicated command to figure out how many variants are in the post-filtering version:
```
zgrep -v "^#" mutect.filtered.vcf.gz | cut -f 7 | egrep "PASS|\." | wc -l
194
```

Much better.

To finish this section up, let's check back in on that Varscan workflow that was running. When the CWL finishes, you can scroll up and see a big block of JSON code describing the outputs.  There should be a line that looks something like this and gives the path to the output file you're looking for - the filtered varscan output.

```
   "varscan_pre_and_post_processing.cwl.filtered_vcf": {
      "format": null,
      "location": "/path/to/output/varscan_pre_and_post_processing.cwl/ead74f9b-811f-4f53-be02-49ae06f6c8a5/call-filter/fp_filter.cwl/97a3088d-efda-479a-a98b-3871d98c12e2/call-hard_filter/execution/varscan_filtered.vcf.gz",
```

Copy that file to your working directory, and also grab the `.tbi` file next to it.  This file contains an index of the VCF that lets tools jump quickly to any point in the file (similar to the index next to our bam files).  

- How many variants are contained in this file, after several filtering steps?


### Merging Variants
To complete this exercise, let's merge these two VCFs into one output file. We're going to accomplish this using the CombineVariants tool bundled with GATK.  Drop into a shell in this container:  `mgibio/gatk-cwl:3.6.0`

Then run the command:

```
/usr/bin/java -Xmx8g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan -o combined.vcf.gz -R /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa --variant:mutect mutect.filtered.vcf.gz --variant:varscan varscan_filtered.vcf.gz
```
Notice that we have to tell the tool what to do with overlaps. If both VCFs contain the same variant call, which line do we keep?  Different tools output different types of auxiliary information, and so we have to choose one. It also adds the "set" attribute, which tells which item came from which VCF (or contains a value like "FilteredInAll" which tells us that it didn't pass anywhere.

- Take a look through this file - how many variants (if any) were called by both callers?

The differences in variants called reflect different underlying statistical models and assumptions.  In this case, it also reflects that we didn't run through the whole Mutect filtering process as implemented in our pipelines. Still, it's important to realize that different callers have different strengths and weaknesses!
 
### Wrapup
If all of this seemed a little complicated, well, that's because it is!  Somatic variant calling, filtering, and prioritization is a difficult problem, given the challenges of purity, ploidy, and heterogeneity that we talked about in the lecture.  Thankfully, we've put together a pretty solid pipeline in the [analysis-workflows](https://github.com/genome/analysis-workflows) repository that handles a lot of the complexity.  You still need to be aware of some of the underlying issues and relevant parameters, as we discussed above, but using a pre-built pipeline can save you a lot of time.




