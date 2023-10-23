## Week 06: Variant calling

This week, we cover how to detect small SNPs/Indels in genomes by using variant calling

- [Lecture Recording]()

- [Slides](bfx_workshop_06_variant_calling.pdf)


## Assignment

### Gather your inputs. 

Start by gathering some data. Pull down a set of input data from (human build38) from this location:
[https://storage.googleapis.com/bfx_workshop_tmp/inputs.tar.gz](https://storage.googleapis.com/bfx_workshop_tmp/inputs.tar.gz)

Unwad the inputs tarball with `tar -xvf` and look at the files:

  - There are two bams, containing reads from chromosome 17 of a human cell line tumor/normal pair.
  - A small reference fasta for chromosome 17 is also included. This is to save time and memory (as opposed to using the whole human reference)
  - You might wonder why we have duplicate indices: `tumor.bai` and `tumor.bam.bai` files here. The short answer is that they contain the same information. Some tools strongly prefer one format, and some the other, so to avoid complications in our pipelines, we just keep both around.

## Detecting Germline Variants

Generate a germline variant VCF from an aligned sequence bam and examine it.

1. Use GATK HaplotypeCaller to produce a VCF file of germline variant calls.  Some hints:

  - Pull and launch docker images to run your commands, just like you did in previous weeks. 
 - Either using a jupyter notebook or running things from the command line directly is fine! 
  - Be aware of what directories you're mounting using the `-v` parameter - all data has to be "underneath" the location you mount!
 - you can find the docker containers for [GATK on Dockerhub](https://hub.docker.com/r/broadinstitute/gatk/).  Use version 4.4.0  (click on the "tags" link at top)
 - you only need three required parameters for `./gatk HaplotypeCaller`
   -  `--input`
   -  `--output`
   -  `--reference`
 - all the files you need are in the inputs directory

2. Examine the ouput VCF by using `less` or your favorite text editor.
    - Note how each field in the VCF has a definition in the header.
    - Find a variant that has at least 100 reads of support for the variant allele
    - Find a site that is an indel (insertion or deletion)
3. Use command line utilities to tell you how many variants were called (don't count the header lines!)


# Detecting Somatic Variants


### Somatic Mutation Calling
We're going to run two different somatic variant calling algorithms today - Mutect and Varscan.  Let's try running the core callers by hand first.

Use docker to fire up the same image we used for germline calling: `broadinstitute/gatk:4.4.0.0`.  This image also contains the Mutect somatic variant caller. Once inside, run Mutect like this, remembering to change the input paths as needed:

```
/gatk/gatk Mutect2 --java-options "-Xmx4g" -O /data/mutect.vcf.gz -R /data/inputs/chr17.fa -I /data/inputs/tumor.bam -tumor "Exome_Tumor" -I /data/inputs/normal.bam -normal "Exome_Normal"
```

Next, let's run Varscan, a different variant caller. For this, we'll use a different docker container: `mgibio/varscan_helper-cwl:1.0.0`  Varscan can be run like this:

```
java -jar /opt/varscan/VarScan.jar somatic <(/opt/samtools/bin/samtools mpileup --no-baq -f /data/inputs/chr17.fa /data/inputs/normal.bam /data/inputs/tumor.bam) varscan.vcf --strand-filter 0 --min-coverage 8 --min-var-freq 0.1 --p-value 0.99 --mpileup 1 --output-vcf
```
Note the piping on the command line with `<()` that takes the input of `samtools mpileup` and feeds it to varscan.

- Which caller was faster?  (and remember, we're only calling variants on a very small portion of the genome!!)

Both callers typically require post-processing to give only high-quality variant calls, but let's take a quick look at the raw outputs. Mutect outputs a gzipped files, so you'll have to either run `gunzip -c mutect.vcf.gz | less`, or for a shortcut, just use `zless mutect.vcf.gz`

- Which caller produced more variants in the VCF? (don't include the header!)
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

### Filtering a VCF file

Filtering is generally needed to get these "raw" calls down to something with a better specificity. There are multiple layers of filtering in our somatic pipeline, but here we'll just run the first-pass, built-in Mutect filter, using the same GATK container as above:

```
gatk FilterMutectCalls -R /data/inputs/chr17.fa -V /data/mutect.vcf.gz -O /data/mutect.filtered.vcf.gz
```

If you run FilterMutectCalls without inputs, you can see the many parameters that can be tweaked and fine-tuned as needed.

After filtering, you'll find some detailed stats in `mutect.filtered.vcf.gz.filteringStats.tsv` afterwards.  Now, compare the number of mutations in these VCFs before and after filtering:

```
zgrep -v "^#" mutect.vcf.gz | wc -l
768

zgrep -v "^#" mutect.filtered.vcf.gz | wc -l
768
```
Huh? I thought we were removing false positives here.  To resolve this mystery, look in the VCF at the FILTER column. If you use `zless` to look at your files, you'll note that in the filtered VCF, there are values filled in - `weak_evidence`, `PASS`, etc.  The key to understanding the FILTER field is that any value besides `.` or `PASS` means that the call is being removed/filtered.  It can be useful to keep those around, though, in case we want to "recover" low VAF variants, or understand why a specific site was removed.

Nothing is marked as filtered in the original VCF, so there are 768 variant calls.  Now let's try a more complicated command to figure out how many variants are in the post-filtering version:
```
zgrep -v "^#" mutect.filtered.vcf.gz | cut -f 7 | egrep "PASS|\." | wc -l
194
```

Much better.


### Combining VCF files

Varscan is a bit more complicated because we have separate files for SNVs and Indels.  Let's merge them to make things easier using these two commands to zip up and index our varscan variant calls:

```
docker run -v $(pwd -P):/data biocontainers/tabix:v1.9-11-deb_cv1 /bin/sh -c "bgzip varscan.vcf.snp;tabix -p vcf varscan.vcf.snp.gz"
docker run -v $(pwd -P):/data biocontainers/tabix:v1.9-11-deb_cv1 /bin/sh -c "bgzip varscan.vcf.indel;tabix -p vcf varscan.vcf.indel.gz"

```

this BCFtools command to merge them: 

```
docker run -v $(pwd -P):/data mgibio/bcftools-cwl:1.12 /opt/bcftools/bin/bcftools concat --allow-overlaps --remove-duplicates --output-type z -o /data/varscan_merged.vcf.gz /data/varscan.vcf.snp.gz /data/varscan.vcf.indel.gz
```

and finally, let's index our resulting merged varscan VCF:

```
docker run -v $(pwd -P):/data biocontainers/tabix:v1.9-11-deb_cv1 tabix /data/varscan_merged.vcf.gz
```

Note that we didn't do any filtering on these Varscan ouputs. That isn't because filtering isn't needed (the number of variants is ludicrously high in the raw outputs!)  We're just omitting it in the interest of keeping this tutorial brief!


### Merging Variants from different callers
To complete this exercise, let's merge these VCFs from two different callers into one combined output file. We're going to accomplish this using the CombineVariants tool bundled with older versions of GATK.  Drop into this container:  `mgibio/gatk-cwl:3.6.0`

Then run the command:

```
/usr/bin/java -Xmx2g -jar /opt/GenomeAnalysisTK.jar -T CombineVariants -genotypeMergeOptions PRIORITIZE --rod_priority_list mutect,varscan -o /data/combined.vcf.gz -R /data/inputs/chr17.fa --variant:mutect /data/mutect.filtered.vcf.gz --variant:varscan /data/varscan_merged.vcf.gz
```

Notice that we have to tell the tool what to do with overlaps. If both VCFs contain the same variant call, which line do we keep?  Different tools output different types of auxiliary information, and so we have to choose one. It also adds the "set" attribute, which tells which item came from which VCF (or contains a value like "FilteredInAll" which tells us that it didn't pass anywhere.

- Take a look through this file - how many variants (if any) were called by both callers?

The differences in variants called reflect different underlying statistical models and assumptions. In this case, it also reflects that we didn't run through the whole Mutect filtering process as implemented in our pipelines. Still, it's important to realize that different callers have different strengths and weaknesses!



 
## Wrap up
If all of this seems a little complicated, well, that's because it is!  Variant calling, filtering, and prioritization is a difficult problem, and doubly so for somatic variants - the challenges of purity, ploidy, and heterogeneity that we talked about in the lecture are difficult to overcome.  

Thankfully, solid pipelines exist for running variant calling. Pipelines help you in two ways: 1) they encapsulate many of these steps, and just allow you specify needed inputs and say "Go".  2) They benefit from years of fine-tuning from folks who have thought deeply about the problem. Many such pipelines exist, for one example, you could look at the CWL pipelines in the [analysis-workflows](https://github.com/genome/analysis-workflows) repository or comparable WDL workflows in the [analysis-wdls](https://github.com/wustl-oncology/analysis-wdls/issues) repository. You still need to be aware of some of the underlying issues and relevant parameters, as we discussed above, but using a pre-built pipeline can save you a lot of time!


## Visualization (optional)

If you'd like to see the results of your hard work, you can take a look at some of these variants in IGV. 

1.  Open IGV, make sure you have genome build hg38 set, then load your `normal.bam`, `tumor.bam`, and the `combined.vcf.gz` using "Open File".  
2. Go to the position of the first variant in your VCF.  Do you have high confidence that this variant is real? (look at the mapping quality of the reads)
3. Now jump to position `chr17:7674764-7674830` in the TP53 gene.  Does this look like a somatic (tumor-specific) variant to you?
4. Click on the name of your VCF file in the lefthand panel, then hit "CTRL-F".  Note how it jumps to the next variant in the file.  Does this one look more convincingly somatic?  
5. Click on the coverage track to view the readcounts and VAF of this variant. What could the VAF indicate about this tumor's purity or ploidy?
6. Finally, zoom out until you can see several exons of TP53, but still see the coverage track. Notice the "bumps" in coverage over each exon that indicate this was Exome sequencing, not whole genome. 