# Germline Variant Calling and Filtering

This week we use the GATK HaplotypeCaller to call germline variants from cell line bams, combine them with public data from the 1000 Genomes Project, and filter the results. For a more in-depth look, see this excellent [GATK tutorial](https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/tutorials/7869-howto-Discover-variants-with-GATK-A-GATK-Workshop-Tutorial), provided by the Broad Institute.

By the end of this exercise you should be able to:

- Call germline variants on a single sample with GATK HaplotypeCaller
- Understand what joint genotyping is, why it matters, and run it on our sample plus public 1000 Genomes Project exomes
- Manipulate and filter VCFs to remove artifacts and identify variants of interest
- Manage long-running jobs at the command line, and know when to hand them off to Slurm instead

> [!NOTE]
> As always, we strongly suggest typing the commands out yourself rather than copying and pasting. The exceptions are the very long GATK commands, where a typo can cost you a lot of time. Copy those, but read them carefully first so you know what each piece does.

We are using [GATK v4.7.0.0](https://github.com/broadinstitute/gatk/releases/tag/4.7.0.0) for this tutorial.

## Contents

1. [Setup](#1-setup)
2. [Working with long-running jobs](#2-working-with-long-running-jobs)
3. [Single-sample calling with GATK HaplotypeCaller](#3-single-sample-calling-with-gatk-haplotypecaller)
4. [Joint genotyping](#4-joint-genotyping)
5. [Filtering your VCF](#5-filtering-your-vcf)
6. [Cleaning up](#6-cleaning-up)
7. [Cheat sheet](#7-cheat-sheet)

---

## 1. Setup

### The container

GATK (along with samtools, bcftools, tabix, and bgzip) is packaged in the Broad Institute's official docker image:

```
broadinstitute/gatk:4.7.0.0
```

> [!NOTE]
> This image is big (about 2.5 GB to download, and more once unpacked). The first time you launch it, it may take several minutes to download and unpack. Be patient, and don't cancel it!

### Set up your workspace

Log in to Compute2, and make a directory for this week:

```bash
ssh <washukey>@c2-login-002.ris.wustl.edu
mkdir -p ~/workshop/week11/germline
```
(as always, replace \<washukey\> with your username)

Just like the RNA-seq weeks, we're going to mount this week's directory inside the container at `/workspace`, so that the paths in the commands below are shorter and the same for everybody.

### Launch an interactive job in the container

```bash
srun -A compute2-workshop -p workshop -c 2 --mem=16G --pty \
  --container-image=broadinstitute/gatk:4.7.0.0 \
  --container-mounts="/storage1/fs1/c.a.miller/Active/bfx-workshop-scratch/$USER/week11:/workspace" \
  /bin/bash
```

This should all look familiar from the last few weeks. We're asking for a little more than usual (16 GB of RAM), because GATK is memory-hungry, and we're going to run a couple of HaplotypeCaller jobs at the same time.

Once you have a prompt, check that GATK is there and that your workspace is mounted:

```bash
gatk --version
ls /workspace
```

You should see a version number ending in `4.7.0.0` (along with some other library versions), and the `germline` directory.

> [!IMPORTANT]
> Unless otherwise noted, **all of the commands for the rest of this assignment are run inside this container**. If you get disconnected or `exit`, re-run the `srun` command above. Everything you've created will still be there, because it was written to your storage directory, not into the container.

### Download the data

Make a directory to hold log files:

```bash
cd /workspace/germline
mkdir -p logs
```

Grab the bam files from the HCC1395 cell line, which have already been aligned, duplicate-marked, and base quality recalibrated:

```bash
wget https://storage.googleapis.com/bfx_workshop_tmp/align.tar
tar -xvf align.tar
ls -lh align/
rm align.tar
```

This one is about 3 GB, so it'll take a few minutes to download. It contains exome and WGS bams for the normal sample (plus an exome bam for the tumor that we won't use this week), each with a `.bai` index.

We'll also need the reference genome:

```bash
mkdir -p /workspace/germline/reference
cd /workspace/germline/reference
wget https://genomedata.org/pmbio-workshop/references/genome/chr6_and_chr17/ref_genome.tar
tar -xvf ref_genome.tar
gunzip ref_genome.fa.gz
rm ref_genome.tar
ls
```

Notice that, along with the FASTA, the tarball came with a `.fai` index and a `.dict` sequence dictionary. GATK requires both of these, and will refuse to run without them.

Rather than calling variants on a whole genome, which takes hours, we'll restrict ourselves to a handful of regions on chromosome 17 (about 5 Mb in total) that contain the loci we'll look at in IGV. We'll provide them in a BED file:

```bash
cd /workspace/germline
wget https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/regions.bed
cat regions.bed
```

---

## 2. Working with long-running jobs

Variant calling on real data takes hours to days. We don't want to sit and stare at a terminal while that happens. Your terminal is capable of running things in the **background** to address that.

Start by running the `sleep` command. It runs in the foreground, meaning that your prompt is blocked until it finishes:

```bash
sleep 300
```

Press **Ctrl-C** to kill it.

Now run it again:

```bash
sleep 300
```

But this time press **Ctrl-Z**. That *suspends* (pauses) the job.

Now list your jobs, then use the `bg` command to resume that job in the background:

```bash
jobs
# [1]+  Stopped         sleep 300

bg 1
# [1]+ sleep 300 &

jobs
# [1]+  Running         sleep 300 &
```

Now it's running in the background. If we want to bring it back under our control, we can use `fg` to move it to the foreground:

```bash
fg 1
```

Now kill the job with **Ctrl-C**.

If you know in advance a job will be slow, you can start it in the background directly with `&`:

```bash
sleep 300 &
```

Background jobs still print to your terminal, though, which can get messy, so we generally make sure to redirect outputs to log files:

```bash
sleep 60 >logs/sleep.log 2>logs/sleep.err &
jobs
```

Reminder: `>` redirects STDOUT to a file, and `2>` sends STDERR to a file.

### Background jobs vs. Slurm jobs

There's an important catch on the cluster. A background job is a child of the shell you started it from, and that shell lives inside your `srun` session. **If you `exit`, lose your connection, or your `srun` job hits its time limit, every background job you started goes down with it.**

So `&` and `bg` are great for keeping your prompt free while a job runs for a few minutes. For anything that needs to run for hours, or that you want to keep running after you log off, submit it to Slurm with `sbatch` instead, as we did in [week 2](../week_02/compute2-slurm-intro.md#7-batch-jobs-sbatch). We'll do both below.

---

## 3. Single-sample calling with GATK HaplotypeCaller

First, let's run [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/4404604697243-HaplotypeCaller) to call germline SNPs and indels. Whenever HaplotypeCaller finds signs of variation, it performs a local *de novo* reassembly of reads. This improves the accuracy of variant calling, especially in challenging regions. We'll also include an option to generate bam output from HaplotypeCaller so that local reassembly/realignment around called variants can be visualized.

We'll run GATK HaplotypeCaller with the following options:

* `--java-options '-Xmx5g'` tells GATK (which is written in Java) to use up to 5 GB of memory
* `HaplotypeCaller` specifies the GATK command to run
* `-R` specifies the path to the reference genome
* `-I` specifies the path to the input bam file for which to call variants
* `-O` specifies the path to the output VCF file to write variants to
* `--bam-output` specifies the path to an optional bam file which will store local realignments around variants that HaplotypeCaller used to make calls
* `-L` is used to limit calling to specific regions. In our case, this is the `regions.bed` file that we pulled down earlier. To run on the whole sample, you could just omit it.

All right, let's run this tool on the exome:

```bash
cd /workspace/germline

# Call variants for exome data
gatk --java-options '-Xmx5g' HaplotypeCaller \
    -R /workspace/germline/reference/ref_genome.fa \
    -I align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    -O /workspace/germline/Exome_Norm_HC_calls.vcf \
    --bam-output /workspace/germline/Exome_Norm_HC_out.bam \
    -L /workspace/germline/regions.bed 2>logs/hc_exome.err >logs/hc_exome.out
```

Looks like that's going to take a while. Let's send it to the background. Hit **Ctrl-Z** to suspend the job, then run `jobs` to find its ID. Finally, use `bg 1` (replacing 1 if needed) to send it to the background.

We're also going to want to create a GVCF (more on that later on). Let's go ahead and get that chugging along in the background too, this time starting it there directly with `&`:

```bash
gatk --java-options '-Xmx5g' HaplotypeCaller \
    -ERC GVCF \
    -R /workspace/germline/reference/ref_genome.fa \
    -I align/Exome_Norm_sorted_mrkdup_bqsr.bam \
    -O /workspace/germline/Exome_Norm_HC_calls.g.vcf \
    --bam-output /workspace/germline/Exome_Norm_HC_GVCF_out.bam \
    -L /workspace/germline/regions.bed 2>logs/hc_gvcf.err >logs/hc_gvcf.out &
```

You can check on their progress at any time by looking at the end of the log files (e.g. `tail logs/hc_exome.err`). GATK prints a line every so often telling you which genomic position it has reached.

### Submitting the WGS calling as a Slurm job

Calling variants on the WGS data is definitely going to take a while. We're not using `-L` this time, so it will look at *all* of chr6 and chr17, and there's a lot more data. This is exactly the kind of job that we don't want to tie to our interactive session, so let's hand it off to Slurm.

`sbatch` isn't available inside the container, so **open a second terminal window** and log in to Compute2 again. (Leave your first window, with the container session, open.) Then, on the login node:

```bash
cd ~/workshop/week11/germline

sbatch -A compute2-workshop -p general-cpu -t 4:00:00 -c 2 --mem=8G -J wgs_hc \
  --container-image=broadinstitute/gatk:4.7.0.0 \
  --container-mounts="/rdcw/fs2/c.a.miller/Active/bfx-workshop-scratch/$USER/week11:/workspace" \
  --container-workdir=/workspace/germline \
  --output=logs/hc_wgs.out --error=logs/hc_wgs.err \
  --wrap="gatk --java-options '-Xmx5g' HaplotypeCaller \
    -R /workspace/germline/reference/ref_genome.fa \
    -I align/WGS_Norm_merged_sorted_mrkdup_bqsr.bam \
    -O /workspace/germline/WGS_Norm_HC_calls.vcf \
    --bam-output /workspace/germline/WGS_Norm_HC_out.bam"
```

That's a mouthful, but every piece of it is something you've seen before:

- `-A`, `-p`, `-c`, `--mem`, `-J`, and `-t` (a 4-hour time limit) are the same options we used in [week 2](../week_02/compute2-slurm-intro.md#7-batch-jobs-sbatch)
- `--container-image`, `--container-mounts`, and `--container-workdir` are the same as in our `srun` command, so the job sees the same files at the same paths that we do
- `--output` and `--error` say where to put STDOUT and STDERR. Unlike everything inside `--wrap`, these are written by Slurm *outside* the container, so they're relative to the directory we submitted from (`~/workshop/week11/germline`)
- `--wrap` holds the command to run, which is the same HaplotypeCaller command as above, minus `-L`

Check that it's running with:

```bash
squeue -u $USER
```

You should see your `wgs_hc` job, alongside the `srun` session in your other window. Now that it's in Slurm's hands, it'll keep going even if you log out. You can close this second window, and head back to your container session.

------------------------------------------------

### Exercise - Part 1

Use `jobs` in your container session to check if the exome jobs are done. (Don't worry about the WGS job yet).  Once they're done, verify that they ran to completion by looking at the logs:

```
tail logs/hc_exome.err
```
You should see something like: 

```
. . . HaplotypeCaller done. Elapsed time: 0.83 minutes.
```

Let's start by looking through the exome VCF file, `Exome_Norm_HC_calls.vcf`:

1. How many variants were called?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    Each variant is one non-header line, so count the lines that don't start with `#`:

    ```bash
    grep -v "^#" Exome_Norm_HC_calls.vcf | wc -l
    ```

    `-v` inverts the match (lines that *don't* start with `#`)

    </td></tr></table>
    </details>

2. For the first variant in the VCF file, how many reads are listed as supporting the reference allele and how many support the variant allele?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    Print the first non-header line:

    ```bash
    grep -v "^#" Exome_Norm_HC_calls.vcf | head -n 1
    ```

    The FORMAT column (second to last) tells you how to read the sample column (last). It looks something like `GT:AD:DP:GQ:PL`. The `AD` (allelic depth) field holds the read counts in the order `REF,ALT`, so a sample value of `0/1:4,3:...` means 4 reads support the reference allele and 3 support the variant.

    </td></tr></table>
    </details>

3. Open IGV on your laptop, and view the bams the same way we did in [week 4](../week_04/alignment_assignment.md#5-viewing-the-data-with-igv):
    * Mount your storage on your laptop (Mac: `open smb://storage1.ris.wustl.edu/c.a.miller/Active/bfx-workshop-scratch/<washukey>`; Windows: `\\storage1.ris.wustl.edu\c.a.miller\Active\bfx-workshop-scratch\<washukey>`)
    * In IGV, choose the **Human (hg38)** genome. (Like our week 4 reference, this one has full-length chr6 and chr17 with GRCh38 coordinates, so it lines up with the standard genome.)
    * Use **File > Load from File...** to load these two bams from the `week11/germline/` folder:
        * `align/Exome_Norm_sorted_mrkdup_bqsr.bam`
        * `Exome_Norm_HC_out.bam`

4. Look up the location of the first variant in the VCF file and then navigate there in IGV. Why doesn't the number of reads match the allelic depth from the VCF?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    - It ignores reads with low mapping quality (MAPQ < 20), bases with low quality, or that fail other QC metrics
    - Duplicate-marked reads are dropped (IGV hides these by default too).

    The HaplotypeCaller bam track is a closer match to `AD`, because it shows the reads the caller didn't filter out

    </td></tr></table>
    </details>

5. Navigate to `chr17:76287089`.

6. Some variant callers might report `chr17:76287089:C>T` as an event. Would that be correct?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    No. In the bwa-mem track the "C>T" mismatches pile up at the ends of reads, in a CA dinucleotide repeat, and right next to soft clips, which is one classic sign of an alignment artifact.

    </td></tr></table>
    </details>

7. What information do the soft-clipped bases give you?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    Soft-clipped bases are part of the read but didn't match the reference where the aligner put it. When many reads clip at the same position, it tells you the sequence beyond that point is different from the reference. Here, the clipped bases match the reference sequence *after* the deleted stretch, so the reads are really spanning a deletion that the aligner didn't model. (Clusters of soft clips are also how structural variant callers find breakpoints.)

    </td></tr></table>
    </details>

8. What happens if you zoom out to a 2000 bp window or so? Why does the bwa-mem alignment have more reads?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    The HaplotypeCaller track only has reads in a small window around the event, while the bwa-mem track has reads everywhere there's coverage. HaplotypeCaller only reassembles "active regions" where it sees signs of variation, and `--bam-output` writes out only the reads it used there. It also leaves out reads it filtered (low MAPQ, duplicates, etc.) and can downsample very deep sites, so even inside the active region it can show fewer reads than bwa-mem.

    </td></tr></table>
    </details>

9. View the alignments in a "squished" format to see more reads at a time. What's your interpretation of this event?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    (Right-click the track and choose `Squished`.) In the HaplotypeCaller track, the reassembled reads line up cleanly with a 28 bp gap, and the scattered mismatches and soft clips from the bwa-mem track disappear. It's a single real deletion in a CA repeat, not a SNV. Repeats like this are hard for read-by-read aligners, and local reassembly is exactly what HaplotypeCaller is designed to fix. Compare the fraction of reads carrying the deletion with the reads that don't, to judge whether it looks heterozygous or homozygous.

    </td></tr></table>
    </details>

10. Find the `chr17:76287089:C>T` site in your VCF file. Does the genotype field (GT) match your expectation?

    <details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
    <table><tr><td markdown="1">

    ```bash
    grep 76287089 Exome_Norm_HC_calls.vcf | grep "chr17" | head
    ```

    There's no C>T SNV at 76287089. Instead, HaplotypeCaller reports the deletion at `chr17:76287061`, with a long `REF` allele and a short `ALT` allele. Its `GT` should agree with what you saw in IGV: `0/1` if roughly half the reads carry the deletion, `1/1` if nearly all do. Check the `AD` field to see the read counts behind that call.

    </td></tr></table>
    </details>

---

## 4. Joint genotyping

So far we've looked at variants from one sample, but often we want to look across a population. To do so, we typically use joint genotyping, where all samples are genotyped *together* at every site where *any* sample shows evidence of variation. That buys us:

- **A complete matrix.** At every variant site, each sample gets a genotype: `0/0` (confidently reference), `0/1`, `1/1`, or `./.` (no data). With single-sample VCFs, a site missing from one file could mean "reference" *or* "not covered," and you can't tell which.
- **Better sensitivity.** A variant with weak support in one sample can be called confidently when other samples show it too.
- **Cohort-level statistics** used by filtering methods such as VQSR.

Why go through GVCFs instead of feeding all the bams to HaplotypeCaller at once? Because of the **N+1 problem**: with multi-bam calling, adding one new sample means redoing *everything*. With GVCFs, each sample's expensive step is done once, and only the cheap joint step is repeated.

### Examining a GVCF

First, make sure your GVCF job has finished (check with `jobs`, or look at the end of `logs/hc_gvcf.err`). Then let's compare the GVCF to the standard VCF you made from the same exome. Start by skipping past the `##` header lines and looking at the first few records of each:

```bash
cd /workspace/germline
grep -v "^##" Exome_Norm_HC_calls.vcf | head -10
grep -v "^##" Exome_Norm_HC_calls.g.vcf | head -30
```

- Most GVCF records have `<NON_REF>` as the ALT allele and an `END=` tag in the INFO field. What do you think those records represent?
- Look at the `GQ` and `MIN_DP` values in those records. Why would GATK split the reference into many separate blocks instead of one big one? (Hint: look up the `--gvcf-gq-bands` option.)

Now count the records in each file:

```bash
grep -v "^#" Exome_Norm_HC_calls.vcf | wc -l
grep -v "^#" Exome_Norm_HC_calls.g.vcf | wc -l
```

- Why does the GVCF have so many more records?

Finally, pull out just the GVCF records that contain a real ALT allele (the ones that aren't only `<NON_REF>`):

```bash
grep -v "^#" Exome_Norm_HC_calls.g.vcf \
    | awk '$5 != "<NON_REF>"' \
    | head -5
```

- Compare one of these to the same site in `Exome_Norm_HC_calls.vcf`. What extra allele is listed in the GVCF's ALT column?
- Suppose another sample has a variant at a position that falls inside one of HCC1395's `<NON_REF>` blocks. What can the GVCF tell us about HCC1395 at that site that the standard VCF can't?

### Running joint genotyping on public data

Let's pull down some data, get joint genotyping running, and then we'll talk about where it came from. In advance, we downloaded exome alignments for 5 females of GBR descent, subset them to the same regions as our cell line data, and then performed germline variant calling using the GATK HaplotypeCaller GVCF workflow commands, as above. The details for this can be found on the [Developer's Notes](https://pmbio.org/module-10-appendix/0010/02/01/Developer_Notes/) page.

Go ahead and grab these GVCFs and their index files (they're about 120 MB each):

```bash
cd /workspace/germline
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00099_HC_calls.g.vcf
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00099_HC_calls.g.vcf.idx
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00102_HC_calls.g.vcf
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00102_HC_calls.g.vcf.idx
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00104_HC_calls.g.vcf
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00104_HC_calls.g.vcf.idx
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00150_HC_calls.g.vcf
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00150_HC_calls.g.vcf.idx
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00158_HC_calls.g.vcf
wget http://genomedata.org/pmbio-workshop/references/1KGP/gvcfs/chr17/HG00158_HC_calls.g.vcf.idx
```

To create a joint genotype GVCF, we'll combine our HCC1395 exome normal with these five using the `CombineGVCFs` command.

GATK CombineGVCFs has the following new options:

* `-V` [multiple] specifies the path to each of the input g.vcf files
* `-O` specifies the path to the output combined g.vcf file

```bash
gatk --java-options '-Xmx8g' CombineGVCFs \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.g.vcf \
    -V /workspace/germline/HG00099_HC_calls.g.vcf \
    -V /workspace/germline/HG00102_HC_calls.g.vcf \
    -V /workspace/germline/HG00104_HC_calls.g.vcf \
    -V /workspace/germline/HG00150_HC_calls.g.vcf \
    -V /workspace/germline/HG00158_HC_calls.g.vcf \
    -O /workspace/germline/Exome_Norm_1KGP_HC_calls_combined.g.vcf \
    -L /workspace/germline/regions.bed
```

This'll take about two minutes or so.

Next, perform joint genotyping with the `GenotypeGVCFs` command. This command is run with the following new options:

* `-V` specifies the path to the combined g.vcf file
* `-O` specifies the path to the output joint-genotyped VCF file

```bash
gatk --java-options '-Xmx12g' GenotypeGVCFs \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_1KGP_HC_calls_combined.g.vcf \
    -O /workspace/germline/Exome_GGVCFs_jointcalls.vcf \
    -L /workspace/germline/regions.bed
```

> [!NOTE]
> We asked Slurm for 16 GB of RAM, so these 8 and 12 GB GATK commands fit comfortably on their own. If you still have a 5 GB HaplotypeCaller job running in the background, wait for it to finish first (check with `jobs`). If a command dies with a message like `Killed` or `oom-kill`, it ran out of memory: exit the container, re-launch it with `--mem=24G`, and try again.

---

### Sidebar: Finding public data

Including additional samples and doing joint genotyping can improve the overall accuracy of variant calling and filtering. In this example, we only had one "normal" sample, our HCC1395BL cell line. Where did we get the others?

The 1000 Genomes Project is an international effort to create a set of publicly available sequence data for a large cohort of ethnically diverse individuals. Recall that the cell line being used in this tutorial ([HCC1395](https://www.atcc.org/Products/All/CRL-2324.aspx)) was derived from a 43-year-old Caucasian female (by a research group in Dallas, Texas). Let's see if we can find some suitable "matching" samples in the 1000 Genomes data.

As you do this, keep in mind the FAIR principles that we talked about earlier in the course (Findable, Accessible, Interoperable, Reusable).

This part happens in a web browser on your laptop. When it asks you to download a list, save it into your mounted `week11/germline/` folder, so that you can look at it from your container session too.

#### Find available samples in the 1000 Genomes cohort

- Hop onto the [1000 Genomes Samples](https://www.internationalgenome.org/data-portal/sample) page and scroll through the list of samples available.
- In a real study, you'd spend a bunch of time thinking about which samples were the right ones to use, to make sure that you're finding disease-related variants and not just population markers. In this case, we'll just assume that a few samples from the *British in England and Scotland* population will be a reasonable match.
- Filter by population: "British in England and Scotland"
- Filter by data collection: "1000 Genomes on GRCh38"
- Filter by technology: "Low Coverage WGS" (exomes aren't available for these samples)
- Download the list (save it as `igsr_samples_GBR.tsv`) and open it. (**Findable**!)

Does this have enough data to get at the variant calls?

#### Get a list of data files for GBR samples

- Click on "Data Collections" at the top
- Choose "1000 Genomes 30x on GRCh38"
- In the Populations box, choose "British in England and Scotland"
- Scroll down to the bottom and filter the file list to "Alignments"
- Click "Download the list" (save it as `igsr_GBR_GRCh38_alignment.tsv`)

Now, back in your container session, find the URL for HG00099's CRAM:

```bash
cd /workspace/germline
grep HG00099 igsr_GBR_GRCh38_alignment.tsv | cut -f 1
```

With this list, we could download the already aligned data for a bunch of 1000 Genomes individuals in CRAM format. (**Accessible**!)

How do we know if this data is suitable to use with our own data, though?

#### Find technical details about the data

This takes a little digging. On the 1000 Genomes site, go to the top, then click through to Help -> FAQ -> "What methods were used for generating alignments?". That directs you back to the [data collection page](https://www.internationalgenome.org/data-portal/data-collection), and for this "1000 Genomes on GRCh38" study, to the [latest 2019 publication here](https://wellcomeopenresearch.org/articles/4-50).

That was quite a journey, but if you scroll down, you'll notice that they do have excellent documentation of the tools, input files, and parameters used. (**Reusable**!)

With this info, we can determine that these CRAM files were created in a way that is generally compatible with the analysis done so far in this tutorial. (**Interoperable**!)

We still need to consider that there could be batch effects related to potentially different sample preparation, library construction, exome capture reagent and protocol, sequencing pipeline, etc. These are complex topics that need to be built into the statistical models, and are a bit beyond the scope of what we're going to talk about in this course!

---

### Exercise - Part 2: Joint genotyping results

- Look at this joint genotyping VCF (`Exome_GGVCFs_jointcalls.vcf`) and compare it to the "standard" VCF you created above (`Exome_Norm_HC_calls.vcf`). What are the big differences?
- Load the VCF in IGV (from the `week11/germline/` folder) and take a look at what pops up.

### Tidy up

Let's clean up this directory a bit to make our lives easier:

```bash
cd /workspace/germline
mkdir -p 1kgenomes
mv HG0* 1kgenomes
```

Now check in on the WGS job. `squeue` doesn't work inside the container, so use the second terminal window (or open a new one) and run `squeue -u $USER` on the login node. If `wgs_hc` is no longer listed, it's finished. Then, back in the container, look at the end of the log:

```bash
tail logs/hc_wgs.err
```

Did it finish successfully? How long did it take? How would this scale if we were running it on the whole genome?

> [!TIP]
> On the login node, `sacct -j <jobid> --format=JobID,JobName,Elapsed,MaxRSS,State` will tell you how long a finished job ran, and how much memory it used at its peak. (`squeue` and `sbatch` print the job ID.)

**Once the job has finished** (not before, or you'll pull the files out from under it), tidy up the WGS data as well:

```bash
mkdir -p wgs
mv WGS* logs/hc_wgs* wgs/
```

---

## 5. Filtering your VCF

The raw output of GATK HaplotypeCaller will include many variants with varying degrees of quality. Doing additional filtering is generally a good idea. The details of VQSR and hard-filtering strategies are covered in some detail on the [PMBio variant filtering lesson](https://pmbio.org/module-04-germline/0004/02/02/Germline_SnvIndel_FilteringAnnotationReview/) and the GATK site: [How to Filter variants either with VQSR or by hard-filtering](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering).

### Filtering with GATK

In the interest of time, we're going to set up a bash script to breeze through these steps. On your own data, you probably would be following a pre-designed workflow anyway!

The container doesn't have `nano` (it does have `vi`, if you know it), so we'll write the script with a heredoc, just like we wrote our environment file in [week 8](../week_08/rnaseq_alignment.md#4-set-up-the-course-environment-variables): everything between `<<'EOF'` and the final `EOF` is written into `filter.sh`. Take a look through the steps of this script, then copy and paste the whole block into your container session:

```bash
cd /workspace/germline
cat > filter.sh <<'EOF'
# This line provides several safeguards that instruct the script to fail
# if things go wrong (if any command fails, if any variables are undefined, etc.)
set -euo pipefail

# Split the data into SNPs and INDELs to be filtered separately:
cd /workspace/germline/
gatk --java-options '-Xmx12g' SelectVariants \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.vcf \
    -select-type SNP \
    -O /workspace/germline/Exome_Norm_HC_calls.snps.vcf
gatk --java-options '-Xmx12g' SelectVariants \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.vcf \
    -select-type INDEL \
    -O /workspace/germline/Exome_Norm_HC_calls.indels.vcf

# Use hard filtering (specific cutoffs) to label variants
gatk --java-options '-Xmx12g' VariantFiltration \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.snps.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
    --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
    --filter-expression "SOR > 3.0" --filter-name "SOR_gt_3" \
    -O /workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf
gatk --java-options '-Xmx12g' VariantFiltration \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.indels.vcf \
    --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
    --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \
    --filter-expression "SOR > 10.0" --filter-name "SOR_gt_10" \
    -O /workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf

# Merge the filtered SNP and INDEL VCFs back together
gatk --java-options '-Xmx12g' MergeVcfs \
    -I /workspace/germline/Exome_Norm_HC_calls.snps.filtered.vcf \
    -I /workspace/germline/Exome_Norm_HC_calls.indels.filtered.vcf \
    -O /workspace/germline/Exome_Norm_HC_calls.filtered.vcf

# Extract PASSing variants only into a separate file
gatk --java-options '-Xmx12g' SelectVariants \
    -R /workspace/germline/reference/ref_genome.fa \
    -V /workspace/germline/Exome_Norm_HC_calls.filtered.vcf \
    -O /workspace/germline/Exome_Norm_HC_calls.filtered.PASS.vcf \
    --exclude-filtered
EOF
```

Check that the file looks right with `less filter.sh`, then run it:

```bash
bash filter.sh 2>logs/filter.err
```

NOTE: There are a number of MIXED type variants (multi-allelic with both SNP and INDEL alleles) that are currently dropped by the above workflow (selecting for only SNPs and INDELs). In the future we could consider converting these with `gatk LeftAlignAndTrimVariants --split-multi-allelics`.

### Exercise - Part 3

- In IGV, load your filtered single-sample VCF file (`week11/germline/Exome_Norm_HC_calls.filtered.PASS.vcf`). If you still have the unfiltered VCF loaded, remove that track (right-click it and choose **Remove Track**) so you don't mix them up.
- Navigate to the TP53 gene, in which four variants have been called. Inspect each of them, and decide whether the site looks like a solid call to you, and why or why not.
- Navigate to `chr17:7,675,931-7,676,472`. Do you think that GATK got the call right here?
- What can we say about the phase of these three variants? `chr17:5,381,336-5,381,470`
- How about these two? `chr17:82,870,385-82,870,452`
- Or these two? `chr17:7,726,842-7,727,112`

### Filtering with bcftools

The script above uses some pretty heavy-duty GATK commands to filter VCFs. In theory, if we wanted to do this in a more lightweight way, we could write a one-liner with grep or use a complicated regular expression to get the fields we care about, but we *strongly encourage* you not to do that. It's too easy to screw it up, and there are better ways!

One additional method to be aware of is `bcftools`, which is a powerful tool for manipulating VCF files.

Conveniently, `bcftools` comes packaged in the same container as GATK, so you can keep working in your current session. (Check with `bcftools --version`.)

As a quick example, let's just subset one of our VCFs to only keep calls with a mapping quality of at least 55:

```bash
cd /workspace/germline
bcftools view --include 'MQ>55' \
    Exome_Norm_HC_calls.filtered.PASS.vcf >high_mq.vcf
```

Now what if we wanted to combine that with accessing an attribute inside the INFO field (depth of at least 100)?

```bash
bcftools view --include 'INFO/DP>=100 & MQ>55' \
    Exome_Norm_HC_calls.filtered.PASS.vcf >high_mq_high_dp.vcf
```

How many variants did those produce compared to the original? Why do we need the `grep` in this command?

```bash
for file in Exome_Norm_HC_calls.filtered.PASS.vcf high_mq.vcf high_mq_high_dp.vcf; do
    grep -v "^#" $file | wc -l
done
```

How about our multi-sample VCF? One handy way to explore it is to convert the genotypes into a matrix using bcftools:

```bash
bcftools query -Hf '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
    Exome_GGVCFs_jointcalls.vcf > genotype_matrix.tsv
```

Take a look at the file with `less -S genotype_matrix.tsv`. Using it as input, can you find out how many sites have at least one missing genotype?

<details markdown="1"><summary><b>Answer</b> (click to expand)</summary>
<table><tr><td markdown="1">

Missing genotypes are written as `./.`, so we can count the lines that contain one:

```bash
grep -c "\./\." genotype_matrix.tsv
```

We need the backslashes because `.` is a special character in grep that matches *any* character. `\.` matches only a literal period.

</td></tr></table>
</details>

These bcftools commands are nice because they adhere closely to the Unix philosophy: small tools that output to the command line and can be piped to additional downstream steps. That can be convenient, especially for ad hoc work and exploring data.

---

## 6. Cleaning up

When you're done, `exit` your container session so that you're not holding on to a compute allocation. Also check `squeue -u $USER` on the login node, to make sure that you don't have any jobs still running that you've forgotten about. Your files will stay in `~/workshop/week11/germline`.

The `align/` folder is the biggest thing in there (about 3 GB). If you're tight on space, you can delete it, but hang on to it for now: we'll use these bams (including the tumor) again.

---

Adapted from the [PMBio Germline calling exercises](https://pmbio.org/module-04-germline/0004/01/01/Intro_to_Germline_WGS_and_Exome_Variant_Analysis/). See full sources and citations there.
