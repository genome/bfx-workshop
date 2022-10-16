# Week 6
## DNA Alignment Continued

This week, we cover sequence data alignment in more detail, including Base Quality Score Recalibration, SAM flags and marking of duplicate fragments (either from PCR or optical flowcell duplicates). We will then discuss streaming multiple steps together as well as an automated pipeline in the form of a Workflow. 

- [Lecture Recording]()

- [Slides](bfx_workshop_06_alignment.pdf)

- [Workflow Tutorial](bfx_workshop_06_alignment.md)

## Homework Assignment
1. Please complete and/or troubleshoot any of the steps in the [tutorial](bfx_workshop_06_alignment.md) we were unable to complete in our 1-hour session.
HINT: Post questions to `bfx-workshop` Slack: https://ictsprecisionhealth.slack.com/archives/C040Q704WS2 
2. Assess the AlignmentSummaryMetrics. What percentage of the reads aligned? 
HINT: `PCT_PF_READS_ALIGNED`
3. Pull the HsMetrics (ie. Exome) coverage numbers to assess the coverage of the sequence data for the regions of interest (ie. interval lists) provided. 
HINT: Look for `MEAN_BAIT_COVERAGE` or `MEAN_TARGET_COVERAGE`. The values may be different depending on the interval list used. For this tutorial, most of the metrics should be similar in value.
4. View the histogram of InsertSizeMetrics  to see a distribution of fragment sizes for each library. What is the average insert size? 
HINT: Look for `MEAN_INSERT_SIZE`
5. Using IGV, visualize the alignments for both the Tumor and Normal sample together in separate tracks. See if you can visually identify any potential variants, ie. Single Nucleotide Variants (SNV), of interest. 
HINT: BRCA1 is one of the genes covered by the example data set. It's not easy or efficient to "call variants" using visualizations of reference alignments. 
Coming soon... Variant Detection!
