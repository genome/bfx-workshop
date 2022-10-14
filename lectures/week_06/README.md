## Week 6: DNA Alignment Continued

This week, we cover sequence data alignment in more detail, including Base Quality Score Recalibration, SAM flags and marking of duplicate fragments (either from PCR or optical flowcell duplicates). We will then discuss streaming multiple steps together as well as an automated pipeline in the form of a Workflow. 

- [Lecture Recording]()

- [Slides](bfx_workshop_06_alignment.pdf)

- [Workflow Tutorial](bfx_workshop_06_alignment.md)

## Homework Assignment

Using Google Cloud Shell and a Cromwell Server hosted on Google, align a Normal HCC1395 sample and the Tumor HCC1395 sample just as we did in Week 4. This time the workflow includes all of the alignment steps in addition to some basic Exome Quality Control. When finished you will copy your BAMs and QC files locally to view them in IGV (for the BAM files) or an appropriate application (for QC outputs including text, PDF, and PNG files). 
