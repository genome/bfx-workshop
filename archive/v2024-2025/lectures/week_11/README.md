## Week 11: Parsing and Filtering VCFs with Python

This week's lecture will be conducted on your local machine using Jupyter Notebook. Please ensure that you have Jupyter Notebook installed. You can reference the [the week 1 lecture notes](../week_01/bfx_workshop_01_overview.ipynb) for instructions.

To follow along with the lecture in real time, download the [notebook for this lecture](../week_11/python_vcf_parsing_and_filtering.ipynb) either by pulling the latest version of the bfx-workshop GitHub repo or downloading the raw ipynb file from GitHub. Then please run the commands listed in the "Setup" section at the top of the notebook. If you run the commands directly in your terminal and not from the notebook itself, do not include the `!` at the beginning of each command.

View the [lecture recording](https://wustl.box.com/s/ic7fn2ttipq5hl8v3qz00fhwfe0kyoo5)

View the [lecture slides](vcfs_vep_annotation.pdf)

## Homework
- Download the normal bams from here: gs://analysis-workflows-example-data/vcf_parsing_files/normal.bam, gs://analysis-workflows-example-data/vcf_parsing_files/normal.bam.bai
- Calculate readcounts for the normal sample and add them to the mutect.filtered.decomposed.pass_vaf_filtered.vcf.gz file
- Filter the output VCF from the above step to exclude variants where the normal VAF is higher than 0.02
- Amend the TSV-creation script to also output the normal VAF on top of the other information

  For credit studnets - send answers to Jenny
