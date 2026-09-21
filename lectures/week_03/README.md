## Week 3: DNA Sequence generation and QC

This week, we cover some basics of how sequence is generated, delve into the details of data and formats, and talk through basic QC of your sequence.

- [Lecture Recording](https://wustl.box.com/s/h89fnuqxid6nxmch25a3m5ttc7icd26i)
- [Slides](week03.pdf)

### Assignment for this week

1. Complete all of the questions in the [Sequence Data and Containers](sequence-data-and-containers.md) exercise.

2. Earlier we saw that `samtools` isn't available on the cluster. Find a Docker image that contains samtools
   (hint: search [Docker Hub](https://hub.docker.com/) for `samtools`), then click on the `tags` page to find the image name and version. Launch an interactive job with it, and run
   `samtools --version`. Then type `exit` to leave the job.

**For-credit students:** Send the following to John as proof of completion:

- a screenshot of one of your FastQC HTML reports

- the command you used to launch the samtools container in step 2 and it's output.

