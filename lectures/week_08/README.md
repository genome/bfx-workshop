# Week 8

## RNAseq Part 2 - Quantification

- Lecture Recording

- [Abundance Estimation Slides](RNASeq_part2_bfx.pdf)

- [RNA-seq Bioinformatics home](https://rnabio.org/course) (rnabio.org)


## Homework Assignments

1) Complete the "Expression Analysis" portion of Module 3 - [Expression](https://rnabio.org/module-03-expression/0003/02/01/Expression/)

Reminder: this should all be done in the same docker image as last week, with your rnaseq folder mounted appropriately:

```
docker run -it -p 8080:8080 -v ~/workshop/rnabio-workspace:/workspace griffithlab/rnaseq-toolkit:latest
```

For-credit students: please find the TPM value of gene CECR7 in sample HBR_Rep1 and send that to John.
