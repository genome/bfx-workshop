## Week 18: When good experiments go bad

This week, we cover examples of common mistakes & pitfalls when analyzing genomic data. Learn about checking for errors and how to fix them.

- [Lecture Recording](https://wustl.box.com/s/ad3j9xo22gn5hb9pw1043ljjqgkmb1me)

- [Slides](week18_bad_experiments.pdf)

## Exercise

Often, diagnosing issues with sequencing data involves actually visualizing the read data in IGV.  For this week's exercise, we are giving you five sequencing runs, and you will need to infer the species and sequencing modality of each.  (The kind of thing you might have to do when untangling sample swaps!)  Think carefully about what kind of coverage you might expect from each, and what a mismatched species or sequencing type would do to your alignments.

- Each sample could be from Human or Mouse, and the data may have been produced via Exome sequencing, Whole Genome Sequencing, or RNA sequencing

- To save you time, each of these five samples (A, B, C, D, and E) have been pre-aligned with BWA-mem to both the hg38 (human) and mm10 (mouse) genomes. They can be accessed for each sample using this pattern: `https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/exp-fail-inputs/dataset_{A,B,C,D,E}-{hg38,mm10}-bwa-mem.bam`

- Example for Sample A, hg38:
`https://storage.googleapis.com/icts-precision-health-bfx-workshop-public-data/exp-fail-inputs/dataset_A-hg38-bwa-mem.bam`

- Though you could certainly download these somewhat large bams and then open them, IGV can load files directly from URLs, using `File > Open From URL`

- Hint: This data is only from a small subset of the genome, and the PPM1F/Ppm1f and PUM3/Pum3 genes are good places to start if you're looking for regions with sequence coverage.

- If you need IGV tips, you can refer back to [Week 5](https://github.com/genome/bfx-workshop/tree/master/lectures/week_05)

Students taking the course for credit, please send in two things:

- A list of samples and your species/seqtype assignments. 

- An IGV screenshot showing an example of one region and a note explaining why it was informative. (just one, not needed for each sample)


