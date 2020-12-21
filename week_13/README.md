### Kallisto mini lecture

If you would like a refresher on Kallisto, we have made a [mini lecture](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2020/mini/RNASeq_MiniLecture_04_01_AlignmentFreeKallisto.pdf) briefly covering the topic.
We have also made a mini lecture  describing the differences between [alignment, assembly, and pseudoalignment](https://github.com/griffithlab/rnabio.org/blob/master/assets/lectures/cshl/2020/mini/RNASeq_MiniLecture_02_02_Alignment_vs_Assembly_vs_Kmer.pdf).


***

For more information on Kallisto, refer to the [Kallisto project page](https://pachterlab.github.io/kallisto/about.html), the [Kallisto manual page](https://pachterlab.github.io/kallisto/about.html) and the [Kallisto manuscript](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html).

***


### Build a Kallisto transcriptome index
Remember that Kallisto does not perform alignment or use a reference genome sequence. Instead it performs pseudoalignment to determine the compatibility of reads with targets (transcript sequences in this case). However, similar to alignment algorithms like Tophat or STAR, Kallisto requires an index to assess this compatibility efficiently and quickly.

```bash
cd ~/rnaseq
mkdir kallisto
cd kallisto
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto index --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index /docker_workspace/refs/chr22_ERCC92_transcripts.clean.fa
```
***

### Generate abundance estimates for all samples using Kallisto
As we did with `StringTie` and `HT-Seq` we will generate transcript abundances for each of our demonstration samples using `Kallisto`.

```bash
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/UHR_Rep1_ERCC-Mix1 --threads=4 --plaintext /docker_workspace/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/UHR_Rep2_ERCC-Mix1 --threads=4 --plaintext /docker_workspace/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/UHR_Rep3_ERCC-Mix1 --threads=4 --plaintext /docker_workspace/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz

docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/HBR_Rep1_ERCC-Mix2 --threads=4 --plaintext /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/HBR_Rep2_ERCC-Mix2 --threads=4 --plaintext /docker_workspace/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
docker run -v ~/rnaseq:/docker_workspace zlskidmore/kallisto:0.46.0 kallisto quant --rf-stranded --index=/docker_workspace/kallisto/chr22_ERCC92_transcripts_kallisto_index --output-dir=/docker_workspace/kallisto/HBR_Rep3_ERCC-Mix2 --threads=4 --plaintext /docker_workspace/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz /docker_workspace/data/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz
```

Create a single TSV file that has the TPM abundance estimates for all six samples.

```bash
paste */abundance.tsv | cut -f 1,2,5,10,15,20,25,30 > transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id\tlength$_\n"' > header.tsv
cat header.tsv transcript_tpms_all_samples.tsv | grep -v "tpm" > transcript_tpms_all_samples.tsv2
mv transcript_tpms_all_samples.tsv2 transcript_tpms_all_samples.tsv
rm -f header.tsv

```

Take a look at the final kallisto result file we created:

```bash
head transcript_tpms_all_samples.tsv
tail transcript_tpms_all_samples.tsv

```

***