version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for QC of human whole-genome or exome sequencing data.
## ## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  String metric_object_file = "~{metrics_filename}.only_metrics"

  command <<<
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OQ=true \
      OUTPUT=~{metrics_filename}

    grep -v '#' ~{metrics_filename} | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..$#k) { say join("\t", $k[$_], $v[$_]); }' > "~{metric_object_file}"
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3 GiB"
    preemptible: preemptible_tries
  }
  output {
    File metrics_file = "~{metrics_filename}"
    Map[String, String] metrics = read_map(metric_object_file)
    String q20_bases = metrics["Q20_BASES"]
    String pf_q20_bases = metrics["PF_Q20_BASES"]
    String q30_bases = metrics["Q30_BASES"]
    String pf_q30_bases = metrics["PF_Q30_BASES"]
  }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
  input {
    File input_bam
    File input_bam_index
    String base_name
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Boolean collect_gc_bias_metrics = true
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  String alignment_summary_metric_object_file = "~{base_name}.alignment_summary_metrics.only_metrics"
  String insert_size_metric_object_file = "~{base_name}.insert_size_metrics.only_metrics"

  command <<<
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{base_name}.gc_bias.detail_metrics \
      ~{base_name}.gc_bias.pdf \
      ~{base_name}.gc_bias.summary_metrics \
      ~{base_name}.insert_size_metrics \
      ~{base_name}.insert_size_histogram.pdf

    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{base_name} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=QualityScoreDistribution \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY

    grep -v '#' '~{base_name}.alignment_summary_metrics' | grep '.\+' | perl -E 'my ($keys, @values) = <>; chomp $keys; chomp for @values; my @k = split("\t", $keys); for my $values (@values[0..2]) { my @v = split("\t", $values); my $category = $v[0]; for(0..($#k-2)) { say join("\t", join("-",$k[$_], $category), $v[$_]); }; }' > ~{alignment_summary_metric_object_file}

    grep -v '#' '~{base_name}.insert_size_metrics' | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..($#k-2)) { say join("\t", $k[$_], $v[$_]); }' > ~{insert_size_metric_object_file}


  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    memory: "7 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File alignment_summary_metrics_file = "~{base_name}.alignment_summary_metrics"
    Map[String, String] alignment_summary_metrics = read_map(alignment_summary_metric_object_file)
    String pct_chimeras = alignment_summary_metrics["PCT_CHIMERAS-PAIR"]
    String read1_pf_mismatch_rate = alignment_summary_metrics["PF_MISMATCH_RATE-FIRST_OF_PAIR"]
    String read2_pf_mismatch_rate = alignment_summary_metrics["PF_MISMATCH_RATE-SECOND_OF_PAIR"]

    File bait_bias_detail_metrics = "~{base_name}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "~{base_name}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "~{base_name}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{base_name}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{base_name}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "~{base_name}.insert_size_histogram.pdf"
    File insert_size_metrics_file = "~{base_name}.insert_size_metrics"
    Map[String, String] insert_size_metrics = read_map(insert_size_metric_object_file)
    String median_insert_size = insert_size_metrics["MEDIAN_INSERT_SIZE"]
    String median_absolute_deviation = insert_size_metrics["MEDIAN_ABSOLUTE_DEVIATION"]

    File pre_adapter_detail_metrics = "~{base_name}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "~{base_name}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "~{base_name}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{base_name}.quality_distribution_metrics"
    File error_summary_metrics = "~{base_name}.error_summary_metrics"
  }
}

task ValidateSamFile {
  input {
    File input_bam
    File input_bam_index
    String report_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Array[String]? ignore
    Boolean? is_outlier_data
    Int preemptible_tries
    Int memory_multiplier = 1
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  Int memory_size = ceil(7 * memory_multiplier)
  Int java_memory_size = (memory_size - 1) * 1000

  command {
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=SUMMARY \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File report = "~{report_filename}"
  }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_fasta_index
    Int read_length
    Int preemptible_tries
    Int memory_multiplier = 1
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  Int memory_size = ceil((if (disk_size < 110) then 5 else 7) * memory_multiplier)
  String java_memory_size = (memory_size - 1) * 1000

  String wgs_metric_object_file = '~{metrics_filename}.only_metrics'

  command <<<
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      MINIMUM_MAPPING_QUALITY=0 \
      MINIMUM_BASE_QUALITY=3 \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}

    sed -i.original 's/picard.analysis.WgsMetrics/picard.analysis.CollectWgsMetrics\$WgsMetrics/' ~{metrics_filename}

    grep -v '#' '~{metrics_filename}' | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..$#k) { say join("\t", $k[$_], $v[$_]); }' > ~{wgs_metric_object_file}

  >>>
  runtime {
    # Using older image due to: https://github.com/broadinstitute/picard/issues/1402
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.20.4"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File metrics_file = "~{metrics_filename}"
    Map[String, String] metrics = read_map(wgs_metric_object_file)
    String mean_coverage = metrics["MEAN_COVERAGE"]
    String pct_10x = metrics["PCT_10X"]
    String pct_20x = metrics["PCT_20X"]
    String pct_30x = metrics["PCT_30X"]
  }
}

task CollectHsMetrics {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String metrics_filename
    File target_interval_list
    File bait_interval_list
    Int preemptible_tries
    Int memory_multiplier = 1
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20
  # Try to fit the input bam into memory, within reason.
  Int rounded_bam_size = ceil(size(input_bam, "GiB") + 0.5)
  Int rounded_memory_size = ceil((if (rounded_bam_size > 10) then 10 else rounded_bam_size) * memory_multiplier)
  Int memory_size = if rounded_memory_size < 7 then 7 else rounded_memory_size
  Int java_memory_size = (memory_size - 1) * 1000

  String hs_metric_object_file = "~{metrics_filename}.only_metrics"

  # There are probably more metrics we want to generate with this tool
  command <<<
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectHsMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      VALIDATION_STRINGENCY=SILENT \
      TARGET_INTERVALS=~{target_interval_list} \
      BAIT_INTERVALS=~{bait_interval_list} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      OUTPUT=~{metrics_filename}

grep -v '#' '~{metrics_filename}' | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..$#k) { say join("\t", $k[$_], $v[$_]); }' > ~{hs_metric_object_file}

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File hs_metrics_file = metrics_filename
    Map[String, String] hs_metrics = read_map(hs_metric_object_file)
    String mean_target_coverage = hs_metrics["MEAN_TARGET_COVERAGE"]
    String pct_target_bases_10x = hs_metrics["PCT_TARGET_BASES_10X"]
    String pct_target_bases_20x = hs_metrics["PCT_TARGET_BASES_20X"]
    String pct_target_bases_30x = hs_metrics["PCT_TARGET_BASES_30X"]
  }
}

# Collect duplicate metrics
task CollectDuplicateMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  String duplication_metric_object_file = '~{output_bam_prefix}.duplication_metrics.metrics_only'

  command <<<
    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectDuplicateMetrics \
      METRICS_FILE=~{output_bam_prefix}.duplication_metrics \
      INPUT=~{input_bam} \
      ASSUME_SORTED=true \
      REFERENCE_SEQUENCE=~{ref_fasta}

    grep -v '#' '~{output_bam_prefix}.duplication_metrics' | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..$#k) { say join("\t", $k[$_], $v[$_]); }'      > ~{duplication_metric_object_file}

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    memory: "7 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File duplication_metrics_file = "~{output_bam_prefix}.duplication_metrics"
    Map[String, String] duplication_metrics = read_map(duplication_metric_object_file)
    String percent_duplication = duplication_metrics["PERCENT_DUPLICATION"]
  }
}

# Build BAM/CRAM index, this should run before all other QC tasks
task BuildBamIndex {
  input {
    File input_bam
    String base_name
    File ref_cache
    Int preemptible_tries
  }

  Float ref_size = size(ref_cache, "GiB") * 5 + 1.0
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20
  String bam_link = sub(basename(input_bam), basename(basename(input_bam, ".bam"), ".cram"), base_name )
  
  command {
    # Link the BAM/CRAM to a harmonized name and path
    cp ~{input_bam} ~{bam_link}

    # build the reference sequence cache
    tar -zxf ~{ref_cache}
    export REF_PATH=./cache/%2s/%2s/%s
    export REF_CACHE=./cache/%2s/%2s/%s

    # index the BAM or CRAM
    /usr/local/bin/samtools index ~{bam_link}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.2-1564506709"
    memory: "1 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File bam = bam_link
    File bam_index = sub(sub(bam_link, "bam$", "bam.bai"), "cram$", "cram.crai")
  }
}

# Collect BAM/CRAM index stats
# NOTE: Samtools idxstats is slow on CRAM inputs since the index does not contain the info needed to summarize per-chromosome alignments
task BamIndexStats {
  input {
    File input_bam
    File input_bam_index
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20
  String output_name = basename(input_bam)

  command {
    samtools idxstats ~{input_bam} > ~{output_name}.idxstats
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/samtools:1.0.0-1624651616"
    memory: "1 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File idxstats = "~{output_name}.idxstats"
  }
}

task RxIdentifier {
  input {
    File idxstats
    String sample_id
    Int preemptible_tries
  }

  Int disk_size = ceil(size(idxstats, "GiB")) + 20
  String output_name = basename(idxstats) + ".Rx"


# Rx_identifier
# # Based on the ratio of X chromosome-derived shotgun sequencing data to the autosomal coverage to establish the probability of an XX or XY karyotype for ancient samples.
# # Author: Chuan-Chao Wang
# # Contact: wang@shh.mpg.de, chuan-chao_wang@hms.harvard.edu; Department of Genetics, Harvard Medical School; Department of Archaeogenetics, Max Planck Institute for the Science of Human History.
# # Date: 30 Jan, 2016
# # https://doi.org/10.1371/journal.pone.0163019.s003
  command <<<
    R <<SCRIPT
        sample_id <- "~{sample_id}"

        idxstats <- read.table("~{idxstats}",header=F,nrows=24,row.names=1)
        c1 <- c(as.numeric(idxstats[,1]))
        c2 <- c(as.numeric(idxstats[,2]))
        total_ref <- sum(c1)
        total_map <- sum(c2)

        LM <- lm(c1~c2)

        Rt1 <- (idxstats[1,2]/total_map)/(idxstats[1,1]/total_ref)
        Rt2 <- (idxstats[2,2]/total_map)/(idxstats[2,1]/total_ref)
        Rt3 <- (idxstats[3,2]/total_map)/(idxstats[3,1]/total_ref)
        Rt4 <- (idxstats[4,2]/total_map)/(idxstats[4,1]/total_ref)
        Rt5 <- (idxstats[5,2]/total_map)/(idxstats[5,1]/total_ref)
        Rt6 <- (idxstats[6,2]/total_map)/(idxstats[6,1]/total_ref)
        Rt7 <- (idxstats[7,2]/total_map)/(idxstats[7,1]/total_ref)
        Rt8 <- (idxstats[8,2]/total_map)/(idxstats[8,1]/total_ref)
        Rt9 <- (idxstats[9,2]/total_map)/(idxstats[9,1]/total_ref)
        Rt10 <- (idxstats[10,2]/total_map)/(idxstats[10,1]/total_ref)
        Rt11 <- (idxstats[11,2]/total_map)/(idxstats[11,1]/total_ref)
        Rt12 <- (idxstats[12,2]/total_map)/(idxstats[12,1]/total_ref)
        Rt13 <- (idxstats[13,2]/total_map)/(idxstats[13,1]/total_ref)
        Rt14 <- (idxstats[14,2]/total_map)/(idxstats[14,1]/total_ref)
        Rt15 <- (idxstats[15,2]/total_map)/(idxstats[15,1]/total_ref)
        Rt16 <- (idxstats[16,2]/total_map)/(idxstats[16,1]/total_ref)
        Rt17 <- (idxstats[17,2]/total_map)/(idxstats[17,1]/total_ref)
        Rt18 <- (idxstats[18,2]/total_map)/(idxstats[18,1]/total_ref)
        Rt19 <- (idxstats[19,2]/total_map)/(idxstats[19,1]/total_ref)
        Rt20 <- (idxstats[20,2]/total_map)/(idxstats[20,1]/total_ref)
        Rt21 <- (idxstats[21,2]/total_map)/(idxstats[21,1]/total_ref)
        Rt22 <- (idxstats[22,2]/total_map)/(idxstats[22,1]/total_ref)
        Rt23 <- (idxstats[23,2]/total_map)/(idxstats[23,1]/total_ref)
        Rt24 <- (idxstats[24,2]/total_map)/(idxstats[24,1]/total_ref)

        tot <- c(Rt23/Rt1,Rt23/Rt2,Rt23/Rt3,Rt23/Rt4,Rt23/Rt5,Rt23/Rt6,Rt23/Rt7,Rt23/Rt8,Rt23/Rt9,Rt23/Rt10,Rt23/Rt11,Rt23/Rt12,Rt23/Rt13,Rt23/Rt14,Rt23/Rt15,Rt23/Rt16,Rt23/Rt17,Rt23/Rt18,Rt23/Rt19,Rt23/Rt20,Rt23/Rt21,Rt23/Rt22)
        Rx <- mean(tot)
        confinterval <- 1.96*(sd(tot)/sqrt(22))
        CI1 <- Rx-confinterval
        CI2 <- Rx+confinterval
        CI <- paste0(CI1,'-',CI2)

        result <- 'NA'
        if (CI1 > 0.8) {
            result <- 'XX'
        } else if (CI2 < 0.6) {
            result <- 'XY'
        } else if (CI1 > 0.6 & CI2 > 0.8) {
            result <- 'consistent with XX but not XY'
        } else if (CI1 < 0.6 & CI2 < 0.8) {
            result <- 'consistent with XY but not XX'
        } else {
            result <- 'Not Assigned'
        }
        write(paste(sample_id,total_map,Rx,confinterval,CI,result, sep = "\t"), file = "~{output_name}", append=FALSE)
    SCRIPT
  >>>

  runtime {
    docker: "us.gcr.io/anvil-gcr-public/anvil-rstudio-bioconductor:3.14.0"
    memory: "1 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File rx_result = "~{output_name}"
    String rx_value = read_tsv(output_name)[0][5]
  }
}


# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Boolean disable_sanity_check = false
  }

  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"]))
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "4 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c1cba76e979904eb69c31520a0d7f5be63c72253-1553018888"
    cpu: "2"
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
    Map[String, String] metrics = { "FREEMIX": read_string(stdout()) }
  }
}

task CalculateChecksum {
  input {
    File input_bam
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  String output_name = basename(input_bam)

  command {
    md5sum ~{input_bam} > ~{output_name}.md5
  }
  runtime {
    preemptible: preemptible_tries
    memory: "2 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
  }
  output {
    File md5 = "~{output_name}.md5"
    String md5_hash = sub(read_string(md5), " .+$", "")
  }
}

task EvaluateMetrics {
  input {
    File thresholds
    Map[String, String] alignment_summary_metrics
    Map[String, String] duplication_metrics
    Map[String, String] insert_size_metrics
    Map[String, String] quality_yield_metrics
    Map[String, String] contamination_metrics
    Map[String, String]? hs_metrics
    Map[String, String]? wgs_metrics

    Int preemptible_tries
  }

  Int disk_size = ceil(size(thresholds, "GiB")) + 20

  String evaluated_metrics_filename = 'evaluated_metrics.tsv'
  String overall_evaluation_file = 'overall_evaluation'

  command <<<
    perl -E '
      use strict; use warnings;

      open(my $threshold_fh, "<", "~{thresholds}");
      my @thresholds = <$threshold_fh>; chomp @thresholds;
      close($threshold_fh);

      my @metrics_files = @ARGV;
      my %metrics;

      for my $f (@metrics_files) {
        open(my $fh, "<", $f);
        my @data = <$fh>; chomp @data;
        close($fh);

        for my $metric (@data) {
          my ($name, $value) = split("\t", $metric);
          $metrics{$name} = $value;
        }
      }

      my %comparisons = (
        "="  => sub { return $_[0] == $_[1]; },
        "<"  => sub { return $_[0] <  $_[1]; },
        ">"  => sub { return $_[0] >  $_[1]; },
        "<=" => sub { return $_[0] <= $_[1]; },
        ">=" => sub { return $_[0] >= $_[1]; },
      );

      for my $test (@thresholds) {
        my ($name, $threshold_value, $operator) = split("\t", $test);

        my $comparator = $comparisons{$operator} or die "unknown operator $operator";
        my $metric_value = $metrics{$name};
        my $status = "n/a";
        if ($metric_value) {
          my $result = $comparator->($metric_value, $threshold_value);
          $status = $result? "PASS" : "FAIL";
        }
        say join("\t", $name, $status);
      }
    ' ~{write_map(alignment_summary_metrics)} \
      ~{write_map(duplication_metrics)} \
      ~{write_map(insert_size_metrics)} \
      ~{write_map(quality_yield_metrics)} \
      ~{write_map(contamination_metrics)} \
      ~{write_map(select_first([hs_metrics, wgs_metrics]))} \
      > ~{evaluated_metrics_filename}

    if [[ $(grep 'FAIL' ~{evaluated_metrics_filename}) ]]; then
      echo 'FAIL' > ~{overall_evaluation_file}
    elif [[ $(grep 'n/a' ~{evaluated_metrics_filename}) ]]; then
      echo 'PASS*' > ~{overall_evaluation_file}
    else
      echo 'PASS' > ~{overall_evaluation_file}
    fi

  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "2 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
  }
  output {
    File evaluated_metrics_file = "~{evaluated_metrics_filename}"
    Map[String, String] evaluated_metrics = read_map(evaluated_metrics_filename)
    String overall_evaluation = read_string(overall_evaluation_file)
  }
}
