version 1.0

## Portions Copyright Broad Institute, 2018
##
## This WDL pipeline implements QC in human whole-genome or exome/targeted sequencing data.
##
## Requirements/expectations
## - Human paired-end sequencing data in aligned BAM or CRAM format
## - Input BAM/CRAM files must additionally comply with the following requirements:
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL open source code license (BSD-3).
## Full license text at https://github.com/openwdl/wdl/blob/master/LICENSE
## Note however that the programs it calls may be subject to different licenses. 
## Users are responsible for checking that they are authorized to run all programs before running this script.
## - [Picard](https://broadinstitute.github.io/picard/)
## - [VerifyBamID2](https://github.com/Griffan/VerifyBamID)

# Git URL import
import "tasks/Qc.wdl" as QC

# WORKFLOW DEFINITION
workflow SingleSampleQc {
  input {
    File input_bam
    File ref_cache
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String base_name
    Int preemptible_tries
    File coverage_interval_list
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    Boolean is_wgs
    Boolean? is_outlier_data

    File evaluation_thresholds
  }

  # Not overridable:
  Int read_length = 250
  
  # Generate a BAM or CRAM index
  call QC.BuildBamIndex as BuildBamIndex {
    input:
      input_bam = input_bam,
      base_name = base_name,
      ref_cache = ref_cache,
      preemptible_tries = preemptible_tries
  }
  
  # Collect BAM/CRAM index stats
  call QC.BamIndexStats as BamIndexStats {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      preemptible_tries = preemptible_tries
  }

  call QC.RxIdentifier as RxIdentifier {
    input:
      idxstats = BamIndexStats.idxstats,
      sample_id = base_name,
      preemptible_tries = preemptible_tries
  }

  # Validate the BAM or CRAM file
  call QC.ValidateSamFile as ValidateSamFile {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      report_filename = base_name + ".validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ignore = ["MISSING_TAG_NM"],
      is_outlier_data = is_outlier_data,
      preemptible_tries = preemptible_tries
   }

  # generate a md5
  call QC.CalculateChecksum as CalculateChecksum {
    input:
      input_bam = BuildBamIndex.bam,
      preemptible_tries = preemptible_tries
  }

  # QC the final BAM some more (no such thing as too much QC)
  call QC.CollectAggregationMetrics as CollectAggregationMetrics {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      base_name = base_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
  }

  # QC the BAM sequence yield
  call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      metrics_filename = base_name + ".quality_yield_metrics",
      preemptible_tries = preemptible_tries
  }
  
  if (is_wgs) {
    # QC the sample raw WGS metrics since this IS WGS data
    call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
      input:
        input_bam = BuildBamIndex.bam,
        input_bam_index = BuildBamIndex.bam_index,
        metrics_filename = base_name + ".raw_wgs_metrics",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        wgs_coverage_interval_list = coverage_interval_list,
        read_length = read_length,
        preemptible_tries = preemptible_tries
    }
  }
  
  if (!is_wgs) {
    # QC the sample Hs/WES metrics since this IS NOT WGS
    call QC.CollectHsMetrics as CollectHsMetrics {
      input:
        input_bam = BuildBamIndex.bam,
        input_bam_index = BuildBamIndex.bam_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        metrics_filename = base_name + ".hs_metrics",
        target_interval_list = coverage_interval_list,
        bait_interval_list = coverage_interval_list,
        preemptible_tries = preemptible_tries
    }
  }
  
  # Estimate level of cross-sample contamination
  call QC.CheckContamination as CheckContamination {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = base_name + ".verify_bam_id",
      preemptible_tries = preemptible_tries,
  }

  # Calculate the duplication rate since MarkDuplicates was already performed
  call QC.CollectDuplicateMetrics as CollectDuplicateMetrics {
    input:
      input_bam = BuildBamIndex.bam,
      input_bam_index = BuildBamIndex.bam_index,
      output_bam_prefix = base_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
  }

  call QC.EvaluateMetrics as EvaluateMetrics {
    input:
      thresholds = evaluation_thresholds,
      alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics,
      duplication_metrics = CollectDuplicateMetrics.duplication_metrics,
      insert_size_metrics = CollectAggregationMetrics.insert_size_metrics,
      quality_yield_metrics = CollectQualityYieldMetrics.metrics,
      contamination_metrics = CheckContamination.metrics,
      hs_metrics = CollectHsMetrics.hs_metrics,
      wgs_metrics = CollectRawWgsMetrics.metrics,
      preemptible_tries = preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {

    File validation_report = ValidateSamFile.report

    File alignment_summary_metrics_file = CollectAggregationMetrics.alignment_summary_metrics_file
    String pct_chimeras = CollectAggregationMetrics.pct_chimeras
    String read1_pf_mismatch_rate = CollectAggregationMetrics.read1_pf_mismatch_rate
    String read2_pf_mismatch_rate = CollectAggregationMetrics.read2_pf_mismatch_rate

    File bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics

    File insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File insert_size_metrics_file = CollectAggregationMetrics.insert_size_metrics_file
    String median_insert_size = CollectAggregationMetrics.median_insert_size
    String median_absolute_deviation = CollectAggregationMetrics.median_absolute_deviation

    File pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics
    File error_summary_metrics = CollectAggregationMetrics.error_summary_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplication_metrics_file = CollectDuplicateMetrics.duplication_metrics_file
    String percent_duplication = CollectDuplicateMetrics.percent_duplication

    File quality_yield_metrics = CollectQualityYieldMetrics.metrics_file
    String q20_bases = CollectQualityYieldMetrics.q20_bases
    String pf_q20_bases = CollectQualityYieldMetrics.pf_q20_bases
    String q30_bases = CollectQualityYieldMetrics.q30_bases
    String pf_q30_bases = CollectQualityYieldMetrics.pf_q30_bases

    File? raw_wgs_metrics = CollectRawWgsMetrics.metrics_file
    String? mean_coverage = CollectRawWgsMetrics.mean_coverage
    String? pct_10x = CollectRawWgsMetrics.pct_10x
    String? pct_20x = CollectRawWgsMetrics.pct_20x
    String? pct_30x = CollectRawWgsMetrics.pct_30x

    File? hs_metrics = CollectHsMetrics.hs_metrics_file
    String? mean_target_coverage = CollectHsMetrics.mean_target_coverage
    String? pct_target_bases_10x = CollectHsMetrics.pct_target_bases_10x
    String? pct_target_bases_20x = CollectHsMetrics.pct_target_bases_20x
    String? pct_target_bases_30x = CollectHsMetrics.pct_target_bases_30x

    File input_bam_md5 = CalculateChecksum.md5
    String input_bam_md5_hash = CalculateChecksum.md5_hash
    File input_bam_index = BuildBamIndex.bam_index
    File input_bam_idxstats = BamIndexStats.idxstats
    File input_bam_rx_result = RxIdentifier.rx_result
    String input_bam_rx_value = RxIdentifier.rx_value

    File evaluated_metrics_file = EvaluateMetrics.evaluated_metrics_file
    Map[String, String] evaluated_metrics = EvaluateMetrics.evaluated_metrics
    String overall_evaluation = EvaluateMetrics.overall_evaluation
  }
}
