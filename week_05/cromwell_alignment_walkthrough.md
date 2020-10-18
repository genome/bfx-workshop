Throughout this walkthrough, "GROUP" is the name of your compute group in compute1.  (If you run `groups`, specifically this is the part of a group that starts with `compute-` without that prefix included.)

To launch a cromwell workflow, start by launching an LSF job in a docker container containing Cromwell.
```bash
LSF_DOCKER_VOLUMES="$HOME:$HOME /scratch1/fs1/$GROUP:/scratch1/fs1/$GROUP /storage1/fs1/$GROUP/Active:/storage1/fs1/$GROUP/Active /storage1/fs1/bga/Active:/storage1/fs1/bga/Active" \
bsub -G compute-$GROUP -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-8)' -M 8000M -R 'select[mem>8000M] rusage[mem=8000M]' -Is -q general-interactive /bin/bash -l
```

To work properly on compute1, Cromwell needs to use a directory in the scratch space to run, so we'll make one now.
```bash
mkdir -p /scratch1/fs1/$GROUP/$USER/bfx-workshop-week5-cromwell
cd /scratch1/fs1/$GROUP/$USER/bfx-workshop-week5-cromwell
```

Then it's a matter of supplying the proper files to Cromwell to launch the workflow.  It needs a few things to run:
- a configuration file
- a file listing the inputs to the workflow
- the main file for the workflow itself

We can copy an example configuration file and edit it for our use.
```bash
cp /storage1/fs1/bga/Active/shared/tbmooney/bfx-workshop-20201019/cromwell.config .
```
<details>
  <summary>Click for cromwell.config</summary>
  
```hocon
include required(classpath("application"))

backend {
  default = "LSF"
  providers {
    LSF {
      actor-factory = "cromwell.backend.impl.sfs.config.ConfigBackendLifecycleActorFactory"
      config {
        runtime-attributes = """
        Int cpu = 1
        Int memory_mb = 4096
        String? docker
        """

        submit-docker = """
        LSF_DOCKER_VOLUMES='${cwd}:${docker_cwd} /scratch1/fs1/GROUP:/scratch1/fs1/GROUP /storage1/fs1/bga/Active/shared:/storage1/fs1/bga/Active/shared' \
        LSF_DOCKER_PRESERVE_ENVIRONMENT=false \
        bsub \
        -J ${job_name} \
        -cwd ${cwd} \
        -o /dev/null \
        -e cromwell-workflow-logs/cromwell-%J.err \
        -q 'general' \
        -g '/USER/cromwell-workers' \
        -G 'compute-GROUP' \
        -a "docker0(${docker})" \
        -M ${memory_mb}M \
        -n ${cpu} \
        -R "span[hosts=1] select[mem>${memory_mb}M] rusage[mem=${memory_mb}M]" \
        /bin/bash ${script}
        """

        kill = "bkill ${job_id}"
        docker-kill = "bkill ${job_id}"
        check-alive = "bjobs -noheader -o stat ${job_id} | /bin/grep 'PEND\\|RUN'"
        job-id-regex = "Job <(\\d+)>.*"
      }
    }
  }
}
```

</details>

Now it needs to be edited to be tailored to our user and group.  Replace USER in this file with your username and GROUP with your compute group.

For the inputs, an example input file can be found at `/storage1/fs1/bga/Active/shared/tbmooney/bfx-workshop-20201019/inputs.yaml`.  It can be used unmodified.

<details>
<summary>Click for inputs.yaml</summary>
  
```yaml
---

bait_intervals:
  class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_bait.interval_list

target_intervals:
  class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_target.interval_list

sequence:
  - sequence:
      bam:
        class: File
        path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499331.bam
    readgroup: "@RG\tID:2895499331\tPU:H7HY2CCXX.3\tSM:H_NJ-HCC1395-HCC1395_BL\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tPL:Illumina\tCN:WUGSC"
  - sequence:
      bam:
        class: File
        path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/unaligned_subset_bams/normal/2895499399.bam
    readgroup: "@RG\tID:2895499399\tPU:H7HY2CCXX.4\tSM:H_NJ-HCC1395-HCC1395_BL\tLB:H_NJ-HCC1395-HCC1395_BL-lg21-lib1\tPL:Illumina\tCN:WUGSC"

bqsr_known_sites:
- class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_known_indels.vcf.gz
- class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_mills.vcf.gz
- dbsnp_vcf:
  class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_dbsnp.vcf.gz

omni_vcf:
  class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_omni.vcf.gz

picard_metric_accumulation_level: LIBRARY

reference:
  class: File
  path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes.fa

bqsr_intervals:
- chr6
- chr17

per_base_intervals:
- file:
    class: File
    path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_target.interval_list
  label: clinvar

per_target_intervals:
- file:
    class: File
    path: /storage1/fs1/bga/Active/shared/analysis-workflows-example-data/somatic_inputs/hla_and_brca_genes_target.interval_list
  label: acmg_genes

summary_intervals: []
```

</details>

For the workflow itself, pull down the latest version of the analysis-workflows repository from github at https://github.com/genome/analysis-workflows/.
```bash
git clone https://github.com/genome/analysis-workflows.git
```

Now we can run the workflow by putting all the pieces together:
```bash
/usr/bin/java -Dconfig.file=cromwell.config -jar /opt/cromwell.jar run -t cwl -i /storage1/fs1/bga/Active/shared/tbmooney/bfx-workshop-20201019/inputs.yaml analysis-workflows/definitions/pipelines/alignment_exome.cwl
```

It will create two directories, `cromwell-executions`, and `cromwell-workflow-logs`.  One holds the work in progress and the results; the other holds the log files.  When the workflow finishes the output should include a block of JSON that points to where to find the output files.  Notably Cromwell does not place them all the in the current directory for you.  The end of the main log from Cromwell will list the paths to the output files.  For example:

<details>
  <summary>Click for outputs.json</summary>
  
```json
{
  "outputs": {
    "alignment_exome.cwl.per_base_coverage_metrics": [{
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_detailed_hs_metrics/hs_metrics.cwl/c4ae069e-2807-435d-9ca1-dd72b1b832d1/call-collect_per_base_hs_metrics/shard-0/execution/final.base-clinvar-PerBaseCoverage.txt",
      "size": 264508,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    }],
    "alignment_exome.cwl.summary_hs_metrics": [],
    "alignment_exome.cwl.per_base_hs_metrics": [{
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_detailed_hs_metrics/hs_metrics.cwl/c4ae069e-2807-435d-9ca1-dd72b1b832d1/call-collect_per_base_hs_metrics/shard-0/execution/final.base-clinvar-HsMetrics.txt",
      "size": 5547,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    }],
    "alignment_exome.cwl.hs_metrics": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_roi_hs_metrics/execution/final.roi-HsMetrics.txt",
      "size": 4995,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.alignment_summary_metrics": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_alignment_summary_metrics/execution/final.AlignmentSummaryMetrics.txt",
      "size": 4204,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.verify_bam_id_metrics": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-verify_bam_id/execution/final.VerifyBamId.selfSM",
      "size": 229,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.per_target_coverage_metrics": [{
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_detailed_hs_metrics/hs_metrics.cwl/c4ae069e-2807-435d-9ca1-dd72b1b832d1/call-collect_per_target_hs_metrics/shard-0/execution/final.target-acmg_genes-PerTargetCoverage.txt",
      "size": 4937,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    }],
    "alignment_exome.cwl.bam": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-alignment/sequence_to_bqsr.cwl/724f4e63-b233-4c06-9f23-cab10687966d/call-index_bam/execution/final.bam",
      "size": 641787,
      "secondaryFiles": [{
        "format": null,
        "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-alignment/sequence_to_bqsr.cwl/724f4e63-b233-4c06-9f23-cab10687966d/call-index_bam/execution/final.bam.bai",
        "size": null,
        "secondaryFiles": [],
        "contents": null,
        "checksum": null,
        "class": "File"
      }, {
        "format": null,
        "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-alignment/sequence_to_bqsr.cwl/724f4e63-b233-4c06-9f23-cab10687966d/call-index_bam/execution/final.bai",
        "size": null,
        "secondaryFiles": [],
        "contents": null,
        "checksum": null,
        "class": "File"
      }],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.verify_bam_id_depth": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-verify_bam_id/execution/final.VerifyBamId.depthSM",
      "size": 544,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.mark_duplicates_metrics": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-alignment/sequence_to_bqsr.cwl/724f4e63-b233-4c06-9f23-cab10687966d/call-mark_duplicates_and_sort/execution/final.merged.NameSorted.mark_dups_metrics.txt",
      "size": 3045,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.flagstats": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-samtools_flagstat/execution/final.bam.flagstat",
      "size": 402,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.insert_size_metrics": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_insert_size_metrics/execution/final.InsertSizeMetrics.txt",
      "size": 5286,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    },
    "alignment_exome.cwl.per_target_hs_metrics": [{
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_detailed_hs_metrics/hs_metrics.cwl/c4ae069e-2807-435d-9ca1-dd72b1b832d1/call-collect_per_target_hs_metrics/shard-0/execution/final.target-acmg_genes-HsMetrics.txt",
      "size": 5571,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    }],
    "alignment_exome.cwl.insert_size_histogram": {
      "format": null,
      "location": "/scratch1/fs1/GROUP/USER/bfx-workshop-week5-cromwell/cromwell-executions/alignment_exome.cwl/1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b/call-qc/qc_exome.cwl/e9ee5780-959d-453b-9c8a-3f1eef78e0f6/call-collect_insert_size_metrics/execution/final.InsertSizeHistogram.pdf",
      "size": 15046,
      "secondaryFiles": [],
      "contents": null,
      "checksum": null,
      "class": "File"
    }
  },
  "id": "1021e5a7-1df7-42d8-bcd8-978a3c7b8c9b"
}
```

</details>
