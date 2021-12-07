# cwltool --user-space-docker-cmd=podman pipeline.cwl chrx.yml
cwlVersion: v1.2
class: Workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  ResourceRequirement:
    coresMax: 16

inputs:
  fastq_1: File
  fastq_2: File
  reference_fasta: File

steps:
  bowtie_build:
    run: bowtie2_build.cwl
    in:
      reference_fasta: reference_fasta
    out: [bowtie2_index]

  bowtie_map:
    run: bowtie2_map.cwl
    in:
      reference_index: bowtie_build/bowtie2_index
      fastq_1: fastq_1
      fastq_2: fastq_2
    out: [bowtie_output_sam_tmp]

  samsort:
    run: samsort.cwl
    in:
      unsort_bam: bowtie_map/bowtie_output_sam_tmp
    out: [bowtie_output_bam]

  samindex:
    run: samindex.cwl
    in:
      sorted_bam: samsort/bowtie_output_bam
    out: [bowtie_output_index]

outputs:
  bowtie2_index:
    type: File
    outputSource: bowtie_build/bowtie2_index
  bowtie_output_bam:
    type: File
    outputSource: samsort/bowtie_output_bam
  bowtie_output_index:
    type: File
    outputSource: samindex/bowtie_output_index
