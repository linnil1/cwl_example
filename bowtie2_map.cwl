# cwltool --user-space-docker-cmd=podman bowtie2_map.cwl test_var.yml
cwlVersion: v1.2
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMax: 16

class: CommandLineTool
baseCommand: bowtie2

arguments:
  - valueFrom: --no-unal
  # xxx.read1.fq  -> xxx.bowtie.sam
  - valueFrom: $( inputs.fastq_1.basename.slice(0, -9) + ".bowtie.sam" )
    prefix: -S

inputs:
  fastq_1:
    type: File
    inputBinding:
      prefix: "-1"

  fastq_2:
    type: File
    inputBinding:
      prefix: "-2"

  reference_index:
    type: File
    secondaryFiles:
      - pattern: $( inputs.reference_index.basename.slice(0, -6) + ".2.bt2" )
      - pattern: $( inputs.reference_index.basename.slice(0, -6) + ".3.bt2" )
      - pattern: $( inputs.reference_index.basename.slice(0, -6) + ".4.bt2" )
      - pattern: $( inputs.reference_index.basename.slice(0, -6) + ".rev.1.bt2" )
      - pattern: $( inputs.reference_index.basename.slice(0, -6) + ".rev.2.bt2" )
    inputBinding:
      prefix: -x
      valueFrom: $( inputs.reference_index.path.slice(0, -6) )

  threads:
    type: int
    inputBinding:
      prefix: --threads
    default: 8
      
outputs:
  bowtie_output_sam_tmp:
    type: File
    outputBinding:
      glob: $( inputs.fastq_1.basename.slice(0, -9) + ".bowtie.sam" )
