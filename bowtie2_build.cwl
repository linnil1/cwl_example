# cwltool --user-space-docker-cmd=podman bowtie2_build.cwl test_var.yml
cwlVersion: v1.2
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0

requirements:
  InlineJavascriptRequirement: {}

class: CommandLineTool
baseCommand: bowtie2-build

arguments:
  - position: 2
    valueFrom: $( inputs.reference_fasta.nameroot )

inputs:
  reference_fasta:
    type: File
    inputBinding:
      position: 1

  threads:
    type: int
    inputBinding:
      prefix: --threads
    default: 16

outputs:
  bowtie2_index:
    type: File
    outputBinding:
      glob: $(inputs.reference_fasta.nameroot + ".1.bt2")
    secondaryFiles:
      - $(inputs.reference_fasta.nameroot + ".2.bt2")
      - $(inputs.reference_fasta.nameroot + ".3.bt2")
      - $(inputs.reference_fasta.nameroot + ".4.bt2")
      - $(inputs.reference_fasta.nameroot + ".rev.1.bt2")
      - $(inputs.reference_fasta.nameroot + ".rev.2.bt2")
