# cwltool --user-space-docker-cmd=podman samsort.cwl test_var.yml
cwlVersion: v1.2
hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.10--h9402c20_2

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMax: 16

class: CommandLineTool
baseCommand:
- samtools 
- sort

arguments:
  # xxx.sam -> xxx.bam
  - valueFrom: $( inputs.unsort_bam.nameroot + ".bam" )
    prefix: -o

inputs:
  unsort_bam:
    type: File
    inputBinding:
      position: 1

  threads:
    type: int
    inputBinding:
      prefix: "-@"
    default: 8

outputs:
  bowtie_output_bam:
    type: File
    outputBinding:
      glob: $( inputs.unsort_bam.nameroot + ".bam" )
