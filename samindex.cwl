# cwltool --user-space-docker-cmd=podman samindex.cwl test_var.yml
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
- index

arguments:
  # xxx.bam -> xxx.bam.bai
  - valueFrom: $( inputs.sorted_bam.basename + ".bai" )
    position: 2

inputs:
  sorted_bam:
    type: File
    inputBinding:
      position: 1

  threads:
    type: int
    inputBinding:
      prefix: "-@"
    default: 8

outputs:
  bowtie_output_index:
    type: File
    outputBinding:
      glob: $( inputs.sorted_bam.basename + ".bai" )
