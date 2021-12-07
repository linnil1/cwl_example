# Writing my own cwl

CWL stands for Common Workflow Language,
which maybe commonly used in the future 
because of its strict spec and the ability to cross platform.

CWL should written in Yaml format, due to cross-platofrm compatibility.
I think it is not that readable and not easy to write well.

The only flexibilty in CWL is that it allow to write arbitrary expression
e.g. `valueFrom: $( inputs.unsort_bam.nameroot + ".bam" )`, though it is javascript style.

## My example
I generate a fasta file and a pair-end fastq file, then running read-mapping by bowtie2.

In short, we can simlify our pipeline(4 + 1 cwl) into these ipython code (10 lines only).

``` python
index = "data/chrx"
fastq = "data/chrx.synthetic"
thread = 30

%alias dk docker run --rm -it -v $PWD:/app -w /app

%dk quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0 \
    bowtie2-build --threads {thread} {index}.fa {index}
%dk quay.io/biocontainers/bowtie2:2.4.4--py39hbb4e92a_0 \
    bowtie2       --threads {thread} -x {index} --no-unal -1 {fastq}.read1.fq -2 {fastq}.read2.fq -S {fastq}.bowtie.sam
fastq = fastq + ".bowtie"
%dk quay.io/biocontainers/samtools:1.10--h9402c20_2 \
    samtools sort         -@{thread} {fastq}.sam -o {fastq}.bam
%dk quay.io/biocontainers/samtools:1.10--h9402c20_2 \
    samtools index        -@{thread} {fastq}.bam    {fastq}.bam.bai
```

## Run
``` bash
pip3 install cwltool
rm *.bt2 *.sam *.bam *.bai

# step by step
cwltool --user-space-docker-cmd=podman bowtie2_build.cwl chrx_test.yml
cwltool --user-space-docker-cmd=podman bowtie2_map.cwl   chrx_test.yml
cwltool --user-space-docker-cmd=podman samsort.cwl       chrx_test.yml
cwltool --user-space-docker-cmd=podman samindex.cwl      chrx_test.yml

# Run all steps in once
cwltool --user-space-docker-cmd=podman pipeline.cwl      chrx.yml
```

## Run with airflow

Airflow is a visualization & DAG manage tool

``` bash
# install (cwl-airflow will install airflow too)
pip3 install cwl-airflow

# init airflow
export AIRFLOW_HOME=$PWD/airflow/
airflow db init
airflow users create -e linnil1 -f linnil1 -l linnil1 -r Admin -u linnil1 -p password

# add my dags and copy input files
mkdir airflow/dags/
sed -e "s|\$PWD|$PWD|g" airflow/bowtie2_cwl.py > airflow/dags/bowtie2_cwl.py
cp -r data/ airflow/cwl_inputs_folder

# run
airflow webserver

# Run this in another window
export AIRFLOW_HOME=$PWD/airflow/
airflow scheduler

# Run this in another window
export AIRFLOW_HOME=$PWD/airflow/
cwl-airflow init
cwl-airflow api

# Open browser localhost:8080
# uername: linnil1
# password: password
# Find bowtie2_cwl
# **Unpause DAG**
# **Trigger with config**
# Paste the output of this
python3 -c "import json; import yaml; print(json.dumps({'job': yaml.load(open('chrx.yml'))}))"
# {
# "job": {
#   "reference_fasta": {"class": "File", "path": "data/chrx.fa"},
#   "fastq_1": {"class": "File", "path": "data/chrx.synthetic.read1.fq"},
#   "fastq_2": {"class": "File", "path": "data/chrx.synthetic.read2.fq"},
#   "threads": 10}
# }                                                                          
```

## TODO?
* [x] Suffix name
* [x] airflow
* [ ] Log
* [ ] Set fastq paired-end reads as secondaryFiles
* [ ] Skip if output exist
* [ ] put file into correct folder, not flatten structure
* [ ] CWL scatter


## Some Tricks

> What is suffix name

Put every steps' short_description into suffix
```
.
├── chrx.fa
├── chrx.1.bt2
├── chrx.2.bt2
├── chrx.3.bt2
├── chrx.4.bt2
├── chrx.rev.1.bt2
├── chrx.rev.2.bt2
├── chrx.synthetic.read1.fq
├── chrx.synthetic.read2.fq
├── chrx.synthetic.bowtie.bam
└── chrx.synthetic.bowtie.bam.bai
```


> Why don't write the default arguments with expression like below?

Because cwl doesn't support this syntax

``` yml
inputs:
  output_bam_file:
    type: string
    inputBinding:
      prefix: -o
    default: $( inputs.unsort_bam.nameroot + ".bam" )
```

Thus, arguments part is used instead.

``` yml
arguments:
  - valueFrom: $( inputs.unsort_bam.nameroot + ".bam" )
    prefix: -o
```

## Result of cwl-airflow

The dependency graph of cwl

![](https://raw.githubusercontent.com/linnil1/cwl_example/main/airflow/cwl-graph.png)

The status of cwl of each triggers

![](https://raw.githubusercontent.com/linnil1/cwl_example/main/airflow/cwl-status.png)


## Reference
* [CWL spec](https://www.commonwl.org/v1.2/Workflow.html)
* [cwltool](https://github.com/common-workflow-language/cwltool)
* [cwl-airflow](https://github.com/Barski-lab/cwl-airflow)
* bowtie2
* samtools
* podman


## LICENSE
MIT
