from cwl_airflow.extensions.cwldag import CWLDAG
dag = CWLDAG(
    workflow="$PWD/pipeline.cwl",
    dag_id="bowtie2_cwl",
    default_args={"user_space_docker_cmd": "podman"}
)
