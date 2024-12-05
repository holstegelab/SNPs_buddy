import os
from constants import *

sample_names = config.get('sample_names')

dir_name = config.get('dir_name', sample_names.split(',')[0])
input_dir = config.get('input_dir', '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')

rule all:
    input: expand(pj('{dir_name}/{sample_name}.vcf'), sample_name=sample_names, dir_name=dir_name)

rule extract_sample_by_sample_name:
    input: pj(input_dir, '{parts}.annotated.vcf.gz')
    output: temp(pj('{dir_name}/{parts}_annotated.vcf.gz'))
    conda: "envs/snp_buddies.yaml"
    params: sn = sample_names
    threads: 2
    resources: n = 2,
                mem_mb = 8096,
                partition = 'short',
                time_min = 120

    shell:
        """
        bcftools view --threads 2 -Oz -o {output}  -s "{params.sn}" {input}
        """

rule gather_vcfs:
    input: expand(pj('{dir_name}/{parts}_annotated.vcf.gz'), parts=parts, allow_missing=True)
    output: pj('{dir_name}/{sample_name}.vcf')
    conda: "envs/snp_buddies.yaml"
    threads: 2
    resources: n = 2,
                mem_mb = 8096,
                partition = 'short',
                time_min = 120
    shell:
        """
        bcftools concat --threads 2 -Ov -o {output} {input}
        """
