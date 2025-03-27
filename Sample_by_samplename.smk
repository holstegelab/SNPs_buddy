import os
from constants import *

sample_names = config.get('sample_names')

dir_name = config.get('dir_name', sample_names.split(',')[0])
input_dir = config.get('input_dir', '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')

rule all:
    input: expand(pj('{dir_name}/{sample_name}.vcf'), sample_name=sample_names, dir_name=dir_name)

rule extract_sample_by_sample_name:
    input: pj(input_dir, '{parts}.annotated.vcf.gz')
    output: temp(ensure(pj('{dir_name}/{parts}_annotated.vcf.gz'), non_empty=True))
    conda: "envs/snp_buddies.yaml"
    params: sn = sample_names
    benchmark: pj('{dir_name}/benchs/{parts}_annotated.benchmark')
    threads: 4
    resources: n = 4,
                mem_mb = 4000,
                partition = 'normal',
                time_min = '00:59:00'

    shell:
        """
        bcftools view --threads {resources.n} -Oz -o {output}  -s "{params.sn}" {input}
        """

rule gather_vcfs:
    input: expand(pj('{dir_name}/{parts}_annotated.vcf.gz'), parts=parts, allow_missing=True)
    output: pj('{dir_name}/{sample_name}.vcf')
    conda: "envs/snp_buddies.yaml"
    threads: 2
    benchmark: pj('{dir_name}/benchs/{sample_name}.benchmark')
    resources: n = 2,
                mem_mb = 2000,
                partition = 'normal',
                time_min = '00:45:00'
    shell:
        """
        bcftools concat --threads {resources.n} -Ov -o {output} {input}
        """
