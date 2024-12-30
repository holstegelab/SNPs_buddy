import os
from constants import *

gene = config.get('gene')
genotype_mode = config.get("genotype_mode", "WES")
gvcf_caller = config.get("caller", "BOTH")
input_dir = config.get('input_dir', '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')

rule all:
    input: expand(pj('{gene}/FILTRED_{gene}.vcf'), gene=gene)


rule extract_per_part:
    input: pj(input_dir, '{parts}.annotated.vcf.gz')
    output: temp(pj('{parts}/{region}_annotated.vcf.gz'))
    conda: "envs/snp_buddies.yaml"
    resources: n = 2,
                mem_mb = 8000,
                partition = 'normal',
                time_min = 720
    shell: """
            bcftools view -Oz -o {output}  --exclude-uncalled --threads 2 --include 'INFO/Gene.ensGene=="{gene}" || INFO/Gene.refGene=="{gene}"'  {input}
            """

rule gather_parts:
    input: expand(pj('{gene}/{region}_annotated.vcf.gz'), region=parts, allow_missing=True)
    output: pj('{gene}/{gene}.vcf')
    conda: "envs/snp_buddies.yaml"
    resources: n = 2,
                mem_mb = 8000,
                partition = 'normal',
                time_min = 720
    shell: """
            bcftools concat --threads 2 -Ov -o {output} {input}
            """

rule quality_check:
    input: pj('{gene}/{gene}.vcf')
    output: pj('{gene}/FILTRED_{gene}.vcf')
    conda: "envs/snp_buddies.yaml"
    resources: n = 2,
                mem_mb = 8000,
                partition = 'normal',
                time_min = 720
    shell: """
            bcftools view -Ov -o {output} --threads 2 --include 'QUAL>20 & FORMAT/DP>10' {input} 
            
            bcftools stats {output} > {output}.stats
            """


