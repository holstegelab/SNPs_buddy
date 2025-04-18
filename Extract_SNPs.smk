import os
from constants import *

gene = config.get('gene')
genotype_mode = config.get("genotype_mode", "WES")
gvcf_caller = config.get("caller", "BOTH")
input_dir = config.get('input_dir', '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')

rule all:
    input: expand(pj('{gene}/FILTRED_{gene}.vcf'), gene=gene),
              expand(pj('{gene}/FILTRED_{gene}.missigness.tsv'), gene=gene),
              expand(pj('{gene}/FILTRED_{gene}.vcf.stats'), gene=gene)


rule extract_per_part:
    input: pj(input_dir, '{parts}.annotated.vcf.gz')
    output: temp(ensure(pj('{gene}/{parts}_annotated.vcf.gz'), non_empty=True))
    conda: "envs/snp_buddies.yaml"
    benchmark: pj('{gene}/benchs/{parts}_extract.benchmark')
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:45:00'
    shell: """
            bcftools view -Oz -o {output}  --exclude-uncalled --threads 2 --include 'INFO/Gene.ensGene=="{gene}" || INFO/Gene.refGene=="{gene}"'  {input}
            """

rule gather_parts:
    input: expand(pj('{gene}/{parts}_annotated.vcf.gz'), parts=parts, allow_missing=True)
    output: pj('{gene}/{gene}.vcf')
    conda: "envs/snp_buddies.yaml"
    benchmark: pj('{gene}/benchs/{gene}.gather.benchmark')
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:20:00'
    shell: """
            bcftools concat --threads 2 -Ov -o {output} {input}
            """

rule quality_check:
    input: pj('{gene}/{gene}.vcf')
    output: pj('{gene}/FILTRED_{gene}.vcf')
    conda: "envs/snp_buddies.yaml"
    benchmark: pj('{gene}/benchs/{gene}.quality.benchmark')
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:20:00'
    shell: """
            bcftools view -Ov -o {output} --threads 2 --include 'QUAL>20 & FORMAT/DP>10' {input} 
            
            bcftools stats {output} > {output}.stats
            """

rule calculate_missigness:
    input: pj('{gene}/FILTRED_{gene}.vcf')
    output: pj('{gene}/FILTRED_{gene}.missigness.tsv')
    conda: "envs/snp_buddies.yaml"
    resources: n = 1,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:20:00'
    shell:
        """
        bcftools query -f '%CHROM\t%POS[\t%DP]\n' {input} | awk 'BEGIN {{ OFS="\t" }}
          NR==1 {{total_samples = NF - 2}}
                 {{count=0
                     for(i=3; i<=NF; i++) {{
                         if($i != "." && $i < 10) count++}}
                     percent = (count / total_samples) * 100
                     printf "%s\t%s\t%d\t%.2f%%\n", $1, $2, count, percent
                 }}' > {output}
        """

