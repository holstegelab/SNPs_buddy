import os
from constants import *

gene = config.get('gene')
bed = config.get('bed')
genotype_mode = config.get("genotype_mode", "WES")
gvcf_caller = config.get("caller", "BOTH")
input_dir = config.get('input_dir', '/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED')


def get_target_name(gene_name, bed_path):
    if gene_name and bed_path:
        raise ValueError("Provide either 'gene' or 'bed', not both.")
    if not gene_name and not bed_path:
        raise ValueError("Provide one of 'gene' or 'bed'.")

    if bed_path:
        bed_name = os.path.basename(str(bed_path))
        if bed_name.endswith('.bed.gz'):
            return bed_name[:-7]
        if bed_name.endswith('.bed'):
            return bed_name[:-4]
        return os.path.splitext(bed_name)[0]

    return str(gene_name)


target_name = get_target_name(gene, bed)

if bed:
    bed = os.path.abspath(str(bed))
    if not os.path.exists(bed):
        raise FileNotFoundError(f"BED file not found: {bed}")

rule all:
    input: pj(target_name, f'FILTRED_{target_name}.vcf'),
           pj(target_name, f'FILTRED_{target_name}.missigness.tsv'),
           pj(target_name, f'FILTRED_{target_name}.vcf.stats')


rule extract_per_part:
    input: pj(input_dir, '{parts}.annotated.vcf.gz')
    output: temp(pj(target_name, '{parts}_annotated.vcf.gz'))
    conda: "envs/snp_buddies.yaml"
    benchmark: pj(target_name, 'benchs/{parts}_extract.benchmark')
    params: gene = gene,
            bed = bed
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:45:00'
    shell: """
            if [ "{params.bed}" = "None" ]; then
                bcftools view -Oz -o {output} --exclude-uncalled --threads 2 --include 'INFO/Gene.ensGene=="{params.gene}" || INFO/Gene.refGene=="{params.gene}"' {input}
            else
                bcftools view -Oz -o {output} --exclude-uncalled --threads 2 -R "{params.bed}" {input}
            fi
            """

rule gather_parts:
    input: expand(pj(target_name, '{parts}_annotated.vcf.gz'), parts=parts, allow_missing=True)
    output: pj(target_name, f'{target_name}.vcf')
    conda: "envs/snp_buddies.yaml"
    benchmark: pj(target_name, f'benchs/{target_name}.gather.benchmark')
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:20:00'
    shell: """
            bcftools concat --threads 2 -Ov -o {output} {input}
            """

rule quality_check:
    input: pj(target_name, f'{target_name}.vcf')
    output: vcf = pj(target_name, f'FILTRED_{target_name}.vcf'),
            stats = pj(target_name, f'FILTRED_{target_name}.vcf.stats')
    conda: "envs/snp_buddies.yaml"
    benchmark: pj(target_name, f'benchs/{target_name}.quality.benchmark')
    resources: n = 2,
                mem_mb = 1000,
                partition = 'normal',
                time_min = '00:20:00'
    shell: """
            bcftools view -Ov -o {output.vcf} --threads 2 --include 'QUAL>20 & FORMAT/DP>10' {input} 
            
            bcftools stats {output.vcf} > {output.vcf}.stats
            """

rule calculate_missigness:
    input: pj(target_name, f'FILTRED_{target_name}.vcf')
    output: pj(target_name, f'FILTRED_{target_name}.missigness.tsv')
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
