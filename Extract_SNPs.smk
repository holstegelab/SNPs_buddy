gene = config.get('gene')
genotype_mode = config.get("genotype_mode", "WES")
gvcf_caller = config.get("caller", "BOTH")


def get_regions(lrange):
    """Converts a region describer (tuple format) to a list of regions (string format).
    E.g. [('A', 1, 3, 1), ('A', 2, 4, 2)] -> ['A3H', 'A04']
    """
    res = []

    for component, level, splitnr, ploidy in lrange:
        if ploidy == 1:
            ploidy = 'H'
        else:
            ploidy = ''
        if level == 0:
            region = f'{component}{ploidy}'
        else:
            region = f'{component}{splitnr:0{level}d}{ploidy}'
        res.append(region)

    return res

level2_range_diploid_only = [('A', 2,x,2) for x in range(0,100)] + \
               [('X', 1,x,2) for x in range(0,5)] + \
               [('Y', 1, x,2) for x in range(0,2)]

level2_regions_diploid = get_regions(level2_range_diploid_only)
parts = level2_regions_diploid

rule all:
    input: expand(pj('{gene}/FILTRED_{genotype_mode}_{glnexus_dir}_annotated.hg38_multianno.vcf'), gene=gene, genotype_mode=genotype_mode, glnexus_dir=gvcf_caller)


rule extract_per_part:
    input: pj('{genotype_mode}_{glnexus_dir}custom/ANNOTATED/{region}_annotated.hg38_multianno.vcf.gz')
    output: temp(pj('{gene}/{region}_annotated.hg38_multianno.vcf.gz'))
    conda: "envs/snp_buddies.yaml"
    shell: """
            bcftool view -Oz -o {output}  --write-index --exclude-uncalled --threads 2 --include 'INFO/Gene.ensGene=="{gene}"' {input}
            """

rule gather_parts:
    input: expand(pj('{gene}/{genotype_mode}_{glnexus_dir}_{region}_annotated.hg38_multianno.vcf.gz'), region=parts, allow_missing=True)
    output: pj('{gene}/{genotype_mode}_{glnexus_dir}_annotated.hg38_multianno.vcf')
    conda: "envs/snp_buddies.yaml"
    shell: """
            bcftools concat --threads 2 -Ov -o {output} {input}
            """

rule quality_check:
    input: pj('{gene}/{genotype_mode}_{glnexus_dir}_annotated.hg38_multianno.vcf')
    output: pj('{gene}/FILTRED_{genotype_mode}_{glnexus_dir}_annotated.hg38_multianno.vcf')
    conda: "envs/snp_buddies.yaml"
    shell: """
            bcftools view --write-index -Ov -o {output} --threads 2 --include 'QUAL>20 & FORMAT/DP>10' {input} 
            
            bcftools stats {output} > {output}.stats
            """


