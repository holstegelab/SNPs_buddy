import os

vcf_file = config.get('vcf')
region = config.get('region')
if region == "None":
    region = vcf_file.split('/')[-1].split('.')[0]
else:
    region = region.split('/')[-1].split('.')[0]


rule get_precentage:
    input: vcf = vcf_file,
            bed = region
    output: "{region_name}.low_dp_counts.tsv"
    conda: "envs/snp_buddies.yaml"
    params: bed = lambda wildcards, input: int(input.bed) if input.bed != "None" else 0
    default_target: True
    shell:
        """
        if [{params.bed} == 0]
        then
            bcftools query -f '%CHROM\t%POS[\t%DP]\n' FILTRED_SORL1_CDSs.vcf | awk 'BEGIN { OFS="\t" }
             NR==1 {
                 total_samples = NF - 2
             }
             {
                 count=0
                 for(i=3; i<=NF; i++) {
                     if($i != "." && $i < 10) count++
                 }
                 percent = (count / total_samples) * 100
                 printf "%s\t%s\t%d\t%.2f%%\n", $1, $2, count, percent
             }' > {output}

        else
            bcftools view {input.vcf} -R {input.bed} > FILTRED_SORL1_CDSs.vcf
    
            bcftools query -f '%CHROM\t%POS[\t%DP]\n' FILTRED_SORL1_CDSs.vcf | awk 'BEGIN { OFS="\t" }
                     NR==1 {
                         total_samples = NF - 2
                     }
                     {
                         count=0
                         for(i=3; i<=NF; i++) {
                             if($i != "." && $i < 10) count++
                         }
                         percent = (count / total_samples) * 100
                         printf "%s\t%s\t%d\t%.2f%%\n", $1, $2, count, percent
                     }' > {output}
            fi
            """