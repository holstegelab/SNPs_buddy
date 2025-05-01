import os

vcf_file = config.get('vcf')
region = config.get('region')
if region is None:
    vcf_file = str(vcf_file)
    region_name = os.path.basename(vcf_file).split('.')[0]
else:
    region = str(region)
    region_name = os.path.basename(region).split('.')[0]


rule all:
    input: expand("{region_name}.low_dp_counts.tsv", region_name = region_name)


rule get_percentage:
    input:
        vcf = vcf_file
    output:
        "{region_name}.low_dp_counts.tsv"
    conda:
        "envs/snp_buddies.yaml"
    params:
        bed = region if region else None
    shell:
        """
        if [ "{params.bed}" = "None" ]; then
            bcftools query -f '%CHROM\t%POS[\t%DP]\n' {input.vcf} | \
            awk 'BEGIN { OFS="\\t" }
                {
                    count = 0
                    if (NR == 1) total_samples = NF - 2
                    for (i = 3; i <= NF; i++) {
                        if ($i != "." && $i < 10) count++
                    }
                    percent = (count / total_samples) * 100
                    printf "%s\\t%s\\t%d\\t%.2f%%\\n", $1, $2, count, percent
                }' > {output}
        else
            bcftools view {input.vcf} -R {params.bed} > tmp_filtered.vcf
            bcftools query -f '%CHROM\t%POS[\t%DP]\n' tmp_filtered.vcf | \
            awk 'BEGIN { OFS="\\t" }
                {
                    count = 0
                    if (NR == 1) total_samples = NF - 2
                    for (i = 3; i <= NF; i++) {
                        if ($i != "." && $i < 10) count++
                    }
                    percent = (count / total_samples) * 100
                    printf "%s\\t%s\\t%d\\t%.2f%%\\n", $1, $2, count, percent
                }' > {output}
        fi
        """
