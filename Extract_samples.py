import vcf

def extract_samples(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    samples = vcf_reader.samples
    return samples

vcf_file = input("Please enter the path to the VCF file: ")

def find_mutations(vcf_file):
    mutations = {}
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    for record in vcf_reader:
        for sample in record.samples:
            sample_id = sample.sample
            if sample_id not in mutations:
                mutations[sample_id] = []
            if sample.gt_type in [1, 2]:  # het or homo mutations
                mutation = {
                    "sample_ID": sample_id,
                    "POS": record.POS,
                    "chr": record.CHROM,
                    "ref": record.REF,
                    "alt": ','.join(map(str, record.ALT)),
                    "GT": sample['GT'],
                    "AD": sample.data.AD,
                    "DP": sample.data.DP,
                    "GQ": sample.data.GQ
                }
                if mutation["DP"] >= 10 and mutation["GQ"] >= 20:
                    mutations[sample_id].append(mutation)
    return mutations

mutations = find_mutations(vcf_file)
# print_mutation_table(mutations)

def print_mutation_table(mutations, output_file):
    with open(output_file, 'w') as f:
        f.write("Sample\tChr\tPos\tRef\tAlt\tGT\tAD\n")
        for sample, mutations_list in mutations.items():
            if mutations_list:
                for mutation in mutations_list:
                    f.write(f"{sample}\t{mutation['chr']}\t{mutation['POS']}\t{mutation['ref']}\t{mutation['alt']}\t{mutation['GT']}\t{mutation['AD']}\n")
            else:
                f.write(f"{sample}\tNone\tNone\tNone\tNone\tNone\tNone\n")

output_file = input("Please enter the path to the output TSV file: ")
print_mutation_table(mutations, output_file)
