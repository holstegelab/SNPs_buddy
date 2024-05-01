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
                mutations[sample_id].append(record.ID)

    return mutations
## get actual GT with ID, ref and alt allele for specific site/sample
## filter sample-wise

def print_mutation_table(mutations):
    print("Sample - Mutations")
    for sample, mutations_list in mutations.items():
        if mutations_list:
            print(f"{sample} - {', '.join(mutations_list)}")
        else:
            print(f"{sample} - None")

mutations = find_mutations(vcf_file)
print_mutation_table(mutations)

def print_mutation_table(mutations, output_file):
    with open(output_file, 'w') as f:
        f.write("Sample - Mutations\n")
        for sample, mutations_list in mutations.items():
            if mutations_list:
                f.write(f"{sample} - {', '.join(mutations_list)}\n")
            else:
                f.write(f"{sample} - None\n")

output_file = input("Please enter the path to the output file: ")
print_mutation_table(mutations, output_file)