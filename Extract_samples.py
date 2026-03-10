import vcf

HIGH_IMPORTANCE_DP = 30
HIGH_CADD_THRESHOLD = 20.0
HIGH_REVEL_THRESHOLD = 0.5

def extract_samples(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    samples = vcf_reader.samples
    return samples

vcf_file = input("Please enter the path to the VCF file: ")

def serialize_info_value(value):
    if isinstance(value, (list, tuple)):
        return ",".join(map(str, value))
    if value is None:
        return "."
    return str(value)

def get_max_numeric_info_value(info_dict, keys):
    best = None
    for key in keys:
        raw_value = info_dict.get(key)
        if raw_value in (None, "."):
            continue
        for token in str(raw_value).replace("|", ",").split(","):
            token = token.strip()
            if not token or token == ".":
                continue
            try:
                value = float(token)
            except ValueError:
                continue
            if best is None or value > best:
                best = value
    return best

def contains_keyword_in_info(info_dict, keys, keywords):
    text = " ".join(str(info_dict.get(key, "")) for key in keys).lower()
    return any(keyword in text for keyword in keywords)

def classify_importance(info_dict, dp):
    dp = dp if dp is not None else 0

    high_impact_keys = [
        "Func.refGene", "Func.ensGene",
        "ExonicFunc.refGene", "ExonicFunc.ensGene",
        "LoF", "LoF_filter", "LoF_flags"
    ]
    high_impact_keywords = ["lof", "frameshift", "splice", "stopgain", "stoploss"]
    is_high_impact = contains_keyword_in_info(info_dict, high_impact_keys, high_impact_keywords)

    if dp >= HIGH_IMPORTANCE_DP and is_high_impact:
        return "High Importnat"

    nonsynonymous_keys = ["ExonicFunc.refGene", "ExonicFunc.ensGene"]
    nonsynonymous_keywords = ["nonsynonymous", "missense"]
    is_nonsynonymous = contains_keyword_in_info(info_dict, nonsynonymous_keys, nonsynonymous_keywords)

    cadd_score = get_max_numeric_info_value(info_dict, ["CADD_phred", "CADD_PHRED", "CADD"])
    revel_score = get_max_numeric_info_value(info_dict, ["REVEL_score", "REVEL", "dbNSFP_REVEL_score"])

    has_high_cadd = cadd_score is not None and cadd_score >= HIGH_CADD_THRESHOLD
    has_high_revel = revel_score is not None and revel_score >= HIGH_REVEL_THRESHOLD

    if is_nonsynonymous and (has_high_cadd or has_high_revel):
        return "Important"

    return "."

def find_mutations(vcf_file):
    mutations = {}
    info_keys = []
    seen_info_keys = set()
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
                    "GQ": sample.data.GQ,
                    "INFO": {k: serialize_info_value(v) for k, v in record.INFO.items()}
                }
                dp = mutation["DP"] if mutation["DP"] is not None else 0
                gq = mutation["GQ"] if mutation["GQ"] is not None else 0
                if dp >= 10 and gq >= 20:
                    mutation["Importance"] = classify_importance(mutation["INFO"], dp)
                    for info_key in mutation["INFO"]:
                        if info_key not in seen_info_keys:
                            seen_info_keys.add(info_key)
                            info_keys.append(info_key)
                    mutations[sample_id].append(mutation)
    return mutations, info_keys

mutations, info_keys = find_mutations(vcf_file)
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

def print_mutation_table_with_info(mutations, info_keys, output_file):
    def importance_rank(value):
        ranks = {"High Importnat": 0, "Important": 1}
        return ranks.get(value, 2)

    def pos_sort_value(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return float("inf")

    with open(output_file, 'w') as f:
        header = ["Importance", "Sample", "Chr", "Pos", "Ref", "Alt", "GT", "AD"] + info_keys
        f.write("\t".join(header) + "\n")

        rows = []
        for sample, mutations_list in mutations.items():
            if mutations_list:
                for mutation in mutations_list:
                    row = [
                        mutation.get("Importance", "."),
                        sample,
                        str(mutation["chr"]),
                        str(mutation["POS"]),
                        str(mutation["ref"]),
                        str(mutation["alt"]),
                        str(mutation["GT"]),
                        str(mutation["AD"])
                    ]
                    row.extend(mutation["INFO"].get(key, ".") for key in info_keys)
                    rows.append(row)
            else:
                row = [".", sample, "None", "None", "None", "None", "None", "None"]
                row.extend("." for _ in info_keys)
                rows.append(row)

        rows.sort(key=lambda row: (importance_rank(row[0]), row[1], row[2], pos_sort_value(row[3])))
        for row in rows:
            f.write("\t".join(row) + "\n")

output_file = input("Please enter the path to the output TSV file: ")
print_mutation_table(mutations, output_file)

output_file_with_info = input("Please enter the path to the 2nd output TSV file (with INFO columns): ")
print_mutation_table_with_info(mutations, info_keys, output_file_with_info)
