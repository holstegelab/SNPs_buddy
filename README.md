# SNPs_buddy

Utilities for extracting and summarizing variants from annotated VCF files.

All commands below assume you are running from the repository root:

```bash
cd /path/to/SNPs_buddy
```

## Requirements

- `conda` or `miniconda`
- input VCF files for the Snakemake workflows must be split into per-part files named `{part}.annotated.vcf.gz`
- by default the Snakemake workflows read those files from `/project/holstegelab/Share/NL_VUMC_joint_calling_splitted/ANNOTATED`

Create the environment once:

```bash
conda env create -f envs/snp_buddies.yaml
conda activate snp_buddies_env
```

The environment defined in [envs/snp_buddies.yaml](envs/snp_buddies.yaml) includes `bcftools`, `htslib`, `gatk4`, `pyvcf`, `python 3.11`, and `snakemake`.

## Extract variants by gene or BED

Use [Extract_SNPs.smk](Extract_SNPs.smk).

Required config:

- exactly one of `gene=GENE_NAME` or `bed=/path/to/regions.bed`

Optional config:

- `input_dir=/path/to/annotated_parts`

Behavior:

- `gene` mode keeps variants where `INFO/Gene.ensGene` or `INFO/Gene.refGene` matches the provided gene name
- `bed` mode keeps variants overlapping the provided BED intervals via `bcftools view -R`
- the output directory name is the gene name or the BED basename with `.bed` or `.bed.gz` removed

Local example:

```bash
snakemake -c 2 \
  --rerun-incomplete \
  --conda-frontend conda \
  --use-conda \
  --snakefile Extract_SNPs.smk \
  --config gene=APP
```

BED example:

```bash
snakemake -c 2 \
  --rerun-incomplete \
  --conda-frontend conda \
  --use-conda \
  --snakefile Extract_SNPs.smk \
  --config bed=/path/to/regions.bed
```

Cluster example:

```bash
snakemake --jobs 20 \
  --cluster "sbatch -n {resources.n} --mem {resources.mem_mb} -p {resources.partition} -t {resources.time_min}" \
  --rerun-incomplete \
  --conda-frontend conda \
  --keep-going \
  --use-conda \
  --snakefile Extract_SNPs.smk \
  --config gene=APP
```

Outputs for target `APP`:

- `APP/APP.vcf`
- `APP/FILTRED_APP.vcf`
- `APP/FILTRED_APP.vcf.stats`
- `APP/FILTRED_APP.missigness.tsv`
- `APP/benchs/`

## Extract variants by sample name

Use [Sample_by_samplename.smk](Sample_by_samplename.smk).

Required config:

- `sample_names=Sample1,Sample2,...`

Optional config:

- `dir_name=output_directory_name`
- `input_dir=/path/to/annotated_parts`

Behavior:

- `sample_names` is passed directly to `bcftools view -s`, so use a comma-separated list with no spaces
- if `dir_name` is not provided, the workflow uses the first sample name before the first comma

Example:

```bash
snakemake -c 2 \
  --rerun-incomplete \
  --conda-frontend conda \
  --use-conda \
  --snakefile Sample_by_samplename.smk \
  --config sample_names=Sample1,Sample2 dir_name=my_samples
```

Cluster example:

```bash
snakemake --jobs 20 \
  --cluster "sbatch -n {resources.n} --mem {resources.mem_mb} -p {resources.partition} -t {resources.time_min}" \
  --rerun-incomplete \
  --conda-frontend conda \
  --keep-going \
  --use-conda \
  --snakefile Sample_by_samplename.smk \
  --config sample_names=Sample1,Sample2 dir_name=my_samples
```

Outputs are written under `dir_name/`, with benchmarks under `dir_name/benchs/`.

## Summarize variants per sample from a VCF

Use [Extract_samples.py](Extract_samples.py).

Run:

```bash
python Extract_samples.py
```

The script interactively asks for:

- the input VCF path
- the first output TSV path
- the second output TSV path

Behavior:

- reads all sample names from the VCF
- keeps sample variants with genotype type `1` or `2`
- requires `DP >= 10` and `GQ >= 20`
- assigns an `Importance` label based on annotation fields and score thresholds in the script

Outputs:

- the first TSV contains `Sample`, `Chr`, `Pos`, `Ref`, `Alt`, `GT`, `AD`
- the second TSV contains `Importance`, the same core columns, and one column for each INFO field observed in the VCF

## Count low-depth calls per site

Use [low_dp_counts.smk](low_dp_counts.smk).

Required config:

- `vcf=/path/to/file.vcf`

Optional config:

- `region=/path/to/regions.bed`

Behavior:

- without `region`, the workflow scans the whole VCF
- with `region`, the workflow bgzips and indexes the VCF, filters to the BED intervals, then computes the same metrics on the filtered records
- for each site, it counts samples with non-missing `DP < 10`

Example:

```bash
snakemake -c 1 \
  --rerun-incomplete \
  --conda-frontend conda \
  --use-conda \
  --snakefile low_dp_counts.smk \
  --config vcf=/path/to/file.vcf
```

Example with BED regions:

```bash
snakemake -c 1 \
  --rerun-incomplete \
  --conda-frontend conda \
  --use-conda \
  --snakefile low_dp_counts.smk \
  --config vcf=/path/to/file.vcf region=/path/to/regions.bed
```

Output naming:

- without `region`: `<vcf_basename>.low_dp_counts.tsv`
- with `region`: `<region_basename>.low_dp_counts.tsv`

Output columns:

- chromosome
- position
- number of samples with `DP < 10`
- percentage of samples with `DP < 10`

## Compare cohort allele frequencies and p-values

Use [cohort_af_pvalues.py](cohort_af_pvalues.py).

Required arguments:

- `--vcf /path/to/file.vcf` or `--vcf /path/to/file.vcf.gz`
- one or more `--cohort NAME=/path/to/cohort.tsv`
- `--output /path/to/results.tsv`

Optional arguments:

- `--sample-column N`
- `--has-header`
- `--missing-value VALUE`
- `--min-dp N`
- `--manhattan-prefix /path/to/plot_prefix`
- `--genome-wide-threshold FLOAT`

Input expectations:

- cohort files must be tab-delimited
- sample IDs are read from the selected column
- `--sample-column` is 1-based and defaults to `2`
- duplicate sample IDs inside a cohort file are ignored with a warning

Example:

```bash
python cohort_af_pvalues.py \
  --vcf /path/to/file.vcf.gz \
  --cohort Cases=/path/to/cases.tsv \
  --cohort Controls=/path/to/controls.tsv \
  --output /path/to/results.tsv
```

Behavior:

- writes one row per SNP ALT allele
- skips indels
- treats sample genotypes with `DP < 30` as missing by default
- if `--manhattan-prefix` is provided, writes one PNG Manhattan plot per pairwise p-value column

Output columns:

- `chrom`
- `pos`
- `snp_id`
- `ref`
- `alt`
- `AF_COHORTNAME` for each cohort
- `missingness_pct_COHORTNAME` for each cohort
- `oddsratio_COHORT1_vs_COHORT2` for each cohort pair
- `pvalue_COHORT1_vs_COHORT2` for each cohort pair
