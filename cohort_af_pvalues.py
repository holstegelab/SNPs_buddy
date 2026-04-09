#!/usr/bin/env python3
"""Compute cohort AF, missingness, odds ratios, and pairwise p-values from a filtered VCF."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import math
import re
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read a filtered VCF and cohort manifest TSV files, then write one TSV "
            "row per SNP ALT allele with cohort AFs, missingness, pairwise odds "
            "ratios, and Fisher exact test p-values."
        )
    )
    parser.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file.")
    parser.add_argument(
        "--cohort",
        action="append",
        required=True,
        metavar="NAME=PATH",
        help=(
            "Cohort manifest in NAME=PATH form. Repeat for each cohort. "
            "Sample IDs are read from the selected column."
        ),
    )
    parser.add_argument("--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "--sample-column",
        type=int,
        default=2,
        help="1-based column index containing the sample ID in cohort TSVs. Default: 2.",
    )
    parser.add_argument(
        "--has-header",
        action="store_true",
        help="Skip the first row in each cohort TSV.",
    )
    parser.add_argument(
        "--missing-value",
        default="NA",
        help="String used when AF or p-value cannot be computed. Default: NA.",
    )
    parser.add_argument(
        "--min-dp",
        type=int,
        default=30,
        help=(
            "Minimum per-sample DP required to count a genotype. Lower-depth "
            "genotypes are treated as missing. Default: 30."
        ),
    )
    parser.add_argument(
        "--manhattan-prefix",
        help=(
            "If set, generate one Manhattan plot PNG per pairwise p-value column "
            "using this path prefix."
        ),
    )
    parser.add_argument(
        "--genome-wide-threshold",
        type=float,
        default=5e-8,
        help="Significance threshold used in Manhattan plots. Default: 5e-8.",
    )
    return parser.parse_args()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def warn(message: str) -> None:
    print(f"Warning: {message}", file=sys.stderr)


def parse_cohort_argument(argument: str) -> tuple[str, Path]:
    if "=" not in argument:
        raise ValueError(f"Invalid --cohort value {argument!r}. Expected NAME=PATH.")
    name, raw_path = argument.split("=", 1)
    name = name.strip()
    if not name:
        raise ValueError(f"Invalid --cohort value {argument!r}. Cohort name is empty.")
    path = Path(raw_path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Cohort file not found: {path}")
    return name, path


def load_cohort_samples(path: Path, sample_column: int, has_header: bool) -> list[str]:
    if sample_column < 1:
        raise ValueError("--sample-column must be at least 1.")

    column_index = sample_column - 1
    samples: list[str] = []
    seen: set[str] = set()

    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        if has_header:
            next(reader, None)

        for line_number, row in enumerate(reader, start=2 if has_header else 1):
            if not row or all(not cell.strip() for cell in row):
                continue
            if len(row) <= column_index:
                raise ValueError(
                    f"{path}: line {line_number} has {len(row)} columns, "
                    f"but sample column {sample_column} was requested."
                )

            sample_id = row[column_index].strip()
            if not sample_id:
                continue
            if sample_id in seen:
                warn(f"{path}: duplicate sample {sample_id!r} ignored.")
                continue
            seen.add(sample_id)
            samples.append(sample_id)

    if not samples:
        raise ValueError(f"No samples found in cohort file: {path}")

    return samples


def load_vcf_samples(vcf_path: Path) -> list[str]:
    with open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"VCF header line not found in {vcf_path}")


def parse_called_alleles(genotype: str) -> list[int]:
    if not genotype or genotype in {".", "./.", ".|."}:
        return []

    genotype = genotype.replace("|", "/")
    alleles = genotype.split("/")
    if any(allele == "." for allele in alleles):
        return []

    try:
        return [int(allele) for allele in alleles]
    except ValueError:
        return []


def extract_format_value(sample_field: str, field_index: int | None) -> str:
    if field_index is None:
        return ""
    parts = sample_field.split(":")
    if field_index >= len(parts):
        return ""
    return parts[field_index]


def parse_dp_value(raw_dp: str) -> int | None:
    if not raw_dp or raw_dp == ".":
        return None
    try:
        return int(raw_dp)
    except ValueError:
        return None


def allele_frequency(alt_count: int, called_alleles: int) -> float | None:
    if called_alleles == 0:
        return None
    return alt_count / called_alleles


def missingness_percentage(missing_samples: int, total_samples: int) -> float | None:
    if total_samples == 0:
        return None
    return 100.0 * missing_samples / total_samples


def log_choose(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeometric_probability(x: int, row1: int, row2: int, col1: int, total: int) -> float:
    return math.exp(log_choose(row1, x) + log_choose(row2, col1 - x) - log_choose(total, col1))


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float | None:
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2

    if row1 == 0 or row2 == 0 or total == 0:
        return None

    observed = hypergeometric_probability(a, row1, row2, col1, total)
    min_x = max(0, col1 - row2)
    max_x = min(row1, col1)
    tolerance = observed * 1e-12 + 1e-15

    p_value = 0.0
    for x in range(min_x, max_x + 1):
        probability = hypergeometric_probability(x, row1, row2, col1, total)
        if probability <= observed + tolerance:
            p_value += probability

    return min(p_value, 1.0)


def odds_ratio(a: int, b: int, c: int, d: int) -> float | None:
    if (a + b) == 0 or (c + d) == 0:
        return None

    numerator = a * d
    denominator = b * c

    if denominator == 0:
        if numerator == 0:
            return None
        return math.inf

    return numerator / denominator


def format_number(value: float | None, missing_value: str) -> str:
    if value is None:
        return missing_value
    return f"{value:.10g}"


def parse_numeric(value: str, missing_value: str) -> float | None:
    if value == missing_value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def build_cohort_indices(
    vcf_samples: list[str],
    cohort_to_samples: list[tuple[str, list[str]]],
) -> list[tuple[str, list[int]]]:
    sample_to_index = {sample: index for index, sample in enumerate(vcf_samples)}
    cohort_indices: list[tuple[str, list[int]]] = []

    for cohort_name, cohort_samples in cohort_to_samples:
        present = [sample_to_index[sample] for sample in cohort_samples if sample in sample_to_index]
        missing = [sample for sample in cohort_samples if sample not in sample_to_index]

        if missing:
            warn(
                f"{cohort_name}: {len(missing)} samples are not present in the VCF "
                f"and will be ignored."
            )
        if not present:
            raise ValueError(f"{cohort_name}: no cohort samples were found in the VCF.")

        cohort_indices.append((cohort_name, present))

    return cohort_indices


def chromosome_sort_key(chrom: str) -> tuple[int, int | str]:
    normalized = chrom.removeprefix("chr").removeprefix("CHR")
    special_order = {"X": 23, "Y": 24, "M": 25, "MT": 25}

    if normalized.isdigit():
        return (0, int(normalized))
    if normalized in special_order:
        return (0, special_order[normalized])
    return (1, normalized)


def sanitize_filename_part(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return cleaned.strip("_") or "plot"


def make_plot_path(prefix: str, pvalue_column: str) -> Path:
    prefix_path = Path(prefix).expanduser()
    filename = f"{prefix_path.name}_{sanitize_filename_part(pvalue_column)}.png"
    if prefix_path.parent == Path("."):
        return Path(filename)
    return prefix_path.parent / filename


def draw_manhattan_plots(
    plot_data: dict[str, list[tuple[str, int, str, float]]],
    output_prefix: str,
    genome_wide_threshold: float,
) -> None:
    if genome_wide_threshold <= 0 or genome_wide_threshold >= 1:
        raise ValueError("--genome-wide-threshold must be between 0 and 1.")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Matplotlib is required for Manhattan plots. Install matplotlib or "
            "run without --manhattan-prefix."
        ) from exc

    threshold_y = -math.log10(genome_wide_threshold)
    colors = ("#1f77b4", "#ff7f0e")

    for pvalue_column, points in plot_data.items():
        valid_points = [point for point in points if point[3] > 0]
        if not valid_points:
            warn(f"{pvalue_column}: no valid p-values available for plotting.")
            continue

        valid_points.sort(key=lambda point: (chromosome_sort_key(point[0]), point[1], point[2]))

        chromosome_order: list[str] = []
        chromosome_max: dict[str, int] = {}
        for chrom, pos, _snp_id, _pvalue in valid_points:
            if chrom not in chromosome_max:
                chromosome_order.append(chrom)
                chromosome_max[chrom] = pos
            else:
                chromosome_max[chrom] = max(chromosome_max[chrom], pos)

        offsets: dict[str, int] = {}
        tick_positions: list[float] = []
        running_offset = 0
        for chrom in chromosome_order:
            offsets[chrom] = running_offset
            tick_positions.append(running_offset + chromosome_max[chrom] / 2)
            running_offset += chromosome_max[chrom]

        fig, ax = plt.subplots(figsize=(14, 6))
        significant_points: list[tuple[float, float, str]] = []

        for chrom_index, chrom in enumerate(chromosome_order):
            chrom_points = [point for point in valid_points if point[0] == chrom]
            xs = [offsets[chrom] + point[1] for point in chrom_points]
            ys = [-math.log10(point[3]) for point in chrom_points]
            ax.scatter(xs, ys, s=12, color=colors[chrom_index % len(colors)], label=chrom)

            for x, y, (_chrom, _pos, snp_id, pvalue) in zip(xs, ys, chrom_points):
                if pvalue < genome_wide_threshold:
                    significant_points.append((x, y, snp_id))

        for x, y, snp_id in significant_points:
            ax.annotate(
                snp_id,
                xy=(x, y),
                xytext=(0, 4),
                textcoords="offset points",
                fontsize=7,
                rotation=30,
                ha="left",
                va="bottom",
            )

        ax.axhline(
            threshold_y,
            color="red",
            linestyle="--",
            linewidth=1,
            label=f"p = {genome_wide_threshold:.1e}",
        )
        ax.set_title(pvalue_column)
        ax.set_xlabel("Genomic position")
        ax.set_ylabel("-log10(p-value)")
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(chromosome_order, rotation=45, ha="right")
        ax.margins(x=0.01)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()

        plot_path = make_plot_path(output_prefix, pvalue_column)
        plot_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(plot_path, dpi=200)
        plt.close(fig)


def iter_variant_rows(
    vcf_path: Path,
    cohorts: list[tuple[str, list[int]]],
    missing_value: str,
    min_dp: int,
):
    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom, pos, _variant_id, ref, alt_field, _qual, _flt, _info, fmt = fields[:9]
            sample_fields = fields[9:]
            alt_alleles = alt_field.split(",")

            format_keys = fmt.split(":")
            if "GT" not in format_keys:
                raise ValueError(f"VCF record {chrom}:{pos} does not contain a GT field.")
            gt_index = format_keys.index("GT")
            dp_index = format_keys.index("DP") if "DP" in format_keys else None
            if min_dp > 0 and dp_index is None:
                raise ValueError(
                    f"VCF record {chrom}:{pos} does not contain a DP field, "
                    f"but --min-dp {min_dp} was requested."
                )

            cohort_counts: dict[str, dict[str, int | list[int]]] = {}
            for cohort_name, sample_indices in cohorts:
                alt_counts = [0] * len(alt_alleles)
                called_alleles = 0
                total_samples = len(sample_indices)
                usable_samples = 0

                for sample_index in sample_indices:
                    if sample_index >= len(sample_fields):
                        continue
                    sample_field = sample_fields[sample_index]
                    genotype = extract_format_value(sample_field, gt_index)
                    called = parse_called_alleles(genotype)
                    if not called:
                        continue

                    dp = parse_dp_value(extract_format_value(sample_field, dp_index))
                    if min_dp > 0 and (dp is None or dp < min_dp):
                        continue

                    usable_samples += 1
                    called_alleles += len(called)
                    for allele_index in called:
                        if allele_index > 0 and allele_index <= len(alt_alleles):
                            alt_counts[allele_index - 1] += 1

                cohort_counts[cohort_name] = {
                    "total_samples": total_samples,
                    "usable_samples": usable_samples,
                    "called_alleles": called_alleles,
                    "alt_counts": alt_counts,
                }

            for alt_index, alt in enumerate(alt_alleles):
                if len(ref) != 1 or len(alt) != 1:
                    continue

                row = {
                    "chrom": chrom,
                    "pos": pos,
                    "snp_id": f"{chrom}_{pos}_{ref}_{alt}",
                    "ref": ref,
                    "alt": alt,
                }

                for cohort_name, _sample_indices in cohorts:
                    total_samples = int(cohort_counts[cohort_name]["total_samples"])
                    usable_samples = int(cohort_counts[cohort_name]["usable_samples"])
                    called_alleles = int(cohort_counts[cohort_name]["called_alleles"])
                    alt_count = cohort_counts[cohort_name]["alt_counts"][alt_index]
                    row[f"AF_{cohort_name}"] = format_number(
                        allele_frequency(alt_count, called_alleles),
                        missing_value,
                    )
                    row[f"missingness_pct_{cohort_name}"] = format_number(
                        missingness_percentage(total_samples - usable_samples, total_samples),
                        missing_value,
                    )

                for (left_name, _), (right_name, _) in itertools.combinations(cohorts, 2):
                    left_called = int(cohort_counts[left_name]["called_alleles"])
                    right_called = int(cohort_counts[right_name]["called_alleles"])
                    left_alt = cohort_counts[left_name]["alt_counts"][alt_index]
                    right_alt = cohort_counts[right_name]["alt_counts"][alt_index]
                    left_ref = left_called - left_alt
                    right_ref = right_called - right_alt

                    row[f"oddsratio_{left_name}_vs_{right_name}"] = format_number(
                        odds_ratio(left_alt, left_ref, right_alt, right_ref),
                        missing_value,
                    )
                    p_value = fisher_exact_two_sided(
                        left_alt,
                        left_ref,
                        right_alt,
                        right_ref,
                    )
                    row[f"pvalue_{left_name}_vs_{right_name}"] = format_number(
                        p_value,
                        missing_value,
                    )

                yield row


def main() -> int:
    args = parse_args()

    vcf_path = Path(args.vcf).expanduser()
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")

    cohort_inputs = [parse_cohort_argument(argument) for argument in args.cohort]
    cohort_samples = [
        (cohort_name, load_cohort_samples(path, args.sample_column, args.has_header))
        for cohort_name, path in cohort_inputs
    ]

    vcf_samples = load_vcf_samples(vcf_path)
    cohorts = build_cohort_indices(vcf_samples, cohort_samples)

    output_path = Path(args.output).expanduser()
    header = ["chrom", "pos", "snp_id", "ref", "alt"]
    header.extend(f"AF_{cohort_name}" for cohort_name, _ in cohorts)
    header.extend(f"missingness_pct_{cohort_name}" for cohort_name, _ in cohorts)
    oddsratio_columns = [
        f"oddsratio_{left_name}_vs_{right_name}"
        for (left_name, _), (right_name, _) in itertools.combinations(cohorts, 2)
    ]
    header.extend(oddsratio_columns)
    pvalue_columns = [
        f"pvalue_{left_name}_vs_{right_name}"
        for (left_name, _), (right_name, _) in itertools.combinations(cohorts, 2)
    ]
    header.extend(pvalue_columns)

    plot_data = {column: [] for column in pvalue_columns}

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in iter_variant_rows(vcf_path, cohorts, args.missing_value, args.min_dp):
            writer.writerow(row)
            if plot_data:
                chrom = row["chrom"]
                pos = int(row["pos"])
                snp_id = row["snp_id"]
                for column in pvalue_columns:
                    pvalue = parse_numeric(row[column], args.missing_value)
                    if pvalue is not None:
                        plot_data[column].append((chrom, pos, snp_id, pvalue))

    if args.manhattan_prefix:
        draw_manhattan_plots(plot_data, args.manhattan_prefix, args.genome_wide_threshold)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
