import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import *

__date__ = "2024.01.06"

if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--pop1", type=str, required=True, help="Name of pop1, e.g. `sativa133`.")
    parser.add_argument("--pop2", type=str, required=True, help="Name of pop2, e.g. `serriola199`")

    parser.add_argument("--root-dir", type=Path, default=".", help="root directory of scripts folders.")
    parser.add_argument("--vcf-in", type=str, default="00.data", help="Input directory of pop `*.vcf.gz` files.")
    parser.add_argument("--filter-out", type=str, default="01.filter", help="Output directory of filter scripts.")
    parser.add_argument("--overlap-out", type=str, default="02.overlap", help="Output directory of overlap scripts.")

    parser.add_argument("--bcftools", type=Path, default="bcftools")
    parser.add_argument("--vcftools", type=Path, default="vcftools")
    parser.add_argument("--tabix", type=Path, default="tabix")

    args = parser.parse_args()

    bcftools: Path = args.bcftools
    vcftools: Path = args.vcftools
    tabix: Path = args.tabix

    root_dir = Path(args.root_dir).resolve()
    pop1_name: str = args.pop1
    pop2_name: str = args.pop2

    # 0. Find input files
    vcf_dir = root_dir.joinpath(args.vcf_in)
    print(f"Finding vcf.gz files in {vcf_dir}")

    pop1_vcf_paths: List[Path] = []
    pop2_vcf_paths: List[Path] = []
    for p in vcf_dir.glob("*.vcf.gz"):
        if not p.is_file():
            continue

        if pop1_name in p.name:
            pop1_vcf_paths.append(p)
        elif pop2_name in p.name:
            pop2_vcf_paths.append(p)
        else:
            raise ValueError(f"Unknown vcf.gz file found: {p}")

    if len(pop1_vcf_paths) != len(pop2_vcf_paths):
        raise ValueError(f"vcf.gz files not equal, pop1: {len(pop1_vcf_paths)}, pop2: {len(pop2_vcf_paths)}")

    if len(pop1_vcf_paths) <= 0:
        raise ValueError("vcf.gz files not found.")

    pop1_vcf_paths.sort()
    pop2_vcf_paths.sort()

    print("="*10, f"Found {len(pop1_vcf_paths)} pop1 vcf.gz files", "="*10)
    print(*pop1_vcf_paths, sep="\n")

    print("="*10, f"Found {len(pop2_vcf_paths)} pop2 vcf.gz files", "="*10)
    print(*pop2_vcf_paths, sep="\n")

    # 1. Filter
    filter_dir = root_dir.joinpath(args.filter_out)
    filter_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output filter scripts to {filter_dir}")

    pop1_filter_paths: List[Path] = []
    pop2_filter_paths: List[Path] = []
    for p1, p2 in zip(pop1_vcf_paths, pop2_vcf_paths):
        name1 = p1.name[:-len(".vcf.gz")]
        name2 = p2.name[:-len(".vcf.gz")]

        vcf_path1 = filter_dir.joinpath(f"{name1}.filter.vcf")
        vcfgz_path1 = filter_dir.joinpath(f"{name1}.filter.vcf.gz")

        with filter_dir.joinpath(f"{name1}.filter.sh").open("w", encoding="utf8") as f:
            print(f"{bcftools} filter -e 'F_MISSING > 0.5 || MAC < 2' {p1} > {vcf_path1}", file=f)
            print(f"{bcftools} view -Oz -o {vcfgz_path1} {vcf_path1}", file=f)
            print(f"{tabix} -p vcf {vcfgz_path1}", file=f)
            print(f"rm {vcf_path1}", file=f)
        pop1_filter_paths.append(vcfgz_path1)

        vcf_path2 = filter_dir.joinpath(f"{name2}.filter.vcf")
        vcfgz_path2 = filter_dir.joinpath(f"{name2}.filter.vcf.gz")

        with filter_dir.joinpath(f"{name2}.filter.sh").open("w", encoding="utf8") as f:
            print(f"{bcftools} filter -e 'F_MISSING > 0.5 || MAC < 2' {p2} > {vcf_path2}", file=f)
            print(f"{bcftools} view -Oz -o {vcfgz_path2} {vcf_path2}", file=f)
            print(f"{tabix} -p vcf {vcfgz_path2}", file=f)
            print(f"rm {vcf_path2}", file=f)
        pop2_filter_paths.append(vcfgz_path2)

    # 2. Overlap
    overlap_dir = root_dir.joinpath(args.overlap_out)
    overlap_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output overlap scripts to {overlap_dir}")

    pop1_overlap_paths: List[Path] = []
    pop2_overlap_paths: List[Path] = []
    for p1, p2 in zip(pop1_filter_paths, pop2_filter_paths):
        name1 = p1.name[:-len(".filter.vcf.gz")]
        name2 = p2.name[:-len(".filter.vcf.gz")]

        list_path1 = overlap_dir.joinpath(f"{name1}.list")
        list_path2 = overlap_dir.joinpath(f"{name2}.list")
        overlap_list = overlap_dir.joinpath(f"{name1}_{name2}.overlap.list")

        overlap_path1 = overlap_dir.joinpath(f"{name1}.overlap.vcf.gz")
        overlap_path2 = overlap_dir.joinpath(f"{name2}.overlap.vcf.gz")

        with overlap_dir.joinpath(f"{name1}_{name2}.overlap.sh").open("w", encoding="utf8") as f:
            print(f"{bcftools} query -f '%CHROM\\t%POS]\\n' {p1} > {list_path1}", file=f)
            print(f"{bcftools} query -f '%CHROM\\t%POS]\\n' {p2} > {list_path2}", file=f)
            print(f"awk 'NR==FNR{{a[$1,$2]; next}} ($1,$2) in a' {list_path1} {list_path2} > {overlap_list}", file=f)

            print(f"{vcftools} --gzvcf {p1} --positions {overlap_list} --recode --out {name1}.overlap", file=f)
            print(f"{bcftools} view -Oz -o {overlap_path1} {name1}.overlap.recode.vcf", file=f)
            print(f"rm {name1}.overlap.recode.vcf", file=f)
            print(f"{tabix} -p vcf {overlap_path1}", file=f)

            print(f"{vcftools} --gzvcf {p2} --positions {overlap_list} --recode --out {name2}.overlap", file=f)
            print(f"{bcftools} view -Oz -o {overlap_path2} {name2}.overlap.recode.vcf", file=f)
            print(f"rm {name2}.overlap.recode.vcf", file=f)
            print(f"{tabix} -p vcf {overlap_path2}", file=f)

        pop1_overlap_paths.append(overlap_path1)
        pop2_overlap_paths.append(overlap_path2)
