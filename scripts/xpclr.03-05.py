import gzip
import sys
from argparse import ArgumentParser
from pathlib import Path
from typing import *
import re

__date__ = "2024.01.06"


def get_chrpos_from_vcfgz(vcfgz_path: Path) -> Dict[str, List[str]]:
    chr_pos = {}
    last_chrid = None
    positions = []
    with gzip.open(vcfgz_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue

            chrid, pos, *_ = line.split("\t")
            if chrid != last_chrid:
                if last_chrid is not None and len(positions) > 0:
                    chr_pos[last_chrid] = positions
                last_chrid = chrid
                positions = [pos]
            else:
                positions.append(pos)

    if last_chrid is not None and len(positions) > 0:
        chr_pos[last_chrid] = positions

    return chr_pos


def get_chrnum(chrid: str) -> str:
    return re.search(r"\d+", chrid).group()


if __name__ == "__main__":
    parser = ArgumentParser()

    parser.add_argument("--pop1", type=str, required=True, help="Name of pop1, e.g. `sativa133`.")
    parser.add_argument("--pop2", type=str, required=True, help="Name of pop2, e.g. `serriola199`")

    parser.add_argument("--root-dir", type=Path, default=".", help="root directory of scripts folders.")
    parser.add_argument("--overlap-in", type=str, default="02.overlap", help="Input directory of overlap `*.overlap.vcf.gz` files.")
    parser.add_argument("--split-out", type=str, default="03.split", help="Output directory of split scripts.")
    parser.add_argument("--genomap-out", type=str, default="04.genomap", help="Output directory of genomap scripts.")
    parser.add_argument("--xpclr-out", type=str, default="05.xpclr", help="Output directory or XPCLR scripts.")

    parser.add_argument("--chr-list", type=Path, required=True, help="List file of chrom ids.")
    parser.add_argument("--interval", type=int, default=200000, help="Invertal of split split.")

    parser.add_argument("--bcftools", type=Path, default="bcftools")
    parser.add_argument("--tabix", type=Path, default="tabix")
    parser.add_argument("--perl", type=Path, default="perl")
    parser.add_argument("--perl-script", type=Path, default="/ldfssz1/ST_EARTH/P18Z10200N0148/P18Z10200N0148_LETTUCE/final_combine/dp2-50/29.xp-clr/vcf2geno.v20200716.pl")
    parser.add_argument("--xpclr", type=Path, default="/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/liuxinjiang/liuxinjiang/APP/software/XPCLR/bin/XPCLR")

    args = parser.parse_args()

    bcftools: Path = args.bcftools
    tabix: Path = args.tabix
    perl: Path = args.perl
    perl_script: Path = args.perl_script
    xpclr: Path = args.xpclr

    root_dir = Path(args.root_dir).resolve()
    pop1_name: str = args.pop1
    pop2_name: str = args.pop2

    # 0. Find input files
    overlap_dir = root_dir.joinpath(args.overlap_in)
    print(f"Finding overlap.vcf.gz files in {overlap_dir}")

    pop1_overlap_paths: List[Path] = []
    pop2_overlap_paths: List[Path] = []
    for p in overlap_dir.glob("*.overlap.vcf.gz"):
        if not p.is_file():
            continue

        if pop1_name in p.name:
            pop1_overlap_paths.append(p)
        elif pop2_name in p.name:
            pop2_overlap_paths.append(p)
        else:
            raise ValueError(f"Unknown overlap.vcf.gz file found: {p}")

    if len(pop1_overlap_paths) != len(pop2_overlap_paths):
        raise ValueError(f"overlap.vcf.gz files not equal, pop1: {len(pop1_overlap_paths)}, pop2: {len(pop2_overlap_paths)}")

    if len(pop1_overlap_paths) <= 0:
        raise ValueError("overlap.vcf.gz files not found.")

    pop1_overlap_paths.sort()
    pop2_overlap_paths.sort()

    print("="*10, f"Found {len(pop1_overlap_paths)} pop1 overlap.vcf.gz files", "="*10)
    print(*pop1_overlap_paths, sep="\n")

    print("="*10, f"Found {len(pop2_overlap_paths)} pop2 overlap.vcf.gz files", "="*10)
    print(*pop2_overlap_paths, sep="\n")

    # 3. Split
    chr_list = Path(args.chr_list).read_text().strip().split("\n")
    interval: int = args.interval

    split_dir = root_dir.joinpath(args.split_out)
    split_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output split scripts to {split_dir}")

    pop1_split_paths: List[Path] = []
    pop2_split_paths: List[Path] = []

    # preprocess different vcf.gz count
    _pop1_overlap_paths = pop1_overlap_paths
    _pop2_overlap_paths = pop2_overlap_paths
    pop1_overlap_paths: List[Path] = []
    pop2_overlap_paths: List[Path] = []

    if len(_pop1_overlap_paths) <= 1:
        pop1_overlap_paths = _pop1_overlap_paths * len(chr_list)
        pop2_overlap_paths = _pop2_overlap_paths * len(chr_list)
    else:
        if len(_pop1_overlap_paths) > len(chr_list):
            print(f"WARNING: overlap.vcf.gz files {len(_pop1_overlap_paths)} more than chrom list row count {len(chr_list)}", file=sys.stderr)

        for chrid in chr_list:
            _p1 = list(filter(lambda x: chrid in x.name, _pop1_overlap_paths))
            _p2 = list(filter(lambda x: chrid in x.name, _pop2_overlap_paths))
            if len(_p1) <= 0 or len(_p1) > 1:
                raise ValueError(f"{chrid} count incorrect in pop1 overlap.vcf.gz files, all pop1 files: {_pop1_overlap_paths}")
            if len(_p2) <= 0 or len(_p2) > 1:
                raise ValueError(f"{chrid} count incorrect in pop2 overlap.vcf.gz files, all pop2 files: {_pop2_overlap_paths}")
            pop1_overlap_paths.append(_p1[0])
            pop2_overlap_paths.append(_p2[0])

    for chrid, p1, p2 in zip(chr_list, pop1_overlap_paths, pop2_overlap_paths):
        name1 = p1.name[:-len(".overlap.vcf.gz")]
        name2 = p2.name[:-len(".overlap.vcf.gz")]
        chrpos1 = get_chrpos_from_vcfgz(p1)
        chrpos2 = get_chrpos_from_vcfgz(p2)

        with split_dir.joinpath(f"{chrid}.{name1}.split.sh").open("w", encoding="utf8") as f:
            positions = chrpos1[chrid]
            for lstart in range(0, len(positions), interval):
                lend = min(lstart + interval - 1, len(positions) - 1)
                pos_start, pos_end = positions[lstart], positions[lend]

                vcf_path = split_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name1}.vcf")
                split_path = split_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name1}.vcf.gz")

                print(f"{bcftools} filter {p1} --regions {chrid}:{pos_start}-{pos_end} > {vcf_path}", file=f)
                print(f"{bcftools} view {vcf_path} -Oz -o {split_path}", file=f)
                print(f"rm {vcf_path}", file=f)
                print(f"{tabix} -p vcf {split_path}", file=f)

                pop1_split_paths.append(split_path)

        with split_dir.joinpath(f"{chrid}.{name2}.split.sh").open("w", encoding="utf8") as f:
            positions = chrpos2[chrid]
            for lstart in range(0, len(positions), interval):
                lend = min(lstart + interval - 1, len(positions) - 1)
                pos_start, pos_end = positions[lstart], positions[lend]

                vcf_path = split_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name2}.vcf")
                split_path = split_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name2}.vcf.gz")

                print(f"{bcftools} filter {p2} --regions {chrid}:{pos_start}-{pos_end} > {vcf_path}", file=f)
                print(f"{bcftools} view {vcf_path} -Oz -o {split_path}", file=f)
                print(f"rm {vcf_path}", file=f)
                print(f"{tabix} -p vcf {split_path}", file=f)

                pop2_split_paths.append(split_path)

    # 4. Geno & Map
    genomap_dir = root_dir.joinpath(args.genomap_out)
    genomap_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output geno & map scripts to {genomap_dir}")

    pop1_geno_paths: List[Path] = []
    pop2_geno_paths: List[Path] = []
    pop1_map_paths: List[Path] = []
    pop2_map_paths: List[Path] = []
    for p1, p2 in zip(pop1_split_paths, pop2_split_paths):
        name1 = p1.name[:-len(".vcf.gz")]
        name2 = p2.name[:-len(".vcf.gz")]

        geno_path1 = genomap_dir.joinpath(f"{name1}.geno")
        geno_path2 = genomap_dir.joinpath(f"{name2}.geno")
        map_path1 = genomap_dir.joinpath(f"{name1}.map")
        map_path2 = genomap_dir.joinpath(f"{name2}.map")

        with genomap_dir.joinpath(f"{name1}.genomap.sh").open("w", encoding="utf8") as f:
            print(f"{perl} {perl_script} {p1} | sed 's/\// /g' | sed 's/\./9/g' > {geno_path1}", file=f)
            print(f"""{bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {p1} | awk '{{print $1"_"$2"\\t9\\t"158.5/204289203*$2"\\t"$2"\\t"$3"\\t"$4}}' > {map_path1}""", file=f)

        with genomap_dir.joinpath(f"{name2}.genomap.sh").open("w", encoding="utf8") as f:
            print(f"{perl} {perl_script} {p2} | sed 's/\// /g' | sed 's/\./9/g' > {geno_path2}", file=f)
            print(f"""{bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {p2} | awk '{{print $1"_"$2"\\t9\\t"158.5/204289203*$2"\\t"$2"\\t"$3"\\t"$4}}'> {map_path2}""", file=f)

        pop1_geno_paths.append(geno_path1)
        pop2_geno_paths.append(geno_path2)
        pop1_map_paths.append(map_path1)
        pop2_map_paths.append(map_path2)

    # 5. XPCLR
    xpclr_dir = root_dir.joinpath(args.xpclr_out)
    xpclr_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output XPCLR scripts to {xpclr_dir}")

    xpclr_paths: List[Path] = []
    for g1, g2, m1 in zip(pop1_geno_paths, pop2_geno_paths, pop1_map_paths):
        name1 = ".".join(g1.name[:-len(".geno")].split(".")[1:])
        name2 = ".".join(g2.name[:-len(".geno")].split(".")[1:])
        chr_interval = m1.name.split(".")[0]
        chrid = chr_interval.split("-")[0]

        xpclr_path = xpclr_dir.joinpath(f"{chr_interval}.{name1}_{name2}.xpclr.txt")

        with xpclr_dir.joinpath(f"{chr_interval}.{name1}_{name2}.xpclr.sh").open("w", encoding="utf8") as f:
            print(f"{xpclr} -xpclr {g1} {g2} {m1} {xpclr_path} -w1 0.005 100 2000 {chrid} -p0 0.7", file=f)

        xpclr_paths.append(xpclr_path)
