import gzip
import re
import sys
from argparse import ArgumentError, ArgumentParser
from concurrent import futures
from pathlib import Path
from typing import *

__version__ = "0.9.1"


__doc__ = f"""XPCLR shell 脚本生成工具.

{Path(__file__).name} {__version__}: 用于生成运行 XPCLR 的 shell 脚本的 Python 脚本工具.

必需参数:
    --pop1 (字符串): pop1 的标识名, 例如: `sativa133`.
    --pop2 (字符串): pop2 的标识名, 例如: `serriola199`.

路径参数:
    --root-dir (路径): 下列子目录的根目录, 默认为当前目录.
    --data-dir (字符串): 原始 `*.vcf.gz` 数据目录名, 默认为 `00.data`.
    --filter-dir (字符串): filter 操作数据目录名, 默认为 `01.filter`.
    --overlap-dir (字符串): overlap 操作数据目录名, 默认为 `02.overlap`.
    --split-dir (字符串): split 操作数据目录名, 默认为 `03.split`.
    --genomap-dir (字符串): geno 和 map 操作数据目录名, 默认为 `04.genomap`.
    --xpclr-dir (字符串): xpclr 操作数据目录名, 默认为 `05.xpclr`.

可选参数:
    --chrid-list (路径): 一个文本文件, 有一列数据, 是需要进行 split 操作的染色体 id, 与 vcf 文件中保持一致. 默认为 `chrid.list`.
    --interval (整数): 进行 split 操作时, 每个划分文件最大的行数. 默认为 `200000`.

工具路径参数:
    --bcftools (路径): bcftools, 默认为 `bcftools`.
    --vcftools (路径): vcftools, 默认为 `vcftools`.
    --tabix (路径): tabix, 默认为 `tabix`.
    --perl (路径): perl, 默认为 `perl`.
    --perl-script (路径): 用于 geno 操作的 perl 脚本, 默认为 `/ldfssz1/ST_EARTH/P18Z10200N0148/P18Z10200N0148_LETTUCE/final_combine/dp2-50/29.xp-clr/vcf2geno.v20200716.pl`.
    --xpclr (路径): XPCLR, 默认为 `/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/liuxinjiang/liuxinjiang/APP/software/XPCLR/bin/XPCLR`

运行参数:
    下列参数指定生成哪些步骤的脚本, 可以组合使用.

    --run-all: 生成所有脚本.
    --run-filter: 生成 filter 操作脚本.
    --run-overlap: 生成 overlap 操作脚本.
    --run-split: 生成 split 操作脚本.
    --run-genomap: 生成 geno 和 map 操作脚本.
    --run-xpclr: 生成 xpclr 操作脚本.
"""


def get_chrpos_from_vcfgz(vcfgz_path: Path) -> Dict[str, List[str]]:
    """
    Returns:
        {"chrid": [1234, 5678, ...]}
    """

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


def get_chrpos_from_vcfgz_files(vcfgz_paths: List[Path]) -> Dict[str, List[str]]:
    """
    Returns:
        {"chrid": [1234, 5678, ...]}
    """

    tasks: List[futures.Future[Dict[str, List[str]]]] = []
    with futures.ProcessPoolExecutor(len(vcfgz_paths)) as executor:
        for p in vcfgz_paths:
            fu = executor.submit(get_chrpos_from_vcfgz, p)
            tasks.append(fu)

    chr_pos = {}
    for t in futures.as_completed(tasks):
        chr_pos.update(t.result())

    return chr_pos


class XPCLR:
    DATA_SUFFIX = ".vcf.gz"
    FILTER_SUFFIX = ".filter.vcf.gz"
    OVERLAP_SUFFIX = ".overlap.vcf.gz"
    SPLIT_SUFFIX = ".vcf.gz"
    GENO_SUFFIX = ".geno"
    MAP_SUFFIX = ".map"
    XPCLR_SUFFIX = ".xpclr.txt"

    def __init__(self, args) -> None:
        self.args = args

        # command line tools
        self.bcftools = Path(args.bcftools)
        self.vcftools = Path(args.vcftools)
        self.tabix = Path(args.tabix)
        self.perl = Path(args.perl)
        self.perl_script = Path(args.perl_script)
        self.xpclr = Path(args.xpclr)

        # run args
        self.pop1 = str(args.pop1)
        self.pop2 = str(args.pop2)

        # directories
        self.root_dir = Path(args.root_dir).resolve()
        self.data_dir = self.root_dir.joinpath(args.data_dir)
        self.filter_dir = self.root_dir.joinpath(args.filter_dir)
        self.overlap_dir = self.root_dir.joinpath(args.overlap_dir)
        self.split_dir = self.root_dir.joinpath(args.split_dir)
        self.genomap_dir = self.root_dir.joinpath(args.genomap_dir)
        self.xpclr_dir = self.root_dir.joinpath(args.xpclr_dir)

    def find_file_pairs(self, folder: Path, pattern: str) -> Tuple[List[Path], List[Path]]:
        """Return sorted file pairs."""

        print(f"Finding {pattern} files in {folder}")

        pop1_paths: List[Path] = []
        pop2_paths: List[Path] = []
        for p in folder.glob(pattern):
            if not p.is_file():
                continue

            if self.pop1 in p.name:
                pop1_paths.append(p)
            elif self.pop2 in p.name:
                pop2_paths.append(p)
            else:
                raise ValueError(f"Unknown file found: {p}")

        if len(pop1_paths) != len(pop2_paths):
            raise ValueError(f"{pattern} files count not equal, pop1: {len(pop1_paths)}, pop2: {len(pop2_paths)}")

        if len(pop1_paths) <= 0:
            raise ValueError(f"{pattern} files not found in {folder}")

        pop1_paths.sort()
        pop2_paths.sort()

        print("="*10, f"Found {len(pop1_paths)} pop1 {pattern} files", "="*10)
        print(*pop1_paths, sep="\n")

        print("="*10, f"Found {len(pop2_paths)} pop2 {pattern} files", "="*10)
        print(*pop2_paths, sep="\n")

    def _generate_script_filter(self, out_dir: Path, in_path: Path) -> Path:
        name = in_path.name[:-len(self.DATA_SUFFIX)]

        tmp_vcf_path = out_dir.joinpath(f"{name}.filter.vcf")
        filter_path = out_dir.joinpath(f"{name}{self.FILTER_SUFFIX}")

        with out_dir.joinpath(f"{name}.filter.sh").open("w", encoding="utf8") as f:
            print(f"{self.bcftools} filter -e 'F_MISSING > 0.5 || MAC < 2' {in_path} > {tmp_vcf_path}", file=f)
            print(f"{self.bcftools} view -Oz -o {filter_path} {tmp_vcf_path}", file=f)
            print(f"{self.tabix} -p vcf {filter_path}", file=f)
            print(f"rm {tmp_vcf_path}", file=f)

        return filter_path

    def generate_scripts_filter(self, pop1_data_paths: List[Path] = None, pop2_data_paths: List[Path] = None) -> Tuple[List[Path], List[Path]]:
        if pop1_data_paths is None or pop2_data_paths is None:
            pop1_data_paths, pop2_data_paths = self.find_file_pairs(self.data_dir, self.DATA_SUFFIX)

        filter_dir = self.filter_dir
        filter_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output filter scripts to {filter_dir}")

        pop1_filter_paths: List[Path] = []
        pop2_filter_paths: List[Path] = []
        for p1, p2 in zip(pop1_data_paths, pop2_data_paths):
            path = self._generate_script_filter(filter_dir, p1)
            pop1_filter_paths.append(path)

            path = self._generate_script_filter(filter_dir, p2)
            pop2_filter_paths.append(path)

        return pop1_filter_paths, pop2_filter_paths

    def generate_scripts_overlap(self, pop1_filter_paths: List[Path] = None, pop2_filter_paths: List[Path] = None) -> Tuple[List[Path], List[Path]]:
        if pop1_filter_paths is None or pop2_filter_paths is None:
            pop1_filter_paths, pop2_filter_paths = self.find_file_pairs(self.filter_dir, self.FILTER_SUFFIX)

        overlap_dir = self.overlap_dir
        overlap_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output overlap scripts to {overlap_dir}")

        pop1_overlap_paths: List[Path] = []
        pop2_overlap_paths: List[Path] = []
        for p1, p2 in zip(pop1_filter_paths, pop2_filter_paths):
            name1 = p1.name[:-len(self.FILTER_SUFFIX)]
            name2 = p2.name[:-len(self.FILTER_SUFFIX)]

            list_path1 = overlap_dir.joinpath(f"{name1}.list")
            list_path2 = overlap_dir.joinpath(f"{name2}.list")
            overlap_list_path = overlap_dir.joinpath(f"{name1}_{name2}.overlap.list")

            tmp_vcf_path1 = overlap_dir.joinpath(f"{name1}.overlap.recode.vcf")
            tmp_vcf_path2 = overlap_dir.joinpath(f"{name2}.overlap.recode.vcf")

            overlap_path1 = overlap_dir.joinpath(f"{name1}{self.OVERLAP_SUFFIX}")
            overlap_path2 = overlap_dir.joinpath(f"{name2}{self.OVERLAP_SUFFIX}")

            with overlap_dir.joinpath(f"{name1}_{name2}.overlap.sh").open("w", encoding="utf8") as f:
                print(f"{self.bcftools} query -f '%CHROM\\t%POS]\\n' {p1} > {list_path1}", file=f)
                print(f"{self.bcftools} query -f '%CHROM\\t%POS]\\n' {p2} > {list_path2}", file=f)
                print(f"awk 'NR==FNR{{a[$1,$2]; next}} ($1,$2) in a' {list_path1} {list_path2} > {overlap_list_path}", file=f)

                print(f"{self.vcftools} --gzvcf {p1} --positions {overlap_list_path} --recode --out {name1}.overlap", file=f)
                print(f"{self.bcftools} view -Oz -o {overlap_path1} {tmp_vcf_path1}", file=f)
                print(f"{self.tabix} -p vcf {overlap_path1}", file=f)
                print(f"rm {tmp_vcf_path1}", file=f)

                print(f"{self.vcftools} --gzvcf {p2} --positions {overlap_list_path} --recode --out {name2}.overlap", file=f)
                print(f"{self.bcftools} view -Oz -o {overlap_path2} {tmp_vcf_path2}", file=f)
                print(f"{self.tabix} -p vcf {overlap_path2}", file=f)
                print(f"rm {tmp_vcf_path2}", file=f)

            pop1_overlap_paths.append(overlap_path1)
            pop2_overlap_paths.append(overlap_path2)

        return pop1_overlap_paths, pop2_overlap_paths

    def _generate_script_split(self, out_dir: Path, in_path: Path, chrid: str, positions: List[int], interval: int) -> List[Path]:
        name = in_path.name[:-len(self.OVERLAP_SUFFIX)]

        split_paths = []
        with out_dir.joinpath(f"{chrid}.{name}.split.sh").open("w", encoding="utf8") as f:
            for lstart in range(0, len(positions), interval):
                lend = min(lstart + interval - 1, len(positions) - 1)
                pos_start, pos_end = positions[lstart], positions[lend]

                tmp_vcf_path = out_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name}.vcf")
                split_path = out_dir.joinpath(f"{chrid}-{lstart}-{lend}.{name}{self.SPLIT_SUFFIX}")

                print(f"{self.bcftools} filter {in_path} --regions {chrid}:{pos_start}-{pos_end} > {tmp_vcf_path}", file=f)
                print(f"{self.bcftools} view {tmp_vcf_path} -Oz -o {split_path}", file=f)
                print(f"{self.tabix} -p vcf {split_path}", file=f)
                print(f"rm {tmp_vcf_path}", file=f)

                split_paths.append(split_path)

        return split_paths

    def generate_scripts_split(self, pop1_overlap_paths: List[Path] = None, pop2_overlap_paths: List[Path] = None) -> Tuple[List[Path], List[Path]]:
        if pop1_overlap_paths is None or pop2_overlap_paths is None:
            pop1_overlap_paths, pop2_overlap_paths = self.find_file_pairs(self.overlap_dir, self.OVERLAP_SUFFIX)

        chrid_list = Path(self.args.chrid_list).read_text().strip().split("\n")
        interval = int(self.args.interval)

        split_dir = self.split_dir
        split_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output split scripts to {split_dir}")

        pop1_split_paths: List[Path] = []
        pop2_split_paths: List[Path] = []

        # extract positions from vcfgz files
        chrpos1 = get_chrpos_from_vcfgz_files(pop1_overlap_paths)
        chrpos2 = get_chrpos_from_vcfgz_files(pop2_overlap_paths)

        # preprocess different vcf.gz count
        _pop1_overlap_paths = pop1_overlap_paths
        _pop2_overlap_paths = pop2_overlap_paths
        pop1_overlap_paths: List[Path] = []
        pop2_overlap_paths: List[Path] = []

        if len(pop1_overlap_paths) <= 1:
            pop1_overlap_paths = pop1_overlap_paths * len(chrid_list)
            pop2_overlap_paths = pop2_overlap_paths * len(chrid_list)
        else:
            if len(_pop1_overlap_paths) > len(chrid_list):
                print(f"WARNING: {self.OVERLAP_SUFFIX} files {len(_pop1_overlap_paths)} more than chrom list row count {len(chrid_list)}", file=sys.stderr)

            for chrid in chrid_list:
                _p1 = list(filter(lambda x: chrid in x.name, _pop1_overlap_paths))
                _p2 = list(filter(lambda x: chrid in x.name, _pop2_overlap_paths))
                if len(_p1) <= 0 or len(_p1) > 1:
                    raise ValueError(f"{chrid} count incorrect in pop1 overlap.vcf.gz files, all pop1 files: {_pop1_overlap_paths}")
                if len(_p2) <= 0 or len(_p2) > 1:
                    raise ValueError(f"{chrid} count incorrect in pop2 overlap.vcf.gz files, all pop2 files: {_pop2_overlap_paths}")
                pop1_overlap_paths.append(_p1[0])
                pop2_overlap_paths.append(_p2[0])

        for chrid, p1, p2 in zip(chrid_list, pop1_overlap_paths, pop2_overlap_paths):
            paths = self._generate_script_split(split_dir, p1, chrid, chrpos1[chrid], interval)
            pop1_split_paths.extend(paths)

            paths = self._generate_script_split(split_dir, p2, chrid, chrpos2[chrid], interval)
            pop2_split_paths.extend(paths)

        return pop1_split_paths, pop2_split_paths

    def _generate_script_genomap(self, out_dir: Path, in_path: Path) -> Tuple[Path, Path]:
        name = in_path.name[:-len(self.SPLIT_SUFFIX)]

        geno_path = out_dir.joinpath(f"{name}.geno")
        map_path = out_dir.joinpath(f"{name}.map")

        with out_dir.joinpath(f"{name}.genomap.sh").open("w", encoding="utf8") as f:
            print(f"{self.perl} {self.perl_script} {in_path} | sed 's/|/\//g' | sed 's/\// /g' | sed 's/\./9/g' > {geno_path}", file=f)
            print(f"""{self.bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {in_path} | awk '{{print $1"_"$2"\\t9\\t"158.5/204289203*$2"\\t"$2"\\t"$3"\\t"$4}}' > {map_path}""", file=f)

        return geno_path, map_path

    def generate_scripts_genomap(self, pop1_split_paths: List[Path] = None, pop2_split_paths: List[Path] = None) -> Tuple[List[Path], List[Path], List[Path], List[Path]]:
        """
        Returns:
            (geno1, geno2, map1, map2)        
        """

        if pop1_split_paths is None or pop2_split_paths is None:
            pop1_split_paths, pop2_split_paths = self.find_file_pairs(self.split_dir, self.SPLIT_SUFFIX)

        genomap_dir = self.genomap_dir
        genomap_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output geno & map scripts to {genomap_dir}")

        pop1_geno_paths: List[Path] = []
        pop2_geno_paths: List[Path] = []
        pop1_map_paths: List[Path] = []
        pop2_map_paths: List[Path] = []

        for p1, p2 in zip(pop1_split_paths, pop2_split_paths):
            geno_path, map_path = self._generate_script_genomap(genomap_dir, p1)

            pop1_geno_paths.append(geno_path)
            pop1_map_paths.append(map_path)

            geno_path, map_path = self._generate_script_genomap(genomap_dir, p2)

            pop2_geno_paths.append(geno_path)
            pop2_map_paths.append(map_path)

        return pop1_geno_paths, pop2_geno_paths, pop1_map_paths, pop2_map_paths

    def generate_scripts_xpclr(self, pop1_geno_paths: List[Path] = None, pop2_geno_paths: List[Path] = None, pop_map_paths: List[Path] = None) -> List[Path]:
        if pop1_geno_paths is None or pop2_geno_paths is None:
            pop1_geno_paths, pop2_geno_paths = self.find_file_pairs(self.genomap_dir, self.GENO_SUFFIX)
        if pop_map_paths is None:
            pop_map_paths, _ = self.find_file_pairs(self.genomap_dir, self.MAP_SUFFIX)

        xpclr_dir = self.xpclr_dir
        xpclr_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output XPCLR scripts to {xpclr_dir}")

        xpclr_paths: List[Path] = []
        for g1, g2, m in zip(pop1_geno_paths, pop2_geno_paths, pop_map_paths):
            name1 = ".".join(g1.name[:-len(self.GENO_SUFFIX)].split(".")[1:])
            name2 = ".".join(g2.name[:-len(self.GENO_SUFFIX)].split(".")[1:])
            chr_interval = m.name.split(".")[0]
            chrnum = int(re.search("\d+", chr_interval.split("-")[0]).group())  # XPCLR need integer chrom number

            # NOTE: XPCLR Only accepts filenames in current directory, so link geno and map file to current directory
            _link = xpclr_dir.joinpath(g1.name)
            if not _link.is_file():
                _link.symlink_to(g1)
            _link = xpclr_dir.joinpath(g2.name)
            if not _link.is_file():
                _link.symlink_to(g2)
            _link = xpclr_dir.joinpath(m.name)
            if not _link.is_file():
                _link.symlink_to(m)

            xpclr_path = xpclr_dir.joinpath(f"{chr_interval}.{name1}_{name2}")  # Suffix ".xpclr.txt" will be auto added by XPCLR

            with xpclr_dir.joinpath(f"{chr_interval}.{name1}_{name2}.xpclr.sh").open("w", encoding="utf8") as f:
                print(f"{self.xpclr} -xpclr {g1.name} {g2.name} {m.name} {xpclr_path.name} -w1 0.005 100 2000 {chrnum} -p0 0.7", file=f)  # Only accept filename, not filepath

            xpclr_paths.append(xpclr_path)

        return xpclr_paths

    def run(self):
        pop1_filter_paths, pop2_filter_paths = None, None
        pop1_overlap_paths, pop2_overlap_paths = None, None
        pop1_split_paths, pop2_split_paths = None, None
        pop1_geno_paths, pop2_geno_paths = None, None
        pop1_map_paths, pop2_map_paths = None, None
        xpclr_paths = None

        if args.run_filter:
            pop1_filter_paths, pop2_filter_paths = self.generate_scripts_filter()
        if args.run_overlap:
            pop1_overlap_paths, pop2_overlap_paths = self.generate_scripts_overlap(pop1_filter_paths, pop2_filter_paths)
        if args.run_split:
            pop1_split_paths, pop2_split_paths = self.generate_scripts_split(pop1_overlap_paths, pop2_overlap_paths)
        if args.run_genomap:
            pop1_geno_paths, pop2_geno_paths, pop1_map_paths, pop2_map_paths = self.generate_scripts_genomap(pop1_split_paths, pop2_split_paths)
        if args.run_xpclr:
            xpclr_paths = self.generate_scripts_xpclr(pop1_geno_paths, pop2_geno_paths, pop1_map_paths)


if __name__ == "__main__":
    parser = ArgumentParser(description=f"{Path(__file__).name} {__version__}: 用于生成运行 XPCLR 的 shell 脚本的 Python 脚本工具.")

    parser.add_argument("--pop1", type=str, required=True, help="pop1 的标识名, 例如: `sativa133`.")
    parser.add_argument("--pop2", type=str, required=True, help="pop2 的标识名, 例如: `serriola199`.")

    parser.add_argument("--root-dir", type=Path, default=".", help="数据子目录的根目录, 默认为当前目录.")
    parser.add_argument("--data-dir", type=str, default="00.data", help="原始 `*.vcf.gz` 数据目录名, 默认为 `%(default)s`.")
    parser.add_argument("--filter-dir", type=str, default="01.filter", help="filter 操作数据目录名, 默认为 `%(default)s`.")
    parser.add_argument("--overlap-dir", type=str, default="02.overlap", help="overlap 操作数据目录名, 默认为 %(default)s`.")
    parser.add_argument("--split-dir", type=str, default="03.split", help="split 操作数据目录名, 默认为 `%(default)s`.")
    parser.add_argument("--genomap-dir", type=str, default="04.genomap", help="geno 和 map 操作数据目录名, 默认为 `%(default)s`.")
    parser.add_argument("--xpclr-dir", type=str, default="05.xpclr", help="xpclr 操作数据目录名, 默认为 `%(default)s`.")

    parser.add_argument("--chrid-list", type=Path, default="chrid.list", help="一个文本文件, 有一列数据, 是需要进行 split 操作的染色体 id, 与 vcf 文件中保持一致. 默认为 `%(default)s`.")
    parser.add_argument("--interval", type=int, default=200000, help="进行 split 操作时, 每个划分文件最大的行数. 默认为 `%(default)s`.")

    parser.add_argument("--bcftools", type=Path, default="bcftools", help="bcftools, 默认为 `%(default)s`.")
    parser.add_argument("--vcftools", type=Path, default="vcftools", help="vcftools, 默认为 `%(default)s`.")
    parser.add_argument("--tabix", type=Path, default="tabix", help="tabix, 默认为 `%(default)s`.")
    parser.add_argument("--perl", type=Path, default="perl", help="perl, 默认为 `%(default)s`.")
    parser.add_argument("--perl-script", type=Path, default="/ldfssz1/ST_EARTH/P18Z10200N0148/P18Z10200N0148_LETTUCE/final_combine/dp2-50/29.xp-clr/vcf2geno.v20200716.pl", help="用于 geno 操作的 perl 脚本, 默认为 `%(default)s`.")
    parser.add_argument("--xpclr", type=Path, default="/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/liuxinjiang/liuxinjiang/APP/software/XPCLR/bin/XPCLR", help="XPCLR, 默认为 `%(default)s`")

    parser.add_argument("--run-all", action="store_true", help="生成所有脚本.")
    parser.add_argument("--run-filter", action="store_true", help="生成 filter 操作脚本.")
    parser.add_argument("--run-overlap", action="store_true", help="生成 overlap 操作脚本.")
    parser.add_argument("--run-split", action="store_true", help="生成 split 操作脚本.")
    parser.add_argument("--run-genomap", action="store_true", help="生成 geno 和 map 操作脚本.")
    parser.add_argument("--run-xpclr", action="store_true", help="生成 xpclr 操作脚本.")

    args = parser.parse_args()

    if args.run_all:
        args.run_filter \
            = args.run_overlap \
            = args.run_split \
            = args.run_genomap \
            = args.run_xpclr \
            = True

    if not (args.run_filter or args.run_overlap or args.run_split or args.run_genomap or args.run_xpclr):
        raise ArgumentError(None, "No run arguments specifed, nothing to do.")

    xpclr = XPCLR(args)
    xpclr.run()
