"""该文件用来解决 `fneighbor` ID 不能超过 10 个字符的问题."""

# -*- coding: UTF-8 -*-

import sys
from argparse import ArgumentParser
from hashlib import md5


def submatrix(mat_path: str, out_path: str, mapping_path: str):
    id_prefix = "_" + md5(mat_path.encode("utf8")).hexdigest().upper()[:4]  # 取前 4 个

    with open(mat_path, "r", encoding="utf8") as f_mat:
        with open(out_path, "w", encoding="utf8") as f_out:
            with open(mapping_path, "w", encoding="utf8") as f_mapping:
                count_line = f_mat.readline()
                f_out.write(count_line)

                count = int(count_line.rstrip("\n"))
                print("Matrix line count: {}".format(count))

                for i in range(count):
                    line = f_mat.readline()
                    sid, nums = line.split("\t")
                    sid2 = "{}{:04d}".format(id_prefix, i)

                    f_out.write(sid2)
                    f_out.write(nums)
                    print(sid2, sid, sep="\t", file=f_mapping)


def restoretree(tree_path: str, mapping_path: str, out_tree_path: str):
    with open(tree_path, "r", encoding="utf8") as f_tree:
        tree_text = f_tree.read()
    with open(mapping_path, "r", encoding="utf8") as f_mapping:
        for line in f_mapping:
            sid2, sid = line.rstrip("\n").split("\t")
            tree_text = tree_text.replace(sid2, sid)
    with open(out_tree_path, "w", encoding="utf8") as f_out:
        f_out.write(tree_text)


if __name__ == "__main__":
    parser = ArgumentParser(description="用来临时替换 njtree 输入文件")

    parser.add_argument("cmd", type=str, choices=["submatrix", "restoretree"],
                        help="`submatrix` 用来替换输入矩阵文件, `restoretree` 用来还原输出的树结构文件")

    parser.add_argument("--mat", type=str, help="要被替换 ID 的矩阵输入文件")
    parser.add_argument("--tree", type=str, help="要被还原 ID 的 njtree 结构输入文件")

    parser.add_argument("--mapping", type=str, default="matrix.cols.mapping", help="临时矩阵文件的 ID 的对应关系文件")
    parser.add_argument("--out-mat", type=str, default="matrix.tmp", help="替换后的临时矩阵输出文件")
    parser.add_argument("--out-tree", type=str, default="njtree.restore", help="还原后的 njtree 结构输出文件")

    args = parser.parse_args()

    cmd: str = args.cmd

    mat_path: str = args.mat
    tree_path: str = args.tree
    out_mat_path: str = args.out_mat
    mapping_path: str = args.mapping
    out_tree_path: str = args.out_tree

    if cmd == "submatrix" and mat_path:
        submatrix(mat_path, out_mat_path, mapping_path)
    elif cmd == "restoretree":
        restoretree(tree_path, mapping_path, out_tree_path)
    else:
        print("无事发生, 请检查 `cmd` 参数和对应的选项参数.", file=sys.stderr)
