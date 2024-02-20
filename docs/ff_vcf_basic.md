---
author: 松下问童子
date: 2023-02-20
---

# 文件格式 | 对VCF文件的常见操作

## 压缩

```bash
bcftools view sativa332.indel.Lsat_1_v5_gn_7_15020.hard.vcf -Oz -o sativa332.indel.Lsat_1_v5_gn_7_15020.hard.vcf.gz
```

- `-Oz`: 指定输出文件为压缩格式 `gz`.
- `-o`: 指定输出文件名。

## 建立索引

```bash
bcftools index -t sativa332.indel.Lsat_1_v5_gn_7_15020.hard.vcf.gz
```

- `-t`: 填压缩后的 `vcf.gz` 文件.

---

```bash
tabix -p vcf view.vcf.gz
```

- `-p`: 填压缩后的 `vcf.gz` 文件.

---

上述命令二选一, 均可完成压缩.

## 提取变异

```bash
gatk SelectVariants -R Lactuca_sativa_Salinas_V8.fa -V Lactuca.snp.vcf.gz -L chr7:18846094-18857694 -L chr8:1884094-188009694 -O sativa332.snp.Lsat_1_v5_gn_7_15020.hard.vcf.gz
```

- `-R`: 参考基因组fasta文件
- `-V`: 输入vcf文件
- `-L`: 需要提取的变异区间,可以为 `chrx:xx-xx` 或直接提取整条染色体 `chrx`, 但格式必须与输入vcf文件相同,如: `chr1`, `Chr1`, `Chr01`, 可以同时提取多个区间 (重复使用该参数)
- `-O`: 输出vcf文件名

---

```bash
bcftools view -S keep.list test.vcf > sub_indv.vcf
```

- `-S`: 需要提取的变异的位置信息,可以是 `染色体\t具体位置` 两列，也可以是 `染色体\t起始\t终止` 三列.

---

上述命令二选一均可完成

## 提取样品

```bash
gatk SelectVariants -R Lactuca_sativa_Salinas_V8.fa -V Lactuca.snp.vcf.gz -sn s001 -sn s002 -sn s003 -O sativa332.snp.Lsat_1_v5_gn_7_15020.hard.vcf
```

- `-R`: 参考基因组fasta文件
- `-V`: 输入vcf文件
- `-O`: 输出vcf文件名
- `-sn`: 需要提取的样品名,可以同时提取多个样品 (重复使用该参数)

---

```bash
bcftools view view.vcf.gz -s NA00001,NA00002  -o subset.vcf 
```

- `-s`: 需要提取的样品名,`,`分割
- `-o`: 输出文件名

---

```bash
bcftools view input.SNP.revhet.recode.vcf.gz -S input_samples.list -Oz -o output_fixed.vcf.gz
```

- `-S`: 需要提取的样品名称列表
  
```bash
vcftools --gzvcf in.vcf.gz --recode --recode-INFO-all --stdout --keep id.txt  > out.vcf.gz
```

- `--gzvcf`: 输入文件名
- `--recode`: 告诉程序使用过滤器写入一个新的vcf文件
- `--recode-INFO-all`: 保留旧vcf文件中的所有INFO标志
- `--stdout`: 配合后面的管道函数，如果输出文件是压缩文件 `--stdout > out.vcf.gz` (也可以不压缩),不压缩时，可以直接 --out out 此处第一个out是命令，第二个out是文件名 (也可以 `--stdout > out.vcf`)
- `--keep`: 所需要的样品名列表

---

选择一个即可完成,较多样品时推荐选择 `vcftools`, 可以直接提供 list

## 去除样品

```bash
bcftools view example.vcf.gz -s ^NA00001,NA00002  -o subset.vcf 
```

- `-s` `^` 后接需要去除的样品名, `,` 分割
- `-o` 输出文件名

```bash
vcftools --gzvcf in.vcf.gz --recode --recode-INFO-all --stdout  --remove id.txt  > out.vcf.gz
```

- `--gzvcf`: 输入文件名
- `--recode`: 告诉程序使用过滤器写入一个新的vcf文件
- `--recode-INFO-all`: 保留旧vcf文件中的所有INFO标志
- `--stdout`: 配合后面的管道函数，如果输出文件是压缩文件 `--stdout > out.vcf.gz` (也可以不压缩), 不压缩时，可以直接 `--out out` 此处第一个 `out` 是命令，第二个 `out` 是文件名 (也可以 `--stout > out.vcf`)
- `--remove`: 所需要的样品名列表

---

## 修改样品名

[bcftools的常用命令](https://www.jianshu.com/p/8ddee965e412)

```bash
bcftools reheader --samples sub.list -o output.vcf.gz input.vcf.gz
```

- `--samples`: 内容为`oldname newname`,可有多行(修改多个样品名)


---

## MAF 过滤

参考: [VCFtools的使用(参数说明)](https://www.jianshu.com/p/c16dfe358596)

```bash
/zfssz5/BC_PUB/Software/03.Soft_ALL/vcftools-0.1.13/bin/vcftools --gzvcf sativa332.hard.Lsat_1_v5_gn_7_15020.snp_indel.vcf.gz --max-missing 0.9 --maf 0.05 --recode --recode-INFO-all --stdout | /ldfssz1/ST_AGRIC/P18Z10200N0148_LETTUCE/USER/bin/01.software/smcpp1.13.1/bin/bgzip -c > sativa332.maf.Lsat_1_v5_gn_7_15020.snp_indel.vcf.gz
```

- `--gzvcf`: 输入文件名
- `--max-missing`: 缺失率，0为接受完全缺失，1为接受数据全都存在
- `--maf`: Minor Allele Frequency二等位基因频率进行过滤，通常保留大于0.05的。
- `--recode` : 告诉程序使用过滤器写入一个新的vcf文件
- `--recode-INFO-all`: 保留旧vcf文件中的所有INFO标志
- `--stdout`: 配合后面的管道函数，如果输出文件是压缩文件 `--stdout > out.vcf.gz` (也可以不压缩), 不压缩时，可以直接 `--out out` 此处第一个 `out` 是命令，第二个 `out` 是文件名 (也可以 `--stout > out.vcf`)
- `bgzip -c >`: 压缩文件

## 提取染色体名称

```bash
bioawk -c vcf ' { print $chrom } ' old.vcf.gz | sort -n | uniq > ChrName.txt
```

- `sort -n` :按数值排序
- `uniq`: 去重

## 更改染色体名称

```bash
bcftools annotate --rename-chrs NewChrName.txt old.xxx.vcf.gz -Oz -o new.xxx.gz 
```

- `rename-chr`: 接新老名称的对应关系,其中NewChrName.txt文件中第一列为旧ID，第二列为对应的新ID
- `-Oz`: 指定输出文件为压缩格式 `gz`.
- `-o`: 指定输出文件名。

## 按染色体拆分vcf

```bash
for i in {chr1..chr9};do echo -e "gatk SelectVariants -R Lactuca_sativa_Salinas_V8.fa -V Lactuca.snp.vcf.gz -L chr7:18846094-18857694 -L $i -O sativa332.snp.$i.hard.vcf.gz";done
```

循环使用提取变异命令

## 按染色体合并vcf (推荐转换成bed文件合并，更快速)

参考: [多款软件进行vcf合并-gatk、vcftools、bcftools](https://www.jianshu.com/p/116c8e4c93f4)

```bash
gatk GatherVcfs -I concat-a.vcf -I concat-b.vcf -O combine_a_b_samesample_diffsites.vcf
```

- `-I`: 输入文件名(可重复使用该参数)
- `-O`: 输出文件名

```bash
gatk MergeVcfs -I concat-a.vcf -I concat-b.vcf -O combine_a_b_diffsample_allsites_gatk.vcf
```

- `-I`: 输入文件名(可重复使用该参数)
- `-O`: 输出文件名

```bash
bcftools concat concat-a.vcf concat-b.vcf -o combine_a_b_samesample_diffsites_bcftools.vcf
```

- `concat`:后接需要合并的vcf名
- `-O`: 输出文件名

## 合并SNP和indel

```bash
bcftools concat -R Lsat_1_v5_gn_7_15020.region -a sativa332.indel.hard.vcf.gz sativa332.snp.vcf.gz -Oz -o sativa332.hard.snp_indel.vcf.gz
```

- `-R`:  只有包含在这个文件中指定的区域的记录才会被合并,格式为 `chr1\t1389\t111333`,可有多行(多个区间)
- `-a`: 需要合并的文件,示例中的两个文件将按顺序合并
- `-Oz`: 指定输出文件为压缩格式 `gz`.
- `-o`: 指定输出文件名。

## vcf的基本信息

```bash
bcftools stats example.vcf.gz >  example.vcf.stats
```

## 提取样品ID信息

```bash
bcftools query -l example.vcf.gz >example.id.list
```

- `-l`: 后接需要查询的vcf文件

---

## 提取变异信息

参考: [bcftools常用命令详解](https://www.jianshu.com/p/266c55c87978)

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n' view.vcf.gz
```

- `%CHROM`: 代表`VCF`文件中染色体那一列，其他的列，比如`POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`也是类似的写法
- `[]`: 对于`FORMAT`字段的信息，必须要中括号括起来
- `%SAMPLE`: 代表样本名称
- `%GT`: 代表`FORMAT`字段中`genotype`的信息
