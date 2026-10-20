# MagCluster <img src='https://raw.githubusercontent.com/RunJiaJi/MagCluster/main/docs/_static/logo_new.svg' align='right' width=220 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.1-orange?style=flat-square)](https://github.com/RunJiaJi/magcluster)

[English](README.md) | [中文说明](README_zh.md)

**MagCluster** 是一个用于在细菌基因组中**识别、注释和可视化磁小体基因簇（MGCs）**的生信工具。

---

## v1.0 做了哪些更新？

- **一条命令搞定**：新版用 `magcluster run` 整合了以前繁琐的多个步骤，一键完成基因预测、HMM 比对、基因簇识别和绘图。
- **不再依赖外部 Prokka**：内置原生 [Pyrodigal](https://github.com/althonos/pyrodigal) 进行快速基因预测，不需要再配置复杂的 Perl 环境和 BLAST 索引。
- **92 个精细 HMM 模型**：覆盖 8 大磁小体蛋白家族（`Mam`、`Mms`、`Mad`、`Man`、`Mae`、`Maq`、`Mca`、`Feo`），使用 [PyHMMER](https://github.com/althonos/pyhmmer) 进行快速检索。
- **支持多种输入格式**：支持带注释的 GenBank 文件（`.gbk`/`.gbff`）、未注释的核酸 Contig（`.fasta`/`.fna`/`.fa`），以及纯蛋白序列（`.faa`），同时支持 `.gz` 压缩文件。
- **每个基因组独立出结果**：每个样本单独生成 `{genome}_mgc_genes.tsv` 和 `{genome}_magnetosome_proteins.faa`，不再混在一起。
- **运行日志**（v1.0.1）：每次运行都会在输出目录保存 `magcluster.log`，事后可以回看参数、软件版本和每个基因组的结果。更新内容见 [CHANGELOG.md](https://github.com/RunJiaJi/MagCluster/blob/main/CHANGELOG.md)。

---

## 安装方法

需要 Python 3.10 或更高版本。

### 方式 A：pip（PyPI）

```bash
pip install magcluster
```

### 方式 B：conda / mamba（Bioconda）

```bash
conda create -n magcluster -c conda-forge -c bioconda magcluster
conda activate magcluster
```

装进已有环境：`conda install -c conda-forge -c bioconda magcluster`。请按上面的顺序把 `conda-forge` 写在 `bioconda` 前面。如果 conda 配置了严格的频道优先级（strict），只写 `conda install -c bioconda magcluster` 可能会无解。

### 方式 C：从 GitHub 安装开发版

```bash
pip install git+https://github.com/RunJiaJi/MagCluster.git
# 或者克隆后本地安装：
git clone https://github.com/RunJiaJi/MagCluster.git
cd MagCluster
pip install -e ".[dev]"
```

安装完成后检查版本：
```bash
magcluster --version
magcluster run --help      # 也可以用 -h
```

---

## 快速上手

### 1. 分析带注释的 GenBank 基因组（`.gbff` / `.gbk`）
```bash
magcluster run genome.gbff -o results/
```

### 2. 分析未注释的 FASTA Contig（`.fasta` / `.fna` / `.fa`）
程序会自动调用 Pyrodigal 进行基因预测，默认调用所有 CPU 核心：
```bash
magcluster run contigs.fasta -o results/ -t 8
```

### 3. 批量处理一个文件夹里的所有基因组
```bash
magcluster run genomes_dir/ -o results/
```

### 4. 分析蛋白质序列文件（`.faa` / `.pep`）
Prodigal 格式的表头（`>id # start # end # strand # ...`）带有基因组坐标，可以正常成簇和绘图。如果输入的是没有基因组坐标的纯蛋白 FASTA，程序会完成磁小体蛋白注释并输出 TSV 表格和 FASTA 序列文件，不进行基因簇绘图：
```bash
magcluster run proteins.faa -o results/
```

### 5. 覆盖已有结果目录
如果输出目录已存在且里面有文件，加上 `--force` 即可：
```bash
magcluster run genome.gbff -o results/ --force
```

---

## 常用命令与参数说明

### `magcluster run`

```bash
magcluster run [选项] 输入文件/文件夹...
```

`输入` 可以是文件或文件夹。文件夹里只读取带基因组扩展名的文件（`.gbk`、`.gbff`、`.gbf`、`.gb`、`.fa`、`.fasta`、`.fna`、`.fsa`、`.fsa_nt`、`.faa`、`.pep`，都可以再加 `.gz`）。其他文件（例如 `README.txt`）和隐藏文件会跳过，并给出警告。直接写在命令行里的文件没有扩展名也能用，程序会根据第一行判断格式。不支持 GFF3，因为 GFF3 里没有序列，请改用对应的 GenBank/GBFF、FASTA 或 Prodigal `.faa` 文件。

| 参数 | 别名 | 默认值 | 作用说明 |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | `magcluster_results` | 输出结果存放目录（不能是已存在的普通文件） |
| `--force` | `-f` | `False` | 允许使用非空的输出目录。会先删除该目录里上一次 MagCluster 运行的输出（表格、FASTA、`MGC_align.html`、`clusters_gbk/`、`figures/`、`magcluster.log`），其他文件保留 |
| `--threads` | `-t` | `0` | CPU 线程数（`0` 表示使用全部核心；超过可用核心数时自动降到可用核心数） |
| `--evalue` | `-e` | `1e-5` | HMM 比对的 E-value 阈值（≥ 0） |
| `--score` | `-s` | 各模型的 `trusted_cutoff`（25.0） | HMM 比对的最小 bit score。不指定时使用 `rules.yaml` 里每个模型的 `trusted_cutoff`（92 个内置模型都是 25.0）。指定后对所有模型统一使用这个值，也可以低于 25 |
| `--max-gap` | | `5` | 基因簇内两个磁小体基因之间允许插入的最大非磁小体基因数（≥ 0）。另有固定的 25 kb 距离上限 |
| `--threshold` | `-th`、`--min-families` | `3` | 判定一个 MGC 所需的最少磁小体家族数（任意非 Feo 家族都算，≥ 1） |
| `--flank` | `-un`、`--unmag-num` | `3` | 基因簇两端各保留的非磁小体基因数（≥ 0） |
| `--meta / --single` | | 自动判断 | 强制指定 Pyrodigal 的宏基因组模式或单菌模式。自动模式下，输入小于 100 kb，或者 contig 多于 200 条且都短于 20 kb 时，使用 meta 模式。`--single` 遇到序列太短（< 20 kb）无法训练时，会自动改用 meta 模式并给出警告。输入里已有 CDS 或蛋白时，这个参数不起作用 |
| `--identity` | `-i` | `0.3` | Clinker 同源比对条带的最小一致度阈值（0.0 ~ 1.0） |
| `--clinker-opts` | | `None` | 透传给 Clinker 的参数（例如 `"-na -ufo"`）。`-mo matrix.tsv` 这类相对路径的输出文件会写到输出目录里 |
| `--clinker-help` | | | 在终端打印 Clinker 常用参数指南并退出 |
| `--static-plot` | | `False` | 额外在 `figures/` 里生成静态 SVG 矢量图 |
| `--no-plot` | | `False` | 跳过 Clinker 绘图，只输出表格、序列和 GenBank 切片 |
| `--log` | | `<outdir>/magcluster.log` | 把运行日志写到指定文件 |
| `--help` | `-h` | | 显示帮助并退出 |

全局选项：`magcluster --version`（`-v`）和 `magcluster --help`（`-h`）。

基因组 ID 用于输出文件名、`genome` 列和预测基因的 gene ID。它取自输入文件名，去掉 `.gz` 和基因组扩展名，例如 `GCF_000513295.1.fna` → `GCF_000513295.1`，`strainX.v1.fna` → `strainX.v1`。空格和特殊字符会替换成 `_`（`my strain.fna` → `my_strain`）。两个输入得到同一个 ID 时（例如 NCBI Datasets 下载的 `genomic.gbff`），后面的会加后缀（`genomic_2`），并给出警告。同一个文件重复输入只处理一次。

退出码：成功为 `0`。有输入文件读不了、Clinker 失败，或者输出目录不可用时为 `1`，其他基因组的结果照常写出。命令行参数不合法时为 `2`。

### `magcluster db`

查看内置参考数据库信息、10 个核心判定家族，以及全部 92 个 HMM 模型的参考序列数与中位长度：

```bash
magcluster db
```

---

## 输出结果解读

运行 `magcluster run <inputs> -o results/` 后，结果文件夹结构如下：

```text
results/
├── {genome}_mgc_genes.tsv              # 磁小体基因注释总表（每个样本一个文件）
├── {genome}_magnetosome_proteins.faa   # 提取的磁小体蛋白序列（每个样本一个文件）
├── MGC_align.html                      # 交互式基因簇可视化图（用浏览器打开；--no-plot 时不生成）
├── magcluster.log                      # 运行日志（见下文）
├── clusters_gbk/                       # 切片好的微型 GenBank 文件
│   └── {cluster_id}.gbk
└── figures/                            # 仅在使用 --static-plot 时生成
    └── {cluster_id}.svg                # 静态基因簇图（pyGenomeViz）
```

没有任何磁小体命中的基因组不会生成表格和 FASTA，但会记录在 `magcluster.log` 里。

### 1. 基因注释表格（`{genome}_mgc_genes.tsv`）

| 列名 | 说明 |
| :--- | :--- |
| `cluster_id` | 所属基因簇 ID（如 `MSR-1_genome_MGC`）；若未成簇则记为 `unclustered` |
| `genome` | 基因组样本名称 |
| `gene_id` | 基因名称或 Locus tag |
| `contig` | 所在的 Contig / 染色体 ID |
| `start`, `end` | 基因在染色体上的起止坐标（1-based） |
| `strand` | 正负链方向（`+` 或 `-`） |
| `profile` | 具体匹配到的 HMM 模型（如 `MamE_full`、`MamE_Cter`、`MamQ_greigite`） |
| `parent_family` | 家族名称（如 `MamE`、`MamQ`） |
| `role` | 基因功能角色：`core`（核心）、`accessory`（辅助）或 `iron_context`（铁转运） |
| `score`, `evalue` | HMM 检索得分与 E 值 |

### 2. 磁小体蛋白序列（`{genome}_magnetosome_proteins.faa`）

提取出的磁小体蛋白氨基酸序列，FASTA 表头附带详细的比对和簇信息：
```text
>MGMSRV2_RS11765 [MSR-1_genome] cluster=MSR-1_genome_MGC family=MamE profile=MamE_full score=1000.4 evalue=1.97e-301
MSERRPVGRTL...
```

### 3. 交互式基因簇对齐图（`MGC_align.html`）

直接用浏览器打开：

- **基因标签**：基因箭头会标出具体的结构域或变体标签（如 `MamE_Cter`、`MamE_Nter`、`MamQ_greigite`）。
- **家族配色板**：
  - 红色：`mam` 家族
  - 黄色：`mms` 家族
  - 绿色：`mad` 家族
  - 蓝色：`man` 家族
  - 暖黄：`feo` 铁转运蛋白
  - 紫色：`mae` 环境新型蛋白
  - 青绿：`maq` 新型蛋白
  - 洋红：`mca` 新型蛋白
  - 灰色：侧翼其他非磁小体基因
- **网页交互**：可在页面上拖动基因、修改颜色、调节对比度，并直接一键导出用于文章发表的高清 **SVG** 矢量图。

### 4. 运行日志（`magcluster.log`）

纯文本日志。每次运行写入同一个目录时会覆盖（例如使用 `--force` 时）。内容包括：

- 开始和结束时间，MagCluster、Python 和各依赖的版本；
- 完整的命令行和所有实际生效的参数；
- 每个输入文件对应的基因组 ID；
- 每个基因组的处理进度、耗时和结果（基因数、命中数、MGC 坐标、各亚簇及其家族、矿物类型）；
- 所有警告和错误，包括被跳过的文件。

用 `--log FILE` 可以写到其他位置。

### 5. 静态图（`figures/{cluster_id}.svg`）

使用 `--static-plot` 时生成，每个 MGC 一个 SVG，每个亚簇一条轨道，用 pyGenomeViz 绘制。

---

## 实用小建议

### 基因组较多时加快 Clinker 绘图
当比对多个较大基因簇时，基因两两全比对可能会耗费较多时间。可以通过 `--clinker-opts` 透传参数提速：

```bash
# 跳过两两相似度比对条带，只绘制基因簇图（速度极快）：
magcluster run *.gbff -o results/ --clinker-opts "-na"

# 保持文件输入原本的上下排列顺序（不自动按相似度重排）：
magcluster run *.gbff -o results/ --clinker-opts "-ufo"

# 组合使用：
magcluster run *.gbff -o results/ --clinker-opts "-na -ufo"
```

如果想查看 Clinker 支持的所有可视化参数，直接输入：
```bash
magcluster run --clinker-help
```

---

## MGC 的判定方法

1. 同一条 contig 上的磁小体 HMM 命中，如果相邻两个命中之间最多隔 `--max-gap`（默认 5）个非磁小体基因、并且距离不超过 25 kb，就归为一组。Feo（铁转运）命中不能开启或延长一个组。
2. 一个组里至少要有 `--threshold`（默认 3）个不同的磁小体家族才保留。**任意非 Feo 家族都算**，既包括核心 Mam 家族，也包括 Mad、Mms、Man、Mae、Maq、Mca 等辅助家族；Feo 家族不算。
3. 保留下来的组两端各延伸 `--flank`（默认 3）个基因。同一个基因组里所有保留的组合并报告为一个 MGC（`{genome}_MGC`），每个组是一个亚簇（`loc1`、`loc2`……）。
4. 在完整的环状复制子上（GenBank 标注为 `circular`，或者 FASTA 表头含有 "complete genome"、"complete sequence" 或 "circular"），如果 MGC 的亚簇分布在原点两侧，报告的坐标会是 `start > end`，并标注 "(spans origin)"。

`magcluster db` 中列出的 10 个核心判定家族如下：

$$\text{MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ}$$

这 10 个家族仅供参考，判定 MGC 时不要求必须包含它们。

内置 `rules.yaml` 里的 `clustering_params` 和每个模型的 `evalue_cutoff` 只是说明性字段，程序不会读取，请使用上面的命令行参数。

---

## 引用文献

如果您在研究中使用了 MagCluster，请引用：

> Runjia Ji, Wensi Zhang, Yongxin Pan, and Wei Lin. (2022).  
> **MagCluster: a Tool for Identification, Annotation, and Visualization of Magnetosome Gene Clusters.**  
> *Microbiology Resource Announcements*, 11:e01031-21.  
> [https://doi.org/10.1128/mra.01031-21](https://doi.org/10.1128/mra.01031-21)

MagCluster 底层核心工具引用：
- **Pyrodigal**: [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **PyHMMER**: [doi:10.1093/bioinformatics/btad211](https://doi.org/10.1093/bioinformatics/btad211)
- **Clinker**: [doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007)

---

## 联系与反馈

如有任何问题或建议，欢迎联系：
- **GitHub Issues**：[提交 Issue](https://github.com/RunJiaJi/magcluster/issues)
- **邮箱**：Runjia Ji (`jirunjia@gmail.com`) 或 Wei Lin (`weilin@mail.iggcas.ac.cn`)
