# MagCluster <img src='docs/_static/logo_new.svg' align='right' width=220 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.0.0-orange?style=flat-square)](https://github.com/RunJiaJi/magcluster)

[English](README.md) | [中文说明](README_zh.md)

**MagCluster** 是一个用于在细菌基因组中**识别、注释和可视化磁小体基因簇（MGCs）**的生信工具。

---

## v1.0 做了哪些更新？

- **一条命令搞定**：新版用 `magcluster run` 整合了以前繁琐的多个步骤，一键完成基因预测、HMM 比对、基因簇识别和绘图。
- **不再依赖外部 Prokka**：内置原生 [Pyrodigal](https://github.com/althonos/pyrodigal) 进行快速基因预测，不需要再配置复杂的 Perl 环境和 BLAST 索引。
- **92 个精细 HMM 模型**：覆盖 8 大磁小体蛋白家族（`Mam`、`Mms`、`Mad`、`Man`、`Mae`、`Maq`、`Mca`、`Feo`），使用 [PyHMMER](https://github.com/althonos/pyhmmer) 进行快速检索。
- **支持多种输入格式**：支持带注释的 GenBank 文件（`.gbk`/`.gbff`）、未注释的核酸 Contig（`.fasta`/`.fna`/`.fa`），以及纯蛋白序列（`.faa`），同时支持 `.gz` 压缩文件。
- **每个基因组独立出结果**：每个样本单独生成 `{genome}_mgc_genes.tsv` 和 `{genome}_magnetosome_proteins.faa`，不再混在一起。

---

## 安装方法

### 1. 新建 Conda 环境（推荐）

```bash
conda create -n magcluster-1.0 -c conda-forge -c bioconda python=3.12 pip
conda activate magcluster-1.0
```

### 2. 安装 MagCluster

直接从 GitHub 安装：
```bash
pip install git+https://github.com/RunJiaJi/magcluster.git
```

或者下载源码本地安装：
```bash
git clone https://github.com/RunJiaJi/magcluster.git
cd MagCluster
pip install -e .
```

安装完成后检查版本：
```bash
magcluster --version
magcluster run --help
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

### 4. 分析蛋白质序列文件（`.faa`）
注意：如果输入的是没有基因组坐标的纯蛋白 FASTA，程序会完成磁小体蛋白注释并输出 TSV 表格和 FASTA 序列文件，不进行基因簇绘图：
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

| 参数 | 缩写 | 默认值 | 作用说明 |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o` | `magcluster_results` | 输出结果存放目录 |
| `--force` | `-f` | `False` | 如果输出目录非空，强制覆盖写入 |
| `--threads` | `-t` | `0` | CPU 线程数（`0` 表示使用全部可用核心） |
| `--evalue` | `-e` | `1e-5` | HMM 比对的 E-value 过滤阈值 |
| `--score` | `-s` | `25.0` | HMM 比对的最小 Bit score 过滤阈值 |
| `--max-gap` | | `5` | 基因簇内两个磁小体基因之间允许插入的最大非磁小体基因数量 |
| `--threshold` | `-th` | `3` | 判定为一个有效磁小体基因簇所需的最少磁小体家族数量 |
| `--flank` | `-un` | `3` | 基因簇两端额外保留的非磁小体基因数量 |
| `--meta / --single` | | 自动判断 | 强制指定 Pyrodigal 的宏基因组模式或单菌模式 |
| `--identity` | `-i` | `0.3` | Clinker 同源比对条带的最小一致度阈值（0.0 ~ 1.0） |
| `--clinker-opts` | | `None` | 透传给 Clinker 的原生参数（例如 `"-na -ufo"`） |
| `--clinker-help` | | | 在终端打印 Clinker 常见参数指南并退出 |
| `--static-plot` | | `False` | 额外生成静态 SVG 格式的基因簇矢量图 |
| `--no-plot` | | `False` | 只输出表格和序列，跳过绘图（适合服务器批量纯文本分析） |

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
├── MGC_align.html                      # 交互式基因簇可视化图（用浏览器双击打开）
└── clusters_gbk/                       # 切片好的微型 GenBank 文件
    └── {cluster_id}.gbk
```

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

## 10 个核心判定家族

判定一个簇是否为真实磁小体基因簇时，MagCluster 要求至少包含 3 个不同的磁小体家族。核心家族定义如下：

$$\text{MamA, MamB, MamE, MamI, MamK, MamL, MamM, MamO, MamP, MamQ}$$

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
