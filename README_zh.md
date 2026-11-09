# MagCluster <img src='https://raw.githubusercontent.com/RunJiaJi/MagCluster/main/docs/_static/logo_new.svg' align='right' width=220 alt="MagCluster Logo" />

[![Python Version](https://img.shields.io/badge/python-3.10%2B-blue?style=flat-square&logo=python)](https://www.python.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-green.svg?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/version-1.1.1-orange?style=flat-square)](https://github.com/RunJiaJi/MagCluster)

[English](https://github.com/RunJiaJi/MagCluster/blob/main/README.md) | [中文说明](https://github.com/RunJiaJi/MagCluster/blob/main/README_zh.md)

**MagCluster** 是一个用于在细菌基因组中**识别、注释和可视化磁小体基因簇（MGCs）**的生信工具。

---

## v1.1 新增了哪些功能？

- **独立的 `magcluster clinker` 绘图子命令**：不再需要为了调配色、调阈值或换顺序而重跑基因预测和 HMM 比对，直接基于已生成的 `clusters_gbk/` 或 `.gbk` 文件即可秒级重绘，自动保留 8 大磁小体蛋白家族配色和结构域标签；支持用 `-l/--list` 指定基因簇排列顺序（例如按进化树顺序），也支持直接传入外部 GenBank 文件进行磁小体基因自动注释与可视化。
- **大规模基因组自动分批绘图（`--batch-size`）**：当检测到的 MGC 数量较多时，默认按输入顺序每 **50 个基因簇** 自动拆分为多个交互式网页（`MGC_align_batch01.html`、`MGC_align_batch02.html`……），避免上百个基因簇塞进单张 HTML 导致的渲染与浏览器打开卡顿；设为 `0` 可关闭分批。
- **全队列汇总表（`mgc_summary.tsv`）**：输出目录新增样本级总表，每个输入基因组占一行，汇总成簇状态、预测基因数、HMM 命中数、亚簇数、各亚簇精确坐标、矿物类型、边界/跨原点标签和 `.gbk` 相对路径，方便大批量筛选统计。
- **超长基因组名零告警兼容**：针对携带完整物种名的长文件名（如 `GCA_002753255.1_Candidatus_...`），保持 `{cluster_id}.gbk` 文件名及 Clinker 轨道标题完整不变，仅将 `.gbk` 内部 `LOCUS` / `ACCESSION` / `VERSION` 头信息智能缩短为标准登录号，彻底消除 Biopython 长度告警。

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

### 6. 对已有基因簇重新绘图（无需重跑预测与比对）
已经生成了 `clusters_gbk/`（或外部 `.gbk` 文件）时，用 `magcluster clinker` 即可快速调整配色、阈值、排序或分批：
```bash
magcluster clinker results/ -o results/MGC_replot.html --clinker-opts "-na -ufo"
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
| `--force` | `-f` | `False` | 允许使用非空的输出目录。会先删除该目录里上一次 MagCluster 运行的输出（`mgc_summary.tsv`、样本表格、FASTA、`MGC_align*.html`、`clusters_gbk/`、`figures/`、`magcluster.log`），其他文件保留 |
| `--threads` | `-t` | `0` | CPU 线程数（`0` 表示自动检测当前可用/分配的 CPU 核心数，兼容 SLURM 等集群 `os.sched_getaffinity` 配额；超过可用核心数时自动降到可用核心数） |
| `--evalue` | `-e` | `1e-5` | HMM 比对的 E-value 阈值（≥ 0） |
| `--score` | `-s` | 各模型的 `trusted_cutoff`（25.0） | HMM 比对的最小 bit score。不指定时使用 `rules.yaml` 里每个模型的 `trusted_cutoff`（92 个内置模型都是 25.0）。指定后对所有模型统一使用这个值，也可以低于 25 |
| `--max-gap` | | `5` | 基因簇内两个磁小体基因之间允许插入的最大非磁小体基因数（≥ 0）。另有固定的 25 kb 距离上限 |
| `--threshold` | `-th`、`--min-families` | `3` | 判定一个 MGC 所需的最少磁小体家族数（任意非 Feo 家族都算，≥ 1） |
| `--flank` | `-un`、`--unmag-num` | `3` | 基因簇两端各保留的非磁小体基因数（≥ 0） |
| `--meta / --single` | | 自动判断 | 强制指定 Pyrodigal 的宏基因组模式或单菌模式。自动模式下，输入小于 100 kb，或者 contig 多于 200 条且都短于 20 kb 时，使用 meta 模式。`--single` 遇到序列太短（< 20 kb）无法训练时，会自动改用 meta 模式并给出警告。输入里已有 CDS 或蛋白时，这个参数不起作用 |
| `--batch-size` | `-b` | `50` | 每个交互式 Clinker HTML 网页最多容纳的基因簇（MGC）数量，按输入顺序自动切分（`1` 个批次输出 `MGC_align.html`，多个批次输出 `MGC_align_batch01.html`、`MGC_align_batch02.html`……；设为 `0` 关闭分批） |
| `--identity` | `-i` | `0.3` | Clinker 同源比对条带的最小一致度阈值（0.0 ~ 1.0） |
| `--clinker-opts` | | `None` | 透传给 Clinker 的参数（例如 `"-na -ufo"`）。`-mo matrix.tsv` 这类相对路径的输出文件会写到输出目录里（多分批时自动添加 `_batchNN` 后缀） |
| `--clinker-help` | | | 在终端打印 Clinker 常用参数指南并退出 |
| `--static-plot` | | `False` | 额外在 `figures/` 里生成静态 SVG 矢量图 |
| `--no-plot` | | `False` | 跳过 Clinker 绘图，只输出表格、序列和 GenBank 切片 |
| `--log` | | `<outdir>/magcluster.log` | 把运行日志写到指定文件 |
| `--help` | `-h` | | 显示帮助并退出 |

全局选项：`magcluster --version`（`-v`）和 `magcluster --help`（`-h`）。

基因组 ID 用于输出文件名和 `mgc_summary.tsv` 汇总表。它取自输入文件名，去掉 `.gz` 和基因组扩展名，例如 `GCF_000513295.1.fna` → `GCF_000513295.1`，`strainX.v1.fna` → `strainX.v1`。空格和特殊字符会替换成 `_`（`my strain.fna` → `my_strain`）。两个输入得到同一个 ID 时（例如 NCBI Datasets 下载的 `genomic.gbff`），后面的会加后缀（`genomic_2`），并给出警告。同一个文件重复输入只处理一次。

退出码：成功为 `0`。有输入文件读不了、Clinker 失败，或者输出目录不可用时为 `1`，其他基因组的结果照常写出。命令行参数不合法时为 `2`。

### `magcluster clinker`

```bash
magcluster clinker [选项] [输入文件/文件夹...]
```

单独运行 Clinker 交互式基因簇可视化，并自动着色 MagCluster 的 8 大磁小体蛋白家族与结构域标签。输入支持：
- 单个或多个 `.gbk` / `.gbff` 文件（支持 `.gz` 压缩）；
- 包含 `.gbk` 文件的文件夹（例如 `results/clusters_gbk/`）；
- `magcluster run` 的输出根目录（会自动识别并读取其中的 `clusters_gbk/` 子目录）。

由 MagCluster 切片生成的 `.gbk` 文件已自带家族与模型注释，绘图时无需重跑 HMM 比对；若传入外部第三方 GenBank 文件（无 MagCluster `/note` 注释），程序会自动对其中的 CDS 运行内置 92 个磁小体 HMM 模型完成家族上色。

| 参数 | 别名 | 默认值 | 作用说明 |
| :--- | :--- | :--- | :--- |
| `--outdir` | `-o`、`--output` | `clinker_results` | 输出目录（写出 `MGC_align.html` 或 `MGC_align_batch01.html`……），或直接指定 `.html` 文件名（例如 `-o custom_mgc.html`） |
| `--list` | `-l` | `None` | 文本清单文件（每行一个 `.gbk` 路径、文件名或基因簇 ID；`#` 开头为注释），用于筛选子集或固定上下排列顺序（配合 `--clinker-opts "-ufo"` 即可严格按进化树顺序展示） |
| `--batch-size` | `-b` | `50` | 每个 HTML 文件最多绘制的基因簇数量（按输入/清单顺序切分；设为 `0` 表示全部绘入单张 HTML） |
| `--identity` | `-i` | `0.3` | Clinker 同源比对条带的最小序列一致度阈值（0.0 ~ 1.0） |
| `--threads` | `-t` | `0` | 传给 Clinker 的 CPU 线程数（`0` 表示自动检测当前可用/分配核心数） |
| `--clinker-opts` | | `None` | 透传给 Clinker 的附加参数（例如 `"-na -ufo"`、`"-asc"`、`"-mo matrix.tsv"`） |
| `--clinker-help` | | | 在终端打印 Clinker 常用参数指南并退出 |
| `--force` | `-f` | `False` | 允许覆盖已存在的目标 HTML 文件（不会删除输出目录中已有的 TSV、FASTA 或 `clusters_gbk/`） |
| `--log` | | `<outdir>/magcluster_clinker.log` | 把绘图运行日志写到指定文件 |
| `--help` | `-h` | | 显示帮助并退出 |

常用示例：
```bash
# 1. 直接对已有结果目录重新绘图（调高一致度阈值并加快速度）：
magcluster clinker results/ -o replot/ -i 0.5 --clinker-opts "-na"

# 2. 先用 --no-plot 跑完大规模基因组注释，再按每批 30 个基因簇分批绘图：
magcluster run genomes_dir/ -o results/ --no-plot
magcluster clinker results/ -o results/ -b 30 --force

# 3. 严格按照进化树叶节点清单顺序（tree_order.txt）输出单张对齐图：
magcluster clinker results/ -l tree_order.txt -o tree_mgc.html --batch-size 0 --clinker-opts "-ufo"
```

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
├── mgc_summary.tsv                     # 全队列汇总表（每个输入基因组一行）
├── {genome}_mgc_genes.tsv              # 磁小体基因命中表（含成簇基因与 unclustered 散落命中）
├── {genome}_magnetosome_proteins.faa   # 提取的磁小体蛋白序列（每个样本一个文件）
├── MGC_align.html                      # 交互式基因簇可视化图（超过 50 个簇时拆分为 MGC_align_batch01.html 等；--no-plot 时不生成）
├── MGC_align_batches.tsv               # 分批索引清单（仅在生成多个 HTML 批次时写出，记录每个基因簇所属批次）
├── magcluster.log                      # 运行日志（见下文）
├── clusters_gbk/                       # 切片好的微型 GenBank 文件
│   └── {cluster_id}.gbk
└── figures/                            # 仅在使用 --static-plot 时生成
    └── {cluster_id}.svg                # 静态基因簇图（pyGenomeViz）
```

没有任何磁小体命中的基因组不会生成单样本表格和 FASTA；仅有散落同源基因但未达到成簇阈值的基因组（`mgc_summary.tsv` 中 `status = no_MGC`）仍会生成 `{genome}_mgc_genes.tsv`（`cluster_id` 列记为 `unclustered`）以保留全部命中信息。所有输入基因组都会记录在 `mgc_summary.tsv` 和 `magcluster.log` 里。

### 0. 全队列汇总表（`mgc_summary.tsv`）

每个输入基因组一行，方便在大批量数据集中快速筛选成簇阳性样本、矿物类型和切片文件路径：

| 列名 | 说明 |
| :--- | :--- |
| `genome`, `input_file` | 标准化后的基因组 ID 与原始输入文件名 |
| `status` | 成簇状态：`MGC_detected`、`no_MGC`、`no_coordinates`、`no_genes` 或 `failed` |
| `n_genes`, `n_hits`, `n_clustered_hits` | 预测/读取的 CDS 总数、磁小体 HMM 命中总数、归入 MGC 的命中数 |
| `cluster_id`, `n_subclusters`, `coordinates` | MGC 编号、共线性亚簇数量，以及各亚簇精确基因组区间坐标 |
| `n_distinct_families`, `families` | MGC 内去重后的非 Feo 磁小体蛋白家族数及排序列表 |
| `mineral_type`, `flags` | 预测矿物类型（`Magnetite` / `Greigite` / `Putative/Mixed`）及质控标签（`Contig_Edge`、`Spans_Origin`） |
| `gbk_file` | 切片 GenBank 文件相对路径（`clusters_gbk/{cluster_id}.gbk`） |

### 1. 基因注释表格（`{genome}_mgc_genes.tsv`）

| 列名 | 说明 |
| :--- | :--- |
| `cluster_id` | 成簇基因记为 `MGC`；若未成簇则记为 `unclustered` |
| `gene_id` | 基因编号或 Locus tag（如 `gene_00405` 或 `MGMSRV2_RS11765`） |
| `contig` | 所在的 Contig / 染色体 ID |
| `start`, `end` | 基因在染色体上的起止坐标（1-based） |
| `strand` | 正负链方向（`+` 或 `-`） |
| `profile` | 具体匹配到的 HMM 模型（如 `MamE_full`、`MamE_Cter`、`MamQ_greigite`） |
| `parent_family` | 家族名称（如 `MamE`、`MamQ`） |
| `role` | 基因功能角色：`core`（核心）、`accessory`（辅助）或 `iron_context`（铁转运） |
| `score`, `evalue` | HMM 检索得分与 E 值 |

### 2. 磁小体蛋白序列（`{genome}_magnetosome_proteins.faa`）

提取出的磁小体蛋白氨基酸序列，FASTA 表头简洁附带比对和簇信息（不再重复冗长的基因组文件名）：
```text
>gene_00405 cluster=MGC family=MamK profile=MamK score=448.8 evalue=2.40e-135
MSEGEGQAKNRL...
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

### 基因组较多时控制 HTML 体积并加快 Clinker 绘图
默认开启蛋白两两比对时，Clinker 会把批次内所有基因簇两两之间的同源连线数据都写入单张 HTML 中，因此文件体积随同源基因簇数量大致按平方级（`O(N^2)`）增长——尤其在比较数十株共享大量保守基因的近缘菌株时更为明显。默认 `--batch-size 50` 适合常规跨类群数据集（分批时会自动生成 `MGC_align_batches.tsv` 方便检索每个基因簇所在的网页）；若对大批量近缘菌株绘图觉得网页偏大，可调小分批数量或透传 `--clinker-opts` 提速：

```bash
# 对高度同源的近缘菌株调小每批簇数（例如每张 HTML 20 个簇）：
magcluster clinker results/ -o plots/ -b 20

# 跳过两两相似度比对条带，只绘制基因簇图（体积缩小约 10 倍，速度极快）：
magcluster run *.gbff -o results/ --clinker-opts "-na"

# 保持文件输入原本的上下排列顺序（不自动按相似度重排）：
magcluster run *.gbff -o results/ --clinker-opts "-ufo"

# 组合使用：
magcluster run *.gbff -o results/ --clinker-opts "-na -ufo"
```

如果想查看 Clinker 支持的所有可视化参数，直接输入：
```bash
magcluster run --clinker-help
# 或
magcluster clinker --clinker-help
```

---

## MGC 的判定方法

1. 同一条 contig 上的磁小体 HMM 命中，如果相邻两个命中之间最多隔 `--max-gap`（默认 5）个非磁小体基因、并且距离不超过 25 kb，就归为一组。Feo（铁转运）命中不能开启或延长一个组。
2. 一个组里至少要有 `--threshold`（默认 3）个不同的磁小体家族才保留。**任意非 Feo 家族都算**，既包括核心 Mam 家族，也包括 Mad、Mms、Man、Mae、Maq、Mca 等辅助家族；Feo 家族不算。
3. 保留下来的组两端各延伸 `--flank`（默认 3）个基因。同一个基因组里所有保留的组合并报告为一个 MGC（`{genome}_MGC`），每个组是一个亚簇（`loc1`、`loc2`……）。
4. 在完整的环状复制子上（GenBank 标注为 `circular`，或者 FASTA 表头含有 "complete genome"、"complete sequence" 或 "circular"），如果 MGC 的亚簇分布在原点两侧，报告的坐标会是 `start > end`，并标注 "(spans origin)"。

`magcluster db` 中列出的 10 个核心判定家族如下：

`MamA`、`MamB`、`MamE`、`MamI`、`MamK`、`MamL`、`MamM`、`MamO`、`MamP`、`MamQ`

这 10 个家族仅供参考，判定 MGC 时不要求必须包含它们。

---

## 引用文献

如果您在研究中使用了 MagCluster，请引用：

> Runjia Ji, Wensi Zhang, Yongxin Pan, and Wei Lin. (2022).  
> **MagCluster: a Tool for Identification, Annotation, and Visualization of Magnetosome Gene Clusters.**  
> *Microbiology Resource Announcements*, 11:e01031-21.  
> [https://doi.org/10.1128/mra.01031-21](https://doi.org/10.1128/mra.01031-21)

MagCluster 底层核心工具引用：
- **Pyrodigal**: [doi:10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
- **PyHMMER**: [doi:10.1093/bioinformatics/btad214](https://doi.org/10.1093/bioinformatics/btad214)
- **Clinker**: [doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007)

---

## 联系与反馈

如有任何问题或建议，欢迎联系：
- **GitHub Issues**：[提交 Issue](https://github.com/RunJiaJi/MagCluster/issues)
- **邮箱**：Runjia Ji (`jirunjia@gmail.com`) 或 Wei Lin (`weilin@mail.iggcas.ac.cn`)
