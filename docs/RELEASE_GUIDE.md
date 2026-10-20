# MagCluster v1.0.1 打包与发布完整指南 (PyPI & Bioconda)

本指南旨在为维护者提供清晰、可落地的 **MagCluster v1.0.1** 发布操作步骤。

---

## 核心发布顺序

由于 Bioconda 的官方 Recipe 默认直接从 PyPI 源下载 `.tar.gz` 源码包并校验 SHA256，因此标准发布顺序为：

$$\text{本地打包与测试} \longrightarrow \text{发布至 PyPI} \longrightarrow \text{向 Bioconda 提交 PR}$$

---

## 第一部分：PyPI 发布

### 1. 产物构建与本地校验

仓库已配置好现代化 `pyproject.toml`，内置数据（92 个 HMM 模型、rules.yaml、profiles_summary.tsv）已正确包含。

在 Conda 环境中执行一键打包脚本：
```bash
# 激活环境
conda activate magcluster2

# 运行一键构建与检查脚本
./scripts/publish_pypi.sh
```

构建完成后，`dist/` 目录下将生成两个核心文件：
- `dist/magcluster-1.0.1.tar.gz` (源码分发包 sdist)
- `dist/magcluster-1.0.1-py3-none-any.whl` (二进制轮子包 wheel)

且 `twine check dist/*` 结果必须为 `PASSED`。

---

### 2. 发布方式 A：本地命令行直接上传（最快捷）

1. **准备 PyPI API Token**：
   - 登录 [pypi.org](https://pypi.org/)；
   - 进入 **Account Settings** $\to$ **API tokens** $\to$ **Add API token**；
   - 范围选择 `Project: magcluster`（或所有项目）；
   - 复制生成的 Token（形如 `pypi-AgEIcHlwaS5vcmc...`）。

2. **执行上传**：
   ```bash
   # 上传正式 PyPI
   twine upload dist/*
   ```
   - **Username**: 输入 `__token__`
   - **Password**: 粘贴刚才复制的 `pypi-...` Token

3. **验证安装**：
   发布成功后约 1~2 分钟，即可通过 pip 测试安装：
   ```bash
   pip install --upgrade magcluster
   magcluster --help
   ```

---

### 3. 发布方式 B：通过 GitHub Actions 自动化发布

仓库中已配置好 `.github/workflows/publish-pypi.yml`：
1. 打开 GitHub 仓库设置：`Settings` $\to$ `Secrets and variables` $\to$ `Actions`；
2. 新增 Repository Secret：
   - Name: `PYPI_API_TOKEN`
   - Secret: 粘贴你的 PyPI API Token（`pypi-...`）；
3. 在 GitHub 仓库页面发布新 Release（例如 Tag 为 `v1.0.1`），GitHub Actions 将自动触发打包、校验并上传至 PyPI。

---

## 第二部分：Bioconda 发布

MagCluster 在早期（v0.2.x）已收录于 Bioconda。本次 v1.0.0 彻底移除了 `prokka` 和 `blast`，替换为原生 `pyrodigal` 和 `pyhmmer`。

我们已在仓库的 `recipes/magcluster/meta.yaml` 中准备好了最新的 Recipe 文件。

### 1. Recipe 关键信息核对

- **Package Name**: `magcluster`
- **Version**: `1.0.1`
- **Source URL**: `https://pypi.org/packages/source/m/magcluster/magcluster-1.0.1.tar.gz`
- **SHA256**: `63e34f2c6257f5e673152a41d1246d94e9f41cbafc2129d727375031be0f7439`
- **核心依赖变化**：
  - ❌ 彻底移除：`prokka`、`blast`、`tbl2asn`
  - ✅ 运行依赖（1.0.1）：`biopython >=1.80`、`pyrodigal >=3.0.0`、`pyhmmer >=0.10.1`、`typer >=0.16.0`、`rich >=13.0.0`、`pyyaml >=6.0`、`pygenomeviz >=1.0.0`、`clinker-py >=0.0.30`（1.0.1 起不再依赖 `pandas` 和 `pyfamsa`）

---

### 2. 提交更新至 Bioconda

1. **Fork 与分支**：
   - 访问 [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)；
   - 点击右上角 **Fork** 到你的个人 GitHub 账号；
   - 在个人 Fork 仓库中基于 `master` 创建新分支（例如 `magcluster-1.0.1`）。

2. **更新 meta.yaml**：
   - 打开 `recipes/magcluster/meta.yaml` 文件；
   - 将本项目中 [recipes/magcluster/meta.yaml](file:///Users/sybil/Downloads/win.DESKTOP-DUDCGV3/Desktop/博士后/11_magcluster_更新/MagCluster/recipes/magcluster/meta.yaml) 的全部内容复制替换进去；
   - 提交 Commit：`magcluster: update to 1.0.1`。

3. **创建 Pull Request**：
   - 从你的分支向 `bioconda/bioconda-recipes:master` 发起 PR；
   - PR 模板中勾选相应选项。Bioconda 的 CI 机器人（Azure Pipelines）会自动拉取 PyPI 包、构建 conda 包并运行测试（`magcluster --help` 等）；
   - CI 全部绿灯后，在 PR 中评论：
     ```text
     @bioconda-bot please merge
     ```
   - 机器人审核通过后将自动合并，合并后约 1~2 小时内，全网即可通过 conda 安装：
     ```bash
     conda install -c bioconda -c conda-forge magcluster=1.0.1
     ```

---

## 常见问题处理

1. **若重新打包导致 SHA256 改变**：
   - 如果在正式上传前修改了代码并重新运行了 `python -m build`，请务必运行：
     ```bash
     shasum -a 256 dist/magcluster-1.0.1.tar.gz
     ```
   - 并同步更新 `recipes/magcluster/meta.yaml` 中的 `sha256` 字段。

2. **Bioconda 依赖名提示**：
   - 注意：在 Conda 体系中，`clinker` 的包名为 `clinker-py`，Recipe 中已正确指定为 `clinker-py >=0.0.30`。
