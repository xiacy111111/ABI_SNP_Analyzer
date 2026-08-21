# ABI-SNP Analyzer v1.0.0 用户操作手册

---

## 目录

1. [软件概述](#1-软件概述)
2. [系统要求](#2-系统要求)
3. [安装与启动](#3-安装与启动)
4. [界面说明](#4-界面说明)
5. [参数设置说明](#5-参数设置说明)
6. [操作流程](#6-操作流程)
7. [输出文件说明](#7-输出文件说明)
8. [种群遗传分析](#8-种群遗传分析)
9. [数据导出格式](#9-数据导出格式)
10. [操作示例](#10-操作示例)
11. [常见问题](#11-常见问题)

---

## 1. 软件概述

ABI-SNP Analyzer 是一款免费、开源的桌面应用程序，用于从 ABI 格式 Sanger 测序色谱图文件（.ab1）自动进行 SNP 基因分型和种群遗传分析。软件提供端到端的图形化工作流程，将 ABIF 文件解析、杂合子检测、多序列比对、SNP 识别、二倍体基因型矩阵生成、参考基因组坐标锚定以及种群遗传分析整合于单一界面中。

**主要功能模块：**

| 模块 | 功能 |
|---|---|
| ABIF 文件解析 | 读取 .ab1 二进制格式，提取碱基序列、质量分数和四通道荧光信号 |
| 杂合子检测 | 基于多通道信号比率分析，自动识别杂合位点并生成 IUPAC 模糊碱基码 |
| 序列比对与 SNP 识别 | MUSCLE 多序列比对、相似性过滤、SNP 位点发现与统计 |
| 参考基因组映射 | BLASTN 将一致性序列比对到参考基因组，获取染色体坐标 |
| 基因型矩阵 | 生成样本×位点二倍体基因型矩阵（Excel 格式） |
| 种群遗传分析 | 群体管理、遗传多样性统计、Hardy-Weinberg 平衡检验、F-统计量、AMOVA |
| 数据导出 | 支持 Genepop、STRUCTURE、FSTAT、Arlequin、POPGENE32 五种格式 |

---

## 2. 系统要求

| 项目 | 最低配置 | 推荐配置 |
|---|---|---|
| 操作系统 | Windows 7（64位） | Windows 10/11（64位） |
| 内存 | 4 GB | 8 GB 及以上 |
| 硬盘空间 | 500 MB | 1 GB 及以上（不含参考基因组） |
| 外部工具 | MUSCLE v5 | BLAST+ 2.12+（可选，用于基因组坐标映射） |

**外部依赖工具（需单独下载）：**
- **MUSCLE v5**：多序列比对工具，下载地址 https://www.drive5.com/muscle/
- **BLAST+ 2.12+**：参考基因组比对工具（可选），下载地址 https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

> 提示：外部工具只需将 `muscle.exe`（和可选的 `blastn.exe`）放置在软件同一目录下，或通过界面指定路径。

---

## 3. 安装与启动

### 3.1 独立可执行文件（推荐）

1. 从 [GitHub Releases](https://github.com/xiacy111111/ABI_SNP_Analyzer/releases) 下载 `ABI_SNP_Analyzer_v1.0.0.exe`。
2. 将 `muscle.exe`（以及可选的 `blastn.exe`）放入与 EXE 相同的文件夹。
3. 双击 `ABI_SNP_Analyzer_v1.0.0.exe` 启动软件。

### 3.2 从源代码运行

```bash
git clone https://github.com/xiacy111111/ABI_SNP_Analyzer.git
cd abi-snp-analyzer && git checkout v1.0.0
python -m venv .venv && .venv\Scripts\activate
pip install -r requirements.txt
python src/abi_to_snp_analyzer_v1.0.0.py
```

---

## 4. 界面说明

软件界面分为四个功能区：

```
┌─────────────────────────────────────────────────────────┐
│  ① 文件选择区（File Selection）                         │
│  ├─ ABI Files: 选择 .ab1 测序文件或文件夹               │
│  ├─ MUSCLE Path: 指定 muscle.exe 路径                   │
│  ├─ Reference Genome: 参考基因组 FASTA 文件（可选）      │
│  ├─ BLASTN Path: 指定 blastn.exe 路径（可选）            │
│  └─ Output Directory: 选择结果输出文件夹                 │
├─────────────────────────────────────────────────────────┤
│  ② 参数设置区（Filter Threshold Settings）              │
│  ├─ Sequence Similarity / MAF Threshold                 │
│  ├─ Min Observed Het / Signal Ratio                     │
│  └─ Consensus Quality / Min Base Quality                │
├─────────────────────────────────────────────────────────┤
│  ③ 运行控制区                                           │
│  ├─ [Start Analysis] [Stop] [Load Existing Data]       │
│  └─ [Population Manager] [Statistical Analysis]        │
│     [Export SNP Formats]                               │
├─────────────────────────────────────────────────────────┤
│  ④ 运行日志区（Run Log）                                │
│     实时显示分析进度、参数信息和错误提示                  │
└─────────────────────────────────────────────────────────┘
```

---

## 5. 参数设置说明

| 参数 | 默认值 | 范围 | 说明 |
|---|---|---|---|
| **Sequence Similarity** | 0.90 | 0–1 | 序列相似性阈值。低于该值的序列将被排除。建议同扩增子样品设为 0.85–0.95 |
| **MAF Threshold** | 0.10 | 0–1 | 最小等位基因频率阈值。仅 MAF ≥ 该值的 SNP 位点被保留 |
| **Min Observed Het** | 0.10 | 0–1 | 最小观测杂合度。低于该值的位点被过滤。设为 0 表示不过滤 |
| **Signal Ratio** | 0.30 | 0–1 | 杂合信号比率阈值（次级峰/主峰）。值越低越敏感（检测更多杂合位点），越高越特异 |
| **Consensus Quality** | 0.20 | 0–1 | 一致性序列质量阈值。低于该阈值的位置将被标记为 N |
| **Min Base Quality** | 5 | 0–100 | Phred 碱基质量阈值。低于该值的碱基不参与基因分型 |

---

## 6. 操作流程

### 6.1 基础流程：SNP 基因分型

**步骤 1 — 准备输入文件**
将所有待分析的 .ab1 文件放入同一个文件夹中。确保文件名唯一，且所有文件来自同一目标扩增区域。

**步骤 2 — 配置文件路径**
- 点击 **Select Files/Folder**，选择包含 .ab1 文件的文件夹（或直接选择多个 .ab1 文件）
- 点击 **MUSCLE Path** 旁的 **Browse**，选择 `muscle.exe` 文件
- （可选）点击 **Reference Genome** 旁的 **Browse**，选择参考基因组 FASTA 文件
- （可选，若选了参考基因组则必须）点击 **BLASTN Path** 旁的 **Browse**，选择 `blastn.exe`
- 点击 **Output Directory** 旁的 **Select**，选择结果输出目录

**步骤 3 — 设置参数**
根据实验需求调整六个参数（详见[第 5 节](#5-参数设置说明)）。默认参数适用于大多数常规分析场景。

**步骤 4 — 运行分析**
点击 **Start Analysis** 按钮。软件将自动执行以下子步骤：

1. **ABIF 解析与杂合子检测**：读取每个 .ab1 文件，提取色谱信号，识别杂合位点，生成 IUPAC 编码的 FASTA 序列
2. **序列相似性过滤**：去除与群体差异过大的序列
3. **MUSCLE 多序列比对**：对所有保留序列进行比对
4. **一致性序列与峰图生成**：生成共有序列，为有杂合位点的样品绘制色谱峰图和杂合子报告
5. **SNP 位点检测与过滤**：基于 MAF 和杂合度阈值筛选 SNP
6. **（可选）BLASTN 基因组坐标映射**：将一致性序列比对到参考基因组
7. **基因型矩阵生成**：输出样本×位点的二倍体基因型 Excel 表格

**步骤 5 — 查看结果**
分析完成后，软件弹出"Analysis completed"提示框。在输出目录中查看生成的文件（详见[第 7 节](#7-输出文件说明)）。

### 6.2 续接分析：加载已有数据

如果之前已完成基因分型，可直接加载已有结果继续进行种群遗传分析：

1. 点击 **Load Existing Data** 按钮
2. 选择之前生成的 `genotype_matrix.xlsx` 文件
3. 数据加载后，**Population Manager** 和 **Export SNP Formats** 按钮将被激活

---

## 7. 输出文件说明

分析完成后，输出目录中会生成以下文件：

| 文件名 | 格式 | 说明 |
|---|---|---|
| `{样品名}_sequence_with_iupac.fasta` | FASTA | 每个样品的 IUPAC 编码序列 |
| `{样品名}_trace_plot_combined.png` | PNG | 色谱信号峰图，标注杂合位点（仅含杂合子的样品） |
| `{样品名}_heterozygote_report.xlsx` | Excel | 每个样品的杂合位点详情（位置、信号比率等） |
| `combined_all_forward_sequences_with_iupac.fasta` | FASTA | 所有样品 IUPAC 序列合并文件 |
| `alignment.aln` | 文本 | MUSCLE 多序列比对结果 |
| `consensus_sequence.fasta` | FASTA | 一致性序列 |
| `consensus_vs_reference.aln` | 文本 | 一致性序列与参考基因组 BLASTN 比对结果（可选） |
| `genotype_matrix.xlsx` | Excel | 样本×SNP 位点二倍体基因型矩阵（sheet: Genotype_Matrix） |
| `genotype_matrix_snp_analysis.xlsx` | Excel | 每个 SNP 位点的统计信息：MAF、杂合度、基因组坐标等（sheet: SNP_Analysis_Results） |
| `analysis_metadata.json` | JSON | 运行参数、输入文件、依赖版本等分析重现性元数据 |

**基因型编码说明：**
- 纯合子：`AA`、`TT`、`CC`、`GG`
- 杂合子：`AG`、`CT`、`AC` 等
- 缺失：`--`

**SNP 位点命名：**
- 无参考基因组：`SNP_586`（一致性序列上的位置）
- 有参考基因组：`chr1:1234567G>A`（染色体:位置参考碱基>替代碱基）

---

## 8. 种群遗传分析

### 8.1 群体管理

1. 基因分型完成后，点击 **Population Manager** 按钮
2. 在弹出的对话框中：
   - 左侧显示所有样品列表
   - 点击 **Add Population** 创建新群体，输入群体名称
   - 从左侧选中样品，点击 **>** 将样品分配到右侧群体中
   - 每个样品只能属于一个群体
   - 至少需定义两个群体才能进行后续分析

### 8.2 统计分析

群体分配完成后，点击 **Statistical Analysis**，软件将自动计算以下统计量：

| 分析项 | 输出内容 |
|---|---|
| 数据概览 | 各群体样本数和位点数 |
| 遗传多样性 | 观测杂合度 (Ho)、期望杂合度 (He)、多态信息含量 (PIC) |
| Hardy-Weinberg 平衡检验 | 每个位点在每个群体中的 HWE 精确检验 p 值 [Guo & Thompson 1992; Wigginton et al. 2005] |
| WC84 F-统计量 | Fis（群体内近交系数）、Fst（群体间分化系数）、Fit（总近交系数），含 999 次置换检验 [Weir & Cockerham 1984] |
| 零等位基因频率 | Brookfield 估计量 [Brookfield 1996] |
| AMOVA | 分子方差分析，含 Phi-统计量和显著性 p 值 [Excoffier et al. 1992] |
| 基因型频率卡方检验 | 群体间基因型频率差异检验，Bonferroni 校正 p 值 |
| 交互式可视化 | 每个位点的基因型频率柱状图（支持 300 DPI PNG/PDF 保存） |
| 结果汇总 | 所有统计结果汇总为 Excel 表格，包含原始 p 值和多重检验校正值 |

---

## 9. 数据导出格式

软件支持将 SNP 基因型矩阵导出为五种常用种群遗传学软件格式：

| 格式 | 扩展名 | 编码方式 | 兼容软件 |
|---|---|---|---|
| **Genepop** | .txt | 2位整数（A=01, T=02, C=03, G=04）；0000=缺失 | Genepop, Fstat |
| **STRUCTURE** | .txt | 每SNP占2列（等位基因1, 等位基因2）；-9=缺失 | STRUCTURE |
| **FSTAT** | .dat | 整数等位基因编码，含位点名称头信息 | FSTAT |
| **Arlequin** | .arp | Arlequin 3.5 项目文件，含 [[Loci]] 和 [Structure] 块 | Arlequin |
| **POPGENE32** | .pop | 两位字母基因型代码；`..`=缺失 | POPGENE 1.32 |

**导出步骤：**
1. 基因分型完成或加载已有数据后，点击 **Export SNP Formats**
2. 在弹出窗口中选择目标格式
3. 勾选"Include population assignment info"可包含群体分组信息
4. 点击 **Export**，选择保存路径

---

## 10. 操作示例

以下以红罗非鱼（*Oreochromis* spp.）耐盐性 QTL 区域验证研究为例，演示完整操作流程。

### 10.1 数据准备
- 88 个 .ab1 Sanger 测序文件（44 个耐盐个体 + 44 个盐敏感个体）
- MUSCLE v5 (`muscle.exe`)
- 尼罗罗非鱼参考基因组 FASTA 文件（约 1 GB）
- NCBI BLAST+ (`blastn.exe`)

### 10.2 参数设置
| 参数 | 设定值 |
|---|---|
| Sequence Similarity | 0.90 |
| MAF Threshold | 0.10 |
| Min Observed Het | 0.10 |
| Signal Ratio | 0.30 |
| Consensus Quality | 0.20 |
| Min Base Quality | 5 |

### 10.3 运行
1. 选择包含 88 个 .ab1 文件的文件夹
2. 配置 MUSCLE、参考基因组和 BLASTN 路径
3. 选择输出目录
4. 设置上述参数
5. 点击 **Start Analysis**，等待约 8 分钟

### 10.4 分析结果
- 一致性序列长度：861 bp
- 检测到 SNP 位点：10 个
- 基因型数据：880 个（88 样本 × 10 位点）
- 基因分型准确率：99.8%（878/880 与人工核实一致）

### 10.5 种群遗传分析
1. 点击 **Population Manager**，创建两个群体：
   - "Salt-Tolerant"（44 个耐盐个体）
   - "Salt-Sensitive"（44 个盐敏感个体）
2. 点击 **Statistical Analysis**，查看群体遗传统计结果
3. 点击 **Export SNP Formats** 导出数据用于其他软件分析

---

## 11. 常见问题

**Q: 软件提示"MUSCLE file not found"怎么办？**
A: 请确认 `muscle.exe` 路径设置正确。可从 https://www.drive5.com/muscle/ 下载 MUSCLE v5。

**Q: 参考基因组映射需要多长时间？**
A: 取决于参考基因组大小。尼罗罗非鱼基因组（约 1 GB）参考下，BLASTN 比对通常需要数秒至数十秒。

**Q: 批量处理 400 个样品需要多长时间？**
A: 在搭载 BLASTN 参考基因组映射的情况下，约需 37 分钟（i7-10700F, 16 GB RAM 测试结果）。

**Q: 杂合子检测的准确性如何？**
A: 默认信号比率阈值 0.30 下，与人工核实的一致性约为 99.8%。可通过调整 Signal Ratio 参数改变检测灵敏度。

**Q: 软件是否支持多倍体物种？**
A: 当前版本仅支持二倍体物种的基因分型和种群遗传分析。

**Q: 输出文件在哪里？**
A: 所有输出文件位于启动分析前选择的 Output Directory 中。可在输出目录中按修改时间排序查找最新生成的文件。

**Q: 如何停止正在运行的分析？**
A: 点击 **Stop** 按钮即可安全中止当前分析。

**Q: 可以在不同电脑上共享分析设置吗？**
A: 每次运行都会生成 `analysis_metadata.json` 文件，记录所有参数和依赖版本，可用于分析重现。

---

## 联系与支持

- **GitHub**: https://github.com/xiacy111111/ABI_SNP_Analyzer
- **问题反馈**: https://github.com/xiacy111111/ABI_SNP_Analyzer/issues
- **邮箱**: 25061300031@stu.xidian.edu.cn

---

*ABI-SNP Analyzer v1.0.0 | MIT License | Copyright (c) 2026 Chengyu Xia, Wenhui Chu, Pengbo Ning*
