# GuideRNA_predict

#### 介绍
用于对非模式物种进行可用的guideRNA预测和特异性检测评估（基于序列特异性、脱靶位置、和Seed位点保守性）
输入目标区域（基因组位置或FASTA序列）和GFF3注释文件，自动扫描NGG/CCN PAM位点并提取20bp gRNA候选序列，通过全基因组BLAST模拟搜索脱靶位点，经两步威胁评估（是否在CDS区 + 是否有PAM），最终按特异性排名输出Top N推荐gRNA及其脱靶位点的基因注释信息。

#### 软件架构

```
grna_design_pipeline.py          # 主流程脚本（纯Python标准库，零第三方依赖）
├── load_genome()                #   加载基因组FASTA
├── parse_gff3_cds()             #   解析GFF3提取CDS区域
├── parse_gff3_features()        #   解析GFF3全部特征（用于脱靶注释）
├── parse_target_regions()       #   解析目标区域文件（支持多种格式）
├── locate_target_sequences()    #   从FASTA搜索目标序列在基因组中的位置
├── scan_pam_sites()             #   扫描NGG/CCN PAM位点提取gRNA候选
├── build_cds_data()             #   构建CDS串联序列（加速搜索）
├── analyze_specificity()        #   全基因组特异性分析（核心算法）
│   ├── 10bp seed基因组搜索 + 5'延伸（精确匹配）
│   ├── 4bp seed CDS区搜索（近匹配补充）
│   └── 两步威胁评估
│       ├── Step 1: 脱靶位点是否在CDS编码区？
│       └── Step 2: 脱靶位点附近是否有PAM（±2bp）？
├── annotate_threats()           #   威胁位点GFF3注释（基因号+CDS编号）
└── write_outputs()              #   输出TXT + JSON结果文件
```


#### 安装教程

1.  确保 Python 版本 >= 3.6：
    ```
    python --version
    ```

2.  下载脚本到本地工作目录即可使用，**无需安装任何第三方Python包**（仅使用 argparse, json, re, sys, time, collections 等标准库模块）。

3.  准备输入文件：
    - 基因组序列文件（FASTA格式，支持多contig）
    - 基因注释文件（GFF3格式，需包含CDS特征行）
    - 目标区域文件（位置文本）或目标序列文件（FASTA格式）

#### 使用说明

##### 1. 基本用法

```bash
python grna_design_pipeline.py \
    --genome <genome.fasta> \
    --gff3 <annotation.gff3> \
    --regions <targets.txt> \
    --out_dir <output_directory> \
    --out_prefix <output_prefix> \
    --top_n 10
```

##### 2. 三种目标区域指定方式

**方式A — 位置文件（`--regions`）**

每行一个区域，Tab或空格分隔，支持以下格式（可混合）：

```
# 完整格式：contig + 起止位置 + 标签
ctg000250  62541730  62541863  CDS8

# 简化格式：contig + 起止位置（标签自动生成为R1, R2, ...）
ctg000250  62541730  62541863

# 冒号格式（带标签）
ctg000250:62541730-62541863  CDS8

# 冒号格式（不带标签，标签自动取坐标值）
ctg000250:62541730-62541863
```

以 `#` 开头的行视为注释，会被忽略。
文件后缀为.txt、.tsv、.bed或无后缀均可。

**方式B — 序列文件（`--target`）**

提供目标序列的FASTA文件，脚本自动在基因组中搜索每条序列的位置（正链和反链均搜索）。FASTA header 作为区域标签：

```fasta
>CDS8
ATGGCAATTGCCTAACCCAGAGAATGGGCGTCGGGGACT...
>CDS9
TGGAGATCTCAGCTACCACATGGATATCTAAATGCATG...
```
文件内容为fasta格式即可，对文件后缀无要求。

**方式C — 位置文件（仅起止，配合 `--contig`）**

如所有目标区域位于同一contig，可通过 `--contig` 指定默认contig名：

```bash
python grna_design_pipeline.py \
    --genome genome.fasta --gff3 annot.gff3 \
    --contig ctg000250 \
    --regions targets.txt
```

targets.txt 中每行只需写 `start end label`：
```
62541730  62541863  CDS8
62542013  62542125  CDS9
```

##### 3. 参数说明

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| `--genome` | 是 | — | 基因组FASTA文件路径 |
| `--gff3` | 是 | — | GFF3注释文件路径 |
| `--regions` | 与--target二选一 | — | 目标区域位置文件 |
| `--target` | 与--regions二选一 | — | 目标序列FASTA文件（.fasta） |
| `--contig` | 否 | — | regions文件省略contig列时的默认contig |
| `--out_dir` | 否 | `.`（当前目录） | 输出目录 |
| `--out_prefix` | 否 | `grna_analysis` | 输出文件前缀 |
| `--top_n` | 否 | `10` | 输出Top N推荐gRNA数量 |

##### 4. 分析流程

```
目标区域 → PAM扫描(NGG/CCN) → 20bp gRNA候选(CCN自动取反向互补)
    ↓
全基因组特异性搜索
    ├── 10bp seed全基因组str.find() → 5'端延伸
    │     接受: exact ≥ 10bp 或 (aligned > 15bp 且 ≤ 5错配)
    └── 4bp seed在CDS串联序列搜索（补充近匹配）
          接受: ≤ 5错配
    ↓
两步威胁评估 (对每个非目标位点)
    ├── Step 1: 是否在CDS编码区？ → 否 = 安全(NoCDS)
    └── Step 2: 附近±2bp是否有PAM(NGG/CCN)？ → 否 = 安全(NoPAM)
                                               → 是 = 威胁
    ↓
排名: total_threats ↑ → exact_threats ↑ → near_threats ↑ → off_target_total ↑
    ↓
威胁注释: 每个威胁位点交叉GFF3 → 基因号 + CDS编号
    ↓
输出: TXT(完整报告) + JSON(结构化数据) + TXT(威胁位点基因注释)
```

核心参数说明：
- **seed_4bp**：gRNA靠近PAM端的最后4bp（3'端），即Cas9切割特异性最关键的区域
- **PAM判定**：正链gRNA的PAM为NGG（位于gRNA 3'端+20bp处）；反链gRNA目标为反向链，正向链上PAM为CCN（位于gRNA起始-3bp处），容差±2bp
- **on-target判定**：命中位点与目标contig相同且位置偏差≤3bp

##### 5. 输出文件

每次运行产生3个文件（以 `--out_prefix my_gene` 为例）：

| 文件 | 内容 |
|------|------|
| `my_gene_results.txt` | 完整报告：Top N详情(含基因注释) + 全排名表 + 每个候选的详细威胁列表 |
| `my_gene_results.json` | JSON格式：所有候选排序结果，每个威胁包含基因/CDS/外显子/UTR注释 |
| `my_gene_threat_annotation.txt` | Top N候选的威胁位点基因注释（基因号 + CDS编号） |

##### 6. 输出示例

```
--- Rank 1 ---
  gRNA (20bp):  AACCCAGAGAATGGGCGTCG
  Seed (4bp):   GTCG
  PAM:          GGG (strand: +)
  Region:       CDS8
  BLAST hits:   179 (OT: 178)
  Step 1 fail (non-CDS): 142
  Step 2 fail (no PAM):  29
  Exact threats:  0
  Near threats:   7
  TOTAL THREATS:  7

  Threat details:
    ctg000050:21560900  exact=12bp  mm=8  PAM=NGG:CGG@21560922
      Gene: ctg000050_000674
      CDS: 5
    ctg000070:111956806  exact=11bp  mm=9  PAM=NGG:CGG@111956824
      Gene: ctg000070_001900
      CDS: 6
```

##### 7. 性能参考

| 数据规模 | 参考时间 | 内存占用 |
|----------|----------|----------|
| ~878 Mb 基因组, 36 contigs, 25 candidates | ~7.5 分钟 | ~1 GB |
| 基因组每增加 100 Mb | 约增加 1 分钟 | — |
| 候选数每增加 10 个 | 约增加 2-3 分钟 | — |

---
