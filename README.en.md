# GuideRNA_predict

#### Description
A genome-wide CRISPR/Cas9 guideRNA prediction and specificity assessment tool for non-model species (based on sequence specificity, off-target position, and seed region conservation).

Given target regions (genomic coordinates or FASTA sequences) and a GFF3 annotation file, the pipeline scans for NGG/CCN PAM sites to extract 20bp gRNA candidates, performs genome-wide BLAST-like off-target search, applies a two-step threat assessment (in CDS? + has PAM?), and outputs specificity-ranked gRNA recommendations with gene-level annotation of off-target sites.

#### Software Architecture

```
grna_design_pipeline.py          # Main pipeline (pure Python stdlib, zero dependencies)
├── load_genome()                #   Load genome FASTA
├── parse_gff3_cds()             #   Parse GFF3: extract CDS regions
├── parse_gff3_features()        #   Parse GFF3: all features (for annotation)
├── parse_target_regions()       #   Parse target regions file (multiple formats)
├── locate_target_sequences()    #   Search target FASTA against genome
├── scan_pam_sites()             #   Scan NGG/CCN PAM sites, extract gRNA candidates
├── build_cds_data()             #   Build concatenated CDS strings (fast search)
├── analyze_specificity()        #   Genome-wide specificity analysis (core)
│   ├── 10bp seed genome search + 5' extension (exact matches)
│   ├── 4bp seed CDS search (near-match supplement)
│   └── Two-step threat assessment
│       ├── Step 1: Is off-target in a CDS region?
│       └── Step 2: Is there a PAM within ±2bp?
├── annotate_threats()           #   Annotate threats with gene ID & CDS number
└── write_outputs()              #   Output TXT + JSON results
```

#### Installation

1.  Ensure Python >= 3.6:
    ```
    python --version
    ```

2.  Download the script to your working directory. **No third-party Python packages are required** — the pipeline uses only standard library modules (argparse, json, re, sys, time, collections).

3.  Prepare input files:
    - Genome sequence (FASTA format, multi-contig supported)
    - Gene annotation (GFF3 format, must contain CDS feature lines)
    - Target regions file (coordinates) or target sequences file (FASTA)

#### Instructions

##### 1. Basic Usage

```bash
python grna_design_pipeline.py \
    --genome <genome.fasta> \
    --gff3 <annotation.gff3> \
    --regions <targets.txt> \
    --out_dir <output_directory> \
    --out_prefix <output_prefix> \
    --top_n 10
```

##### 2. Three Ways to Specify Target Regions

**Method A — Regions file (`--regions`)**

One region per line (tab or space separated). The following formats are supported (can be mixed):

```
# Full: contig + start + end + label
ctg000250  62541730  62541863  CDS8

# Simplified: contig + start + end (labels auto-generated: R1, R2, ...)
ctg000250  62541730  62541863

# Colon format (with label)
ctg000250:62541730-62541863  CDS8

# Colon format (without label, coordinates used as label)
ctg000250:62541730-62541863
```

Lines starting with `#` are treated as comments.

**Method B — Target FASTA (`--target`)**

Provide target sequences in FASTA format. The script automatically searches each sequence against the genome (both strands) to locate its position. FASTA headers are used as region labels:

```fasta
>CDS8
ATGGCAATTGCCTAACCCAGAGAATGGGCGTCGGGGACT...
>CDS9
TGGAGATCTCAGCTACCACATGGATATCTAAATGCATG...
```

**Method C — Regions with `--contig`**

When all target regions are on the same contig, use `--contig` to set a default:

```bash
python grna_design_pipeline.py \
    --genome genome.fasta --gff3 annot.gff3 \
    --contig ctg000250 \
    --regions targets.txt
```

Each line in targets.txt only needs `start end label`:
```
62541730  62541863  CDS8
62542013  62542125  CDS9
```

##### 3. Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--genome` | Yes | — | Genome FASTA file path |
| `--gff3` | Yes | — | GFF3 annotation file path |
| `--regions` | Either | — | Target regions file (coordinates) |
| `--target` | Either | — | Target sequences FASTA file |
| `--contig` | No | — | Default contig when regions omit contig column |
| `--out_dir` | No | `.` (current dir) | Output directory |
| `--out_prefix` | No | `grna_analysis` | Output file prefix |
| `--top_n` | No | `10` | Number of top gRNA candidates to report |

##### 4. Pipeline Overview

```
Target regions → PAM scan (NGG/CCN) → 20bp gRNA candidates (CCN: auto revcomp)
    ↓
Genome-wide specificity search
    ├── 10bp seed str.find() genome-wide → 5' extension
    │     Accept: exact ≥ 10bp OR (aligned > 15bp AND ≤ 5 mismatches)
    └── 4bp seed in concatenated CDS (near-match supplement)
          Accept: ≤ 5 mismatches
    ↓
Two-step threat assessment (for each off-target hit)
    ├── Step 1: In CDS region? → No = Safe (NoCDS)
    └── Step 2: PAM (NGG/CCN) within ±2bp? → No = Safe (NoPAM)
                                           → Yes = THREAT
    ↓
Ranking: total_threats ↑ → exact_threats ↑ → near_threats ↑ → off_target_total ↑
    ↓
Threat annotation: cross-reference each threat position against GFF3 → gene ID + CDS number
    ↓
Output: TXT (full report) + JSON (structured data) + TXT (threat annotation)
```

Key parameters:
- **seed_4bp**: The 4bp at the PAM-proximal end of the gRNA (3' end), the most critical region for Cas9 cleavage specificity
- **PAM detection**: Forward-strand gRNA expects NGG at gRNA 3' end (+20bp); reverse-strand gRNA targets the reverse strand (CCN on forward strand at gRNA start −3bp), tolerance ±2bp
- **On-target determination**: Hit must be on the intended contig with position deviation ≤ 3bp

##### 5. Output Files

Each run produces 3 files (example with `--out_prefix my_gene`):

| File | Content |
|------|---------|
| `my_gene_results.txt` | Full report: Top N details (with gene annotation) + full ranking table + per-candidate threat details |
| `my_gene_results.json` | JSON format: all ranked candidates with gene/CDS/exon/UTR annotation per threat |
| `my_gene_threat_annotation.txt` | Top N threat site annotation (gene ID + CDS number) |

##### 6. Example Output

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

##### 7. Performance

| Data Scale | Approximate Time | Memory |
|------------|-----------------|--------|
| ~878 Mb genome, 36 contigs, 25 candidates | ~7.5 min | ~1 GB |
| Per additional 100 Mb genome | ~1 min | — |
| Per additional 10 candidates | ~2–3 min | — |
