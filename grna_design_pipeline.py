#!/usr/bin/env python3
"""
CRISPR/Cas9 gRNA Design & Specificity Analysis Pipeline
========================================================
Two ways to specify target regions:

  A) By position (--regions):
       Regions file (one per line):
         contig  start  end  label       ← full format
         contig:start-end  label          ← colon format
         start  end  label                ← contig from --contig
         start-end  label                 ← contig from --contig

  B) By sequence (--target):
       FASTA file with target sequences. Each sequence is searched against
       the genome to locate its position. Header is used as label.

Parameters:
  --genome      Genome FASTA file
  --gff3        GFF3 annotation file
  --regions     Target regions file (positions)
  --target      Target sequences FASTA (alternative to --regions)
  --contig      Default contig for regions without contig column
  --out_dir     Output directory (default: current directory)
  --out_prefix  Output file prefix (default: "grna_analysis")
  --top_n       Number of top candidates to highlight (default: 10)

Pipeline:
  1. Scan target regions for NGG/CCN PAM sites, extract 20bp gRNA candidates
  2. Genome-wide specificity analysis:
     - Exact matches: 10-20bp PAM-proximal match (genome-wide, 10bp seed + extend)
     - Near-matches: >15bp aligned with <=5 mismatches (4bp seed in CDS)
  3. 2-step threat assessment:
     - Step 1: Is off-target in a CDS region?
     - Step 2: Is there a PAM (NGG/CCN) nearby?
  4. Rank gRNAs by specificity, annotate threats with gene/CDS info
  5. Output: summary table + detailed threat annotation (TXT + JSON)
"""

import argparse
import json
import re
import sys
import time
from collections import defaultdict


# ============================================================
# Utility functions
# ============================================================
def revcomp(s):
    return s.translate(str.maketrans("ATGCatgc", "TACGtacg"))[::-1]


def hamming_fast(s1, s2, max_mm=5):
    """Hamming distance with early exit."""
    c = 0
    for a, b in zip(s1, s2):
        if a != b:
            c += 1
            if c > max_mm:
                return c
    return c


# ============================================================
# Data loading
# ============================================================
def load_genome(fasta_path):
    """Load genome sequences from FASTA file. Returns dict: contig_id -> sequence."""
    print(f"[1/6] Loading genome: {fasta_path}")
    contigs = {}
    current_id = None
    current_seq = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    contigs[current_id] = "".join(current_seq)
                current_id = line.strip()[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            contigs[current_id] = "".join(current_seq)
    total_bp = sum(len(s) for s in contigs.values())
    print(f"  Loaded {len(contigs)} contigs, {total_bp:,} bp")
    return contigs


def parse_gff3_cds(gff_path):
    """Parse GFF3 to extract CDS regions per contig. Returns dict: contig -> [(start, end), ...]."""
    print(f"[2/6] Parsing GFF3: {gff_path}")
    cds_map = defaultdict(list)
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] == "CDS":
                cds_map[parts[0]].append((int(parts[3]), int(parts[4])))
    # Merge overlapping/adjacent CDS regions
    for ctg in cds_map:
        regions = sorted(cds_map[ctg])
        merged = []
        for r in regions:
            if merged and r[0] <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], r[1]))
            else:
                merged.append(r)
        cds_map[ctg] = merged
    total_regions = sum(len(v) for v in cds_map.values())
    print(f"  Parsed {total_regions} merged CDS regions across {len(cds_map)} contigs")
    return cds_map


def parse_gff3_features(gff_path):
    """Parse GFF3 to extract ALL features per contig for annotation."""
    features_by_ctg = defaultdict(list)
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            ctg = parts[0]
            ftype = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = {}
            for a in parts[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k] = v
            fid = attrs.get("ID", "")
            parent = attrs.get("Parent", "")
            features_by_ctg[ctg].append((start, end, ftype, fid, parent))

            if ftype == "gene":
                m = re.search(r'(ctg\d+)_(\d+)$', fid)
                if m:
                    genes[fid] = {"ctg": ctg, "start": start, "end": end,
                                  "short": f"{m.group(1)}_{m.group(2)}"}
    return features_by_ctg, genes


def parse_target_regions(regions_file, default_contig=None):
    """Parse target regions file.

    Accepted formats (one per line):
      contig  start  end  label       ← 4 columns
      contig  start  end              ← 3 columns, label auto-generated (R1, R2...)
      contig:start-end  label         ← colon format with label
      contig:start-end                ← colon format, label from coordinates
      start  end  label               ← 3 columns, contig from --contig
      start-end  label                ← compact, contig from --contig
    """
    print(f"[3/6] Parsing target regions: {regions_file}")
    if default_contig:
        print(f"  Default contig: {default_contig}")

    regions = []
    auto_idx = 0
    with open(regions_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) == 1:
                parts = line.split()

            matched = None
            auto_idx += 1

            def _is_int(s):
                try:
                    int(s)
                    return True
                except ValueError:
                    return False

            # ── 4 columns: "contig start end label" ──
            if len(parts) >= 4 and matched is None:
                try:
                    matched = {
                        "contig": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "label": " ".join(parts[3:]),
                    }
                except ValueError:
                    pass

            # ── 3 columns ──
            if matched is None and len(parts) == 3:
                a, b, c = parts[0], parts[1], parts[2]
                # "contig start end" → col1 is non-numeric (contig name), col2/3 are ints
                if not _is_int(a) and _is_int(b) and _is_int(c):
                    matched = {
                        "contig": a,
                        "start": int(b),
                        "end": int(c),
                        "label": f"R{auto_idx}",
                    }
                # "start end label" → col1/2 are ints, col3 is label, contig from --contig
                elif _is_int(a) and _is_int(b):
                    matched = {
                        "contig": default_contig,
                        "start": int(a),
                        "end": int(b),
                        "label": c,
                    }

            # ── 2 columns: "contig:start-end label" or "start-end label" ──
            if matched is None and len(parts) == 2:
                # "contig:start-end label"
                m = re.match(r'(\S+):(\d+)-(\d+)', parts[0])
                if m:
                    matched = {
                        "contig": m.group(1),
                        "start": int(m.group(2)),
                        "end": int(m.group(3)),
                        "label": parts[1],
                    }
                else:
                    # "start-end label" (compact, no contig)
                    m = re.match(r'(\d+)-(\d+)', parts[0])
                    if m:
                        matched = {
                            "contig": default_contig,
                            "start": int(m.group(1)),
                            "end": int(m.group(2)),
                            "label": parts[1],
                        }

            # ── 1 column: "contig:start-end" (label from coordinates) ──
            if matched is None and len(parts) == 1:
                m = re.match(r'(\S+):(\d+)-(\d+)', line)
                if m:
                    matched = {
                        "contig": m.group(1),
                        "start": int(m.group(2)),
                        "end": int(m.group(3)),
                        "label": f"{m.group(2)}-{m.group(3)}",
                    }

            if matched:
                if matched["contig"] is None:
                    print(f"  ERROR: no contig for line: {line}")
                    print(f"  Use --contig to specify a default contig.")
                    sys.exit(1)
                regions.append(matched)
            else:
                print(f"  WARNING: could not parse line: {line}")

    print(f"  Found {len(regions)} target regions")
    for r in regions:
        print(f"    {r['contig']}:{r['start']}-{r['end']}  [{r['label']}]  ({r['end']-r['start']+1}bp)")
    return regions


def locate_target_sequences(contigs, target_fasta):
    """Search target FASTA sequences against genome to find their positions.
    Each sequence's header is used as the region label.
    Returns list of target regions (same format as parse_target_regions).
    """
    print(f"[3/6] Locating target sequences from: {target_fasta}")

    # Load target sequences
    targets = []
    current_id = None
    current_seq = []
    with open(target_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if current_id and current_seq:
                    targets.append((current_id, "".join(current_seq).upper()))
                current_id = line.strip()[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id and current_seq:
            targets.append((current_id, "".join(current_seq).upper()))

    print(f"  Loaded {len(targets)} target sequences")

    regions = []
    for label, seq in targets:
        found = False
        for ctg, ctg_seq in contigs.items():
            pos = ctg_seq.upper().find(seq)
            if pos != -1:
                regions.append({
                    "contig": ctg,
                    "start": pos + 1,
                    "end": pos + len(seq),
                    "label": label,
                })
                print(f"    {label}: found on {ctg}:{pos+1}-{pos+len(seq)} ({len(seq)}bp)")
                found = True
                break
        if not found:
            # Try reverse complement
            seq_rc = revcomp(seq)
            for ctg, ctg_seq in contigs.items():
                pos = ctg_seq.upper().find(seq_rc)
                if pos != -1:
                    regions.append({
                        "contig": ctg,
                        "start": pos + 1,
                        "end": pos + len(seq),
                        "label": label,
                    })
                    print(f"    {label}: found (revcomp) on {ctg}:{pos+1}-{pos+len(seq)} ({len(seq)}bp)")
                    found = True
                    break
        if not found:
            print(f"    WARNING: {label} not found in genome!")

    return regions


# ============================================================
# Step 1: PAM scanning
# ============================================================
def scan_pam_sites(contigs, target_regions):
    """Scan target regions for NGG (forward) and CCN (reverse) PAM sites.
    Returns list of gRNA candidates.
    """
    print("  Scanning for PAM sites...")
    candidates = []

    for region in target_regions:
        ctg = region["contig"]
        if ctg not in contigs:
            print(f"    WARNING: contig {ctg} not in genome, skipping")
            continue

        ctg_seq = contigs[ctg]
        r_start = region["start"]
        r_end = region["end"]
        label = region["label"]

        # Extract region sequence (1-based to 0-based)
        region_seq = ctg_seq[r_start - 1:r_end]

        # Scan for NGG (forward strand)
        for i in range(len(region_seq) - 2):
            triplet = region_seq[i:i+3].upper()
            if triplet[1:3] == "GG" and triplet[0] in "ATGC":
                pam_start = r_start + i
                # gRNA is 20bp upstream of PAM
                grna_start = pam_start - 20
                if grna_start >= 1:
                    grna_seq = ctg_seq[grna_start - 1:pam_start - 1]
                    if len(grna_seq) == 20:
                        candidates.append({
                            "label": label,
                            "strand": "+",
                            "pam": triplet,
                            "pam_start": pam_start,
                            "pam_end": pam_start + 2,
                            "grna_start": grna_start,
                            "grna_end": pam_start - 1,
                            "grna_seq": grna_seq.upper(),
                        })

        # Scan for CCN (reverse strand) — CCN on forward strand means NGG on reverse
        for i in range(len(region_seq) - 2):
            triplet = region_seq[i:i+3].upper()
            if triplet[0:2] == "CC" and triplet[2] in "ATGN":
                pam_start = r_start + i
                # gRNA is 20bp downstream of CCN on forward strand
                # but since Cas9 cuts reverse strand, gRNA = revcomp(downstream 20bp)
                grna_start = pam_start + 3
                if grna_start + 19 <= r_end:
                    grna_fwd = ctg_seq[grna_start - 1:grna_start + 19]
                    if len(grna_fwd) == 20:
                        candidates.append({
                            "label": label,
                            "strand": "-",
                            "pam": triplet,
                            "pam_start": pam_start,
                            "pam_end": pam_start + 2,
                            "grna_start": grna_start,
                            "grna_end": grna_start + 19,
                            "grna_seq": grna_fwd.upper(),  # will be revcomp'd later
                        })

    # Deduplicate (same gRNA sequence in same region)
    seen = set()
    unique = []
    for c in candidates:
        key = (c["grna_seq"], c["strand"], c["label"])
        if key not in seen:
            seen.add(key)
            unique.append(c)

    print(f"  Found {len(candidates)} raw PAM sites, {len(unique)} unique candidates")
    return unique


# ============================================================
# Step 2: Build CDS data structures
# ============================================================
def build_cds_data(contigs, cds_map):
    """Build concatenated CDS strings for each contig (for efficient seed search)."""
    cds_data = {}
    for ctg, cds_list in cds_map.items():
        if ctg not in contigs:
            continue
        parts = []
        pos_map = []
        offset = 0
        for s, e in cds_list:
            seq = contigs[ctg][s-1:e]
            pos_map.append((offset, s, e))
            parts.append(seq)
            offset += len(seq)
        cds_data[ctg] = ("".join(parts), pos_map)
    return cds_data


def offset_to_genomic(ctg, offset, cds_data):
    """Convert CDS concatenated offset back to genomic position."""
    if ctg not in cds_data:
        return None
    for concat_start, gen_start, gen_end in cds_data[ctg][1]:
        cds_len = gen_end - gen_start + 1
        if concat_start <= offset < concat_start + cds_len:
            return gen_start + (offset - concat_start)
    return None


def is_in_cds(contig, pos, cds_map):
    """Check if a genomic position falls within any CDS region."""
    if contig not in cds_map:
        return False
    for s, e in cds_map[contig]:
        if s <= pos <= e:
            return True
    return False


# ============================================================
# Step 3: PAM checking
# ============================================================
def check_pam(ctg_seq, pos_0, strand):
    """Check for PAM at expected position (±2bp tolerance).
    + strand: NGG at pos_0+20 (3' of gRNA binding site)
    - strand: CCN at pos_0-3 (5' of gRNA binding site on fwd strand)
    """
    n = len(ctg_seq)
    if strand == "+":
        for off in [0, -1, 1, -2, 2]:
            p = pos_0 + 20 + off
            if 0 <= p <= n - 3:
                tri = ctg_seq[p:p+3].upper()
                if tri[1:3] == "GG" and tri[0] in "ATGC":
                    return True, f"NGG:{tri}@{p+1}"
    else:
        for off in [0, -1, 1, -2, 2]:
            p = pos_0 - 3 + off
            if 0 <= p <= n - 3:
                tri = ctg_seq[p:p+3].upper()
                if tri[0:2] == "CC" and tri[2] in "ATGN":
                    return True, f"CCN:{tri}@{p+1}"
    return False, None


# ============================================================
# Step 4: Match extension (for 10bp seed hits)
# ============================================================
def try_extend_match(genome_seq, pos_0, gRNA_20):
    """Given a hit of seed_10bp (last 10bp of gRNA) at genome[pos_0+10:pos_0+20],
    extend towards 5' end. Returns (exact_len, total_aligned, mm)."""
    target = genome_seq[pos_0:pos_0+20]
    exact_len = 10
    mm = 0
    for i in range(9, -1, -1):
        if i >= len(target):
            break
        if target[i].upper() == gRNA_20[i].upper():
            exact_len += 1
        else:
            mm += 1

    total_aligned = 20
    mm_5prime = 0
    for i in range(10):
        if i >= len(target):
            total_aligned = i
            break
        if target[i].upper() != gRNA_20[i].upper():
            mm_5prime += 1
        if mm_5prime > 5:
            total_aligned = 10 + i
            break

    return exact_len, total_aligned, mm


def try_extend_match_rev(genome_seq, pos_0, gRNA_20):
    """Same but target = revcomp of genome window (for reverse strand hits)."""
    target = revcomp(genome_seq[pos_0:pos_0+20])
    exact_len = 10
    mm = 0
    for i in range(9, -1, -1):
        if i >= len(target):
            break
        if target[i].upper() == gRNA_20[i].upper():
            exact_len += 1
        else:
            mm += 1

    total_aligned = 20
    mm_5prime = 0
    for i in range(10):
        if i >= len(target):
            total_aligned = i
            break
        if target[i].upper() != gRNA_20[i].upper():
            mm_5prime += 1
        if mm_5prime > 5:
            total_aligned = 10 + i
            break

    return exact_len, total_aligned, mm


# ============================================================
# Step 5: Per-candidate specificity analysis
# ============================================================
def analyze_specificity(candidates, contigs, cds_map, cds_data, intended_contig, top_n=10):
    """Run genome-wide specificity analysis for each gRNA candidate."""
    print(f"[4/6] Analyzing specificity for {len(candidates)} candidates...")

    # Fix CCN gRNA and compute seeds
    for c in candidates:
        if c["strand"] == "-":
            c["grna_seq_corrected"] = revcomp(c["grna_seq"])
        else:
            c["grna_seq_corrected"] = c["grna_seq"]
        c["seed_4bp"] = c["grna_seq_corrected"][-4:]
        c["seed_10bp"] = c["grna_seq_corrected"][-10:]

    all_results = []

    for ci, cand in enumerate(candidates):
        t1 = time.time()
        grna_20 = cand["grna_seq_corrected"]
        strand = cand["strand"]
        label = cand["label"]
        seed_10 = cand["seed_10bp"]
        seed_4 = cand["seed_4bp"]
        intended_pos = cand["grna_start"]

        all_hits = {}

        # ---- Exact + Near-match search via 10bp seed (genome-wide) ----
        seed_10_rc = revcomp(seed_10)

        for ctg, ctg_seq in contigs.items():
            # Forward strand: find seed_10
            pos = 0
            while True:
                pos = ctg_seq.find(seed_10, pos)
                if pos == -1:
                    break
                g_start = pos - 10
                if g_start >= 0:
                    exact_len, total_aligned, mm = try_extend_match(ctg_seq, g_start, grna_20)
                    if exact_len >= 10 or (total_aligned > 15 and mm <= 5):
                        key = (ctg, g_start + 1, "+")
                        if key not in all_hits or all_hits[key]["exact_len"] < exact_len:
                            all_hits[key] = {
                                "contig": ctg, "pos_0": g_start, "pos_1": g_start + 1,
                                "exact_len": exact_len, "total_aligned": total_aligned,
                                "mm": mm, "strand_hit": "+",
                            }
                pos += 1

            # Reverse strand: find seed_10_rc
            pos = 0
            while True:
                pos = ctg_seq.find(seed_10_rc, pos)
                if pos == -1:
                    break
                g_start = pos - 10
                if g_start >= 0:
                    exact_len, total_aligned, mm = try_extend_match_rev(ctg_seq, g_start, grna_20)
                    if exact_len >= 10 or (total_aligned > 15 and mm <= 5):
                        key = (ctg, g_start + 1, "-")
                        if key not in all_hits or all_hits[key]["exact_len"] < exact_len:
                            all_hits[key] = {
                                "contig": ctg, "pos_0": g_start, "pos_1": g_start + 1,
                                "exact_len": exact_len, "total_aligned": total_aligned,
                                "mm": mm, "strand_hit": "-",
                            }
                pos += 1

        # ---- Additional near-match search via 4bp seed in CDS ----
        seed_4_rc = revcomp(seed_4)
        for ctg, (cds_seq, pos_map) in cds_data.items():
            n_cds = len(cds_seq)

            # Forward
            spos = 0
            while True:
                spos = cds_seq.find(seed_4, spos)
                if spos == -1:
                    break
                g_start = spos - 16
                if g_start >= 0 and g_start + 20 <= n_cds:
                    window = cds_seq[g_start:g_start+20]
                    mm = hamming_fast(grna_20, window, 5)
                    if mm <= 5:
                        # Compute exact_len from 3' end
                        exact_len = 0
                        for i in range(19, -1, -1):
                            if window[i].upper() == grna_20[i].upper():
                                exact_len += 1
                            else:
                                break
                        gen_pos = offset_to_genomic(ctg, g_start, cds_data)
                        if gen_pos is not None:
                            key = (ctg, gen_pos, "+")
                            if key not in all_hits:
                                all_hits[key] = {
                                    "contig": ctg, "pos_0": gen_pos - 1, "pos_1": gen_pos,
                                    "exact_len": exact_len, "total_aligned": 20,
                                    "mm": mm, "strand_hit": "+",
                                }
                spos += 1

            # Reverse
            spos = 0
            while True:
                spos = cds_seq.find(seed_4_rc, spos)
                if spos == -1:
                    break
                g_start = spos - 16
                if g_start >= 0 and g_start + 20 <= n_cds:
                    window = revcomp(cds_seq[g_start:g_start+20])
                    mm = hamming_fast(grna_20, window, 5)
                    if mm <= 5:
                        exact_len = 0
                        for i in range(19, -1, -1):
                            if window[i].upper() == grna_20[i].upper():
                                exact_len += 1
                            else:
                                break
                        gen_pos = offset_to_genomic(ctg, g_start, cds_data)
                        if gen_pos is not None:
                            key = (ctg, gen_pos, "-")
                            if key not in all_hits:
                                all_hits[key] = {
                                    "contig": ctg, "pos_0": gen_pos - 1, "pos_1": gen_pos,
                                    "exact_len": exact_len, "total_aligned": 20,
                                    "mm": mm, "strand_hit": "-",
                                }
                spos += 1

        # ---- 2-step threat assessment ----
        on_target = 0
        threats_exact = []
        threats_near = []
        safe_no_cds = []
        safe_no_pam = []

        for key, hit in all_hits.items():
            ctg, pos_1, str_hit = key
            pos_0 = hit["pos_0"]
            exact_len = hit["exact_len"]
            mm = hit["mm"]

            # Check if on-target
            if ctg == intended_contig:
                # For + strand, grna_start should match pos_1
                # For - strand, grna_start also matches pos_1 (genomic position)
                if abs(pos_1 - intended_pos) <= 3:
                    on_target += 1
                    continue

            # Step 1: Is match in CDS?
            in_cds = False
            for gp in range(pos_1, pos_1 + 20):
                if is_in_cds(ctg, gp, cds_map):
                    in_cds = True
                    break

            if not in_cds:
                safe_no_cds.append(f"{ctg}:{pos_1}")
                continue

            # Step 2: Has PAM nearby?
            ctg_seq = contigs[ctg]
            has_pam, pam_detail = check_pam(ctg_seq, pos_0, strand)
            if not has_pam:
                safe_no_pam.append(f"{ctg}:{pos_1}")
                continue

            threat = {
                "contig": ctg, "pos": pos_1, "strand_hit": str_hit,
                "exact_len": exact_len, "total_aligned": hit["total_aligned"],
                "mismatches": mm, "pam_detail": pam_detail,
            }

            if exact_len >= 20 and mm == 0:
                threats_exact.append(threat)
            else:
                threats_near.append(threat)

        total_hits = len(all_hits)
        total_ot = total_hits - on_target
        total_threats = len(threats_exact) + len(threats_near)

        result = {
            "idx": ci + 1,
            "label": label,
            "strand": strand,
            "pam": cand["pam"],
            "grna_seq": grna_20,
            "seed_4bp": seed_4,
            "grna_start": cand["grna_start"],
            "pam_start": cand["pam_start"],
            "total_hits": total_hits,
            "on_target": on_target,
            "off_target_total": total_ot,
            "safe_no_cds": len(safe_no_cds),
            "safe_no_pam": len(safe_no_pam),
            "exact_threats": len(threats_exact),
            "near_threats": len(threats_near),
            "total_threats": total_threats,
            "threat_details": threats_exact + threats_near,
        }
        all_results.append(result)

        elapsed = time.time() - t1
        bar = "#" if ci < 5 else "."
        print(f"  {bar} [{ci+1:2d}/{len(candidates)}] {grna_20} [{label}] "
              f"hits={total_hits} OT={total_ot} noCDS={len(safe_no_cds)} "
              f"noPAM={len(safe_no_pam)} threats={total_threats} ({elapsed:.1f}s)")

    return all_results


# ============================================================
# Step 6: Threat annotation against GFF3
# ============================================================
def annotate_threats(all_results, features_by_ctg, genes, top_n=10):
    """Annotate threat positions with gene ID, CDS number, exon info."""
    print(f"[5/6] Annotating threat positions against GFF3...")

    def find_overlapping(pos, ctg):
        results = []
        if ctg not in features_by_ctg:
            return results
        for start, end, ftype, fid, parent in features_by_ctg[ctg]:
            if start <= pos <= end:
                results.append((ftype, fid, parent, start, end))
        return results

    annotated = []
    for r in all_results:
        for t in r["threat_details"]:
            overlaps = find_overlapping(t["pos"], t["contig"])
            genes_found = [o for o in overlaps if o[0] == "gene"]
            cdss_found = [o for o in overlaps if o[0] == "CDS"]
            exons_found = [o for o in overlaps if o[0] == "exon"]
            utrs_found = [o for o in overlaps if "UTR" in o[0]]

            entry = {
                "rank": 0,  # will be filled after sorting
                "grna_seq": r["grna_seq"],
                "pam": r["pam"],
                "cds_label": r["label"],
                "contig": t["contig"],
                "position": t["pos"],
                "exact_len": t["exact_len"],
                "total_aligned": t["total_aligned"],
                "mm": t["mismatches"],
                "pam_detail": t["pam_detail"],
                "annotation": {
                    "genes": [{"id": g[1], "short": genes.get(g[1], {}).get("short", g[1]),
                               "range": [g[3], g[4]]} for g in genes_found],
                    "cds": [{"id": c[1],
                             "number": re.search(r'CDS\.(\d+)', c[1]).group(1) if re.search(r'CDS\.(\d+)', c[1]) else "?",
                             "range": [c[3], c[4]]} for c in cdss_found],
                    "exons": [{"id": e[1],
                               "number": re.search(r'exon\.(\d+)', e[1]).group(1) if re.search(r'exon\.(\d+)', e[1]) else "?",
                               "range": [e[3], e[4]]} for e in exons_found],
                    "utrs": [{"type": u[0], "id": u[1], "range": [u[3], u[4]]} for u in utrs_found],
                }
            }
            annotated.append(entry)

    return annotated


# ============================================================
# Step 7: Output
# ============================================================
def write_outputs(all_results, annotated_threats, features_by_ctg, genes,
                  out_dir, out_prefix, top_n, target_regions):
    """Write all output files."""
    print(f"[6/6] Writing output files...")

    import os
    os.makedirs(out_dir, exist_ok=True)

    # Sort: by total_threats, exact_threats, near_threats, off_target_total
    all_results.sort(key=lambda r: (r["total_threats"], r["exact_threats"],
                                     r["near_threats"], r["off_target_total"]))

    # Build rank lookup for annotation
    grna_to_rank = {}
    for rank_i, r in enumerate(all_results):
        grna_to_rank[r["grna_seq"]] = rank_i + 1

    for a in annotated_threats:
        a["rank"] = grna_to_rank.get(a["grna_seq"], 0)

    # Group annotations by gRNA
    ann_by_grna = defaultdict(list)
    for a in annotated_threats:
        ann_by_grna[a["grna_seq"]].append(a)

    # ============================
    # Output 1: Main results (TXT)
    # ============================
    txt_path = f"{out_dir}/{out_prefix}_results.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write("=" * 100 + "\n")
        f.write("CRISPR/Cas9 gRNA SPECIFICITY ANALYSIS\n")
        f.write("=" * 100 + "\n")
        f.write(f"Target regions:\n")
        for r in target_regions:
            f.write(f"  {r['contig']}:{r['start']}-{r['end']}  [{r['label']}]\n")
        f.write(f"\nMethod:\n")
        f.write(f"  PAM: NGG (forward strand) / CCN (reverse strand)\n")
        f.write(f"  gRNA: 20bp adjacent to PAM\n")
        f.write(f"  CCN fix: gRNA = revcomp(downstream 20bp)\n")
        f.write(f"  Exact search: 10-20bp PAM-proximal match (genome-wide)\n")
        f.write(f"  Near-match search: >15bp aligned, <=5 mismatches (4bp seed in CDS)\n")
        f.write(f"  Step 1: Is off-target in CDS?\n")
        f.write(f"  Step 2: Has PAM (NGG/CCN) within +/-2bp?\n")
        f.write(f"\nTotal candidates analyzed: {len(all_results)}\n\n")

        # Top N candidates
        f.write("=" * 100 + "\n")
        f.write(f"TOP {top_n} SPECIFIC gRNA CANDIDATES\n")
        f.write("=" * 100 + "\n")

        for rank_i, r in enumerate(all_results[:top_n]):
            f.write(f"\n--- Rank {rank_i + 1} ---\n")
            f.write(f"  gRNA (20bp):  {r['grna_seq']}\n")
            f.write(f"  Seed (4bp):   {r['seed_4bp']}\n")
            f.write(f"  PAM:          {r['pam']} (strand: {r['strand']})\n")
            f.write(f"  Region:       {r['label']}\n")
            f.write(f"  Position:     gRNA {r['grna_start']}-{r['grna_start']+19}, PAM {r['pam_start']}-{r['pam_start']+2}\n")
            f.write(f"  BLAST hits:   {r['total_hits']} (OT: {r['off_target_total']})\n")
            f.write(f"  Step 1 fail (non-CDS): {r['safe_no_cds']}\n")
            f.write(f"  Step 2 fail (no PAM):  {r['safe_no_pam']}\n")
            f.write(f"  Exact threats:  {r['exact_threats']}\n")
            f.write(f"  Near threats:   {r['near_threats']}\n")
            f.write(f"  TOTAL THREATS:  {r['total_threats']}\n")

            if r["threat_details"]:
                f.write(f"\n  Threat details:\n")
                for t in r["threat_details"]:
                    f.write(f"    {t['contig']}:{t['pos']}  exact={t['exact_len']}bp  mm={t['mismatches']}  "
                            f"PAM={t['pam_detail']}\n")
                    # Add annotation
                    key = r["grna_seq"]
                    ann_matches = [a for a in ann_by_grna.get(key, [])
                                   if a["contig"] == t["contig"] and a["position"] == t["pos"]]
                    for am in ann_matches:
                        for g in am["annotation"]["genes"]:
                            f.write(f"      Gene: {g['short']}  [{g['range'][0]}-{g['range'][1]}]\n")
                        for cds in am["annotation"]["cds"]:
                            f.write(f"      CDS: {cds['number']}  [{cds['range'][0]}-{cds['range'][1]}]\n")
            else:
                f.write(f"\n  VERDICT: NO THREATS - HIGHLY SPECIFIC\n")

        # Full ranking table
        f.write(f"\n\n{'=' * 110}\n")
        f.write("FULL RANKING\n")
        f.write(f"{'=' * 110}\n")
        hdr = (f"{'Rank':<5} {'gRNA':<25} {'PAM':<5} {'Hits':<6} {'OT':<5} "
               f"{'NoCDS':<7} {'NoPAM':<7} {'Th-Ex':<7} {'Th-Nr':<7} {'CDS':<8}")
        f.write(hdr + "\n")
        f.write("-" * 110 + "\n")
        for rank_i, r in enumerate(all_results):
            f.write(f"{rank_i+1:<5} {r['grna_seq']:<25} {r['pam']:<5} {r['total_hits']:<6} "
                    f"{r['off_target_total']:<5} {r['safe_no_cds']:<7} {r['safe_no_pam']:<7} "
                    f"{r['exact_threats']:<7} {r['near_threats']:<7} {r['label']:<8}\n")

        # Detailed per-candidate section
        f.write(f"\n\n{'=' * 110}\n")
        f.write("DETAILED (all candidates with threat lines)\n")
        f.write(f"{'=' * 110}\n")
        for r in all_results:
            f.write(f"\n#{r['idx']} {r['grna_seq']} | PAM={r['pam']}({r['strand']}) | {r['label']}\n")
            f.write(f"  Seed4: {r['seed_4bp']}, Hits: {r['total_hits']}, OT: {r['off_target_total']}, "
                    f"NoCDS: {r['safe_no_cds']}, NoPAM: {r['safe_no_pam']}, "
                    f"Th-ex: {r['exact_threats']}, Th-nr: {r['near_threats']}\n")
            if r["threat_details"]:
                for t in r["threat_details"]:
                    f.write(f"  !! THREAT: {t['contig']}:{t['pos']} "
                            f"exact={t['exact_len']}bp mm={t['mismatches']} "
                            f"PAM={t['pam_detail']}\n")

        # Threat annotation summary
        f.write(f"\n\n{'=' * 100}\n")
        f.write("THREAT SITE ANNOTATION SUMMARY\n")
        f.write(f"{'=' * 100}\n")
        for rank_i, r in enumerate(all_results[:top_n]):
            threats = ann_by_grna.get(r["grna_seq"], [])
            f.write(f"\nRank {rank_i+1}: {r['grna_seq']} ({r['label']}) — {len(threats)} threats\n")
            for a in threats:
                genes_str = ", ".join(g["short"] for g in a["annotation"]["genes"]) or "INTERGENIC"
                cds_parts = [f"CDS.{c['number']}" for c in a["annotation"]["cds"]]
                cds_str = ", ".join(cds_parts) if cds_parts else (
                    "UTR" if a["annotation"]["utrs"] else "INTERGENIC")
                f.write(f"  {a['contig']}:{a['position']} -> {genes_str} | {cds_str}\n")

    print(f"  Saved: {txt_path}")

    # ============================
    # Output 2: JSON
    # ============================
    json_path = f"{out_dir}/{out_prefix}_results.json"
    with open(json_path, "w", encoding="utf-8") as f:
        # Build output with rankings
        ranked_results = []
        for rank_i, r in enumerate(all_results):
            entry = dict(r)
            entry["rank"] = rank_i + 1
            # Add annotations to each threat
            key = r["grna_seq"]
            ann_map = {}
            for a in ann_by_grna.get(key, []):
                ann_map[(a["contig"], a["position"])] = a["annotation"]
            for t in entry["threat_details"]:
                t["annotation"] = ann_map.get((t["contig"], t["pos"]), {})
            ranked_results.append(entry)

        json.dump({
            "parameters": {
                "method": "10-20bp exact + >15bp near-match (<=5mm), 4bp seed, 2-step assessment",
                "pam_forward": "NGG",
                "pam_reverse": "CCN",
                "ccn_fix": "gRNA = revcomp(downstream 20bp)",
            },
            "target_regions": target_regions,
            "total_candidates": len(all_results),
            "top_n": top_n,
            "results": ranked_results,
        }, f, indent=2, ensure_ascii=False)
    print(f"  Saved: {json_path}")

    # ============================
    # Output 3: Threat annotation only (TXT)
    # ============================
    ann_path = f"{out_dir}/{out_prefix}_threat_annotation.txt"
    with open(ann_path, "w", encoding="utf-8") as f:
        f.write("=" * 100 + "\n")
        f.write("gRNA THREAT SITE ANNOTATION\n")
        f.write("=" * 100 + "\n\n")

        for rank_i, r in enumerate(all_results[:top_n]):
            threats = ann_by_grna.get(r["grna_seq"], [])
            f.write(f"Rank {rank_i+1}: {r['grna_seq']} | PAM={r['pam']} | {r['label']}\n")
            f.write(f"  Total threats: {len(threats)}\n\n")

            for i, a in enumerate(threats, 1):
                f.write(f"  Threat #{i}: {a['contig']}:{a['position']}  "
                        f"exact={a['exact_len']}bp  mm={a['mm']}  PAM={a['pam_detail']}\n")
                for g in a["annotation"]["genes"]:
                    f.write(f"    Gene: {g['short']}  ({g['id']})  [{g['range'][0]}-{g['range'][1]}]\n")
                if not a["annotation"]["genes"]:
                    f.write(f"    Gene: [NONE - intergenic]\n")
                for cds in a["annotation"]["cds"]:
                    f.write(f"    CDS: {cds['number']}  [{cds['range'][0]}-{cds['range'][1]}]\n")
                for ex in a["annotation"]["exons"]:
                    f.write(f"    Exon: {ex['number']}  [{ex['range'][0]}-{ex['range'][1]}]\n")
                for u in a["annotation"]["utrs"]:
                    f.write(f"    {u['type']}  [{u['range'][0]}-{u['range'][1]}]\n")
                f.write("\n")
            f.write("\n")

        # Summary table
        f.write("=" * 100 + "\n")
        f.write("SUMMARY\n")
        f.write("=" * 100 + "\n\n")
        for rank_i, r in enumerate(all_results[:top_n]):
            threats = ann_by_grna.get(r["grna_seq"], [])
            f.write(f"Rank {rank_i+1}: {r['grna_seq']} ({r['label']}) — {len(threats)} threats\n")
            for a in threats:
                genes_str = ", ".join(g["short"] for g in a["annotation"]["genes"]) or "INTERGENIC"
                cds_parts = [f"CDS.{c['number']}" for c in a["annotation"]["cds"]]
                cds_str = ", ".join(cds_parts) if cds_parts else (
                    "UTR" if a["annotation"]["utrs"] else "INTERGENIC")
                f.write(f"  {a['contig']}:{a['position']} -> {genes_str} | {cds_str}\n")
            f.write("\n")

    print(f"  Saved: {ann_path}")

    return txt_path, json_path, ann_path


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description="CRISPR/Cas9 gRNA Design & Specificity Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # By position (full format with label)
  python grna_design_pipeline.py --genome genome.fasta --gff3 annot.gff3 \\
      --regions targets.txt --out_dir ./results --out_prefix my_gene

  # By position (contig + start/end, label auto-generated)
  python grna_design_pipeline.py --genome genome.fasta --gff3 annot.gff3 \\
      --regions targets.txt --out_dir ./results

  # By sequence (FASTA)
  python grna_design_pipeline.py --genome genome.fasta --gff3 annot.gff3 \\
      --target cds_sequences.fasta --out_dir ./results --out_prefix my_gene

  # Regions file formats (one per line):
      ctg000250  62541730  62541863  CDS8       ← full: contig start end label
      ctg000250  62541730  62541863             ← simplified: contig start end (auto R1,R2...)
      ctg000250:62541730-62541863  CDS8          ← colon with label
      ctg000250:62541730-62541863               ← colon, label from coords
        """
    )
    parser.add_argument("--genome", required=True, help="Genome FASTA file")
    parser.add_argument("--gff3", required=True, help="GFF3 annotation file")
    parser.add_argument("--regions", default=None,
                        help="Target regions file (positions). Use --target for sequence-based input.")
    parser.add_argument("--target", default=None,
                        help="Target sequences FASTA. Each sequence is searched against genome to locate positions.")
    parser.add_argument("--contig", default=None,
                        help="Default contig for --regions lines that omit contig column.")
    parser.add_argument("--out_dir", default=".", help="Output directory (default: .)")
    parser.add_argument("--out_prefix", default="grna_analysis", help="Output file prefix (default: grna_analysis)")
    parser.add_argument("--top_n", type=int, default=10, help="Number of top candidates to highlight (default: 10)")

    args = parser.parse_args()

    if not args.regions and not args.target:
        print("ERROR: Either --regions or --target must be specified.")
        sys.exit(1)

    t_total = time.time()

    # Load data
    contigs = load_genome(args.genome)
    cds_map = parse_gff3_cds(args.gff3)
    features_by_ctg, genes = parse_gff3_features(args.gff3)

    # Get target regions (by position or by sequence search)
    if args.target:
        target_regions = locate_target_sequences(contigs, args.target)
        if not target_regions:
            print("ERROR: No target sequences found in genome!")
            sys.exit(1)
    else:
        target_regions = parse_target_regions(args.regions, args.contig)

    # Scan PAM sites
    print(f"[3/6] Scanning target regions for PAM sites...")
    candidates = scan_pam_sites(contigs, target_regions)
    if not candidates:
        print("ERROR: No PAM sites found in target regions!")
        sys.exit(1)

    # Build CDS data
    cds_data = build_cds_data(contigs, cds_map)

    # Determine intended contig
    intended_contig = target_regions[0]["contig"]

    # Specificity analysis
    all_results = analyze_specificity(candidates, contigs, cds_map, cds_data,
                                       intended_contig, args.top_n)

    # Annotate threats
    annotated = annotate_threats(all_results, features_by_ctg, genes, args.top_n)

    # Write outputs
    txt_path, json_path, ann_path = write_outputs(
        all_results, annotated, features_by_ctg, genes,
        args.out_dir, args.out_prefix, args.top_n, target_regions
    )

    total_elapsed = time.time() - t_total
    print(f"\n{'=' * 60}")
    print(f"Analysis complete! Total time: {total_elapsed:.1f}s")
    print(f"Output files:")
    print(f"  {txt_path}")
    print(f"  {json_path}")
    print(f"  {ann_path}")


if __name__ == "__main__":
    main()
