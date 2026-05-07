#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ==============================================================================
# [Script]: 1.1.1.1_prepare_GWAS_LAVA_inputs_template.py
# [Method]: LAVA (Local Analysis of [co]Variant Associations)
# [Step]: Data Preprocessing 
# 
# [Function]:
# A unified wrapper to preprocess GWAS summary statistics for LAVA. 
#       Performs strict MAF filtering, palindromic SNP removal, and allele alignment 
#       against a reference LD panel, followed by LDSC munge_sumstats formatting.
#       MAF、SNP、, 
#        LDSC  munge_sumstats.py . 
# 
# [Data Availability / ]:
# Ensure you have the corresponding GWAS summary statistics and LAVA UKB LD reference.
# ==============================================================================

from __future__ import annotations
import argparse
import csv
import datetime as dt
import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# -----------------------------
# Environment / 
# -----------------------------
def sanitize_env() -> dict:
    """
    EN: Create a sanitized environment for subprocess calls.
    ZH: . 
    """
    keep_keys = {"PATH", "HOME", "USER", "SHELL", "TMPDIR", "LANG", "LC_ALL"}
    env = {k: v for k, v in os.environ.items() if k in keep_keys}
    env.setdefault("LANG", "en_US.UTF-8")
    env.setdefault("LC_ALL", env["LANG"])
    return env

def init_output_structure(out_base: Path, clean_name: str) -> dict:
    """
    EN: Initialize timestamped folder with standard subdirectories.
    ZH: . 
    """
    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    root = out_base / f"{clean_name}_{ts}"
    dirs = {
        "root": root,
        "analysis": root / "analysis",
        "logs": root / "logs",
        "reports": root / "reports",
        "readme": root / "readme"
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    return dirs

# -----------------------------
# Column mapping / 
# -----------------------------
EXPECTED_MAP = {
    "SNP": ["rsid", "hm_rsid", "SNP", "snp", "MarkerName", "rsids"],
    "A1": ["effect_allele", "hm_effect_allele", "assessed_allele", "A1", "Allele1", "alt"],
    "A2": ["other_allele", "hm_other_allele", "A2", "Allele2", "ref"],
    "BETA": ["beta", "hm_beta", "BETA", "Effect"],
    "SE": ["se", "standard_error", "hm_se", "SE", "StdErr", "sebeta"],
    "P": ["p_value", "p", "pval", "Pvalue", "P", "P_value"],
    "N": ["n", "N", "sample_size", "N_effective"],
    "EAF": ["eaf", "EAF", "frq", "Freq1", "EAF_HRC", "effect_allele_frequency", "af_alt"],
    "MAF": ["maf", "MAF", "MAF_ref"]
}

def read_header(input_path: Path) -> list[str]:
    open_func = gzip.open if str(input_path).endswith(".gz") else open
    with open_func(input_path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
    return header

def detect_mapping(header: list[str]) -> dict:
    hset = {h.strip(): h.strip() for h in header}
    mapping = {}
    missing = []
    for std, cands in EXPECTED_MAP.items():
        found = next((hset[cand] for cand in cands if cand in hset), None)
        if found:
            mapping[std] = found
        elif std in ["N", "EAF", "MAF", "SE"]:
            mapping[std] = None
        else:
            missing.append(std)
    return {"mapping": mapping, "missing": missing}

def write_colname_report(reports_dir: Path, clean_name: str, header: list[str], mapping: dict):
    out = reports_dir / f"colname_mapping_{clean_name}.tsv"
    with out.open("wt", encoding="utf-8") as w:
        w.write("original\tstandard\n")
        inv = {v: k for k, v in mapping.items() if v}
        for col in header:
            w.write(f"{col}\t{inv.get(col, '')}\n")
    return out

# -----------------------------
# Pandas Preprocessing / 
# -----------------------------
def preprocess_and_align(
    input_path: Path, snplist_path: Path, mapping: dict, 
    analysis_dir: Path, logger: logging.Logger, maf_threshold: float = 0.01
) -> Path:
    """
    EN: Filter MAF, remove palindromic SNPs, and strictly align alleles to reference.
    ZH: MAF, SNP, . 
    """
    logger.info("Starting strict alignment and filtering...")
    use_cols = {v: k for k, v in mapping.items() if v}
    comp = 'gzip' if str(input_path).endswith('.gz') else None
    
    df = pd.read_csv(input_path, sep='\t', usecols=use_cols.keys(), compression=comp)
    df.rename(columns=use_cols, inplace=True)
    
    initial_n = len(df)
    logger.info(f"Loaded {initial_n} SNPs.")
    df.dropna(subset=["P"], inplace=True)
    
    # MAF Filter
    if "MAF" in df.columns:
        df = df[df["MAF"] >= maf_threshold]
    elif "EAF" in df.columns:
        df["MAF_calc"] = np.where(df["EAF"] > 0.5, 1 - df["EAF"], df["EAF"])
        df = df[df["MAF_calc"] >= maf_threshold]
        df.drop(columns=["MAF_calc"], inplace=True)
    logger.info(f"After MAF >= {maf_threshold} filter: {len(df)} SNPs")

    # Palindromic Filter
    df["A1"], df["A2"] = df["A1"].str.upper(), df["A2"].str.upper()
    palindromic = (
        ((df["A1"] == "A") & (df["A2"] == "T")) | ((df["A1"] == "T") & (df["A2"] == "A")) |
        ((df["A1"] == "C") & (df["A2"] == "G")) | ((df["A1"] == "G") & (df["A2"] == "C"))
    )
    df = df[~palindromic]
    logger.info(f"Removed palindromic SNPs. Remaining: {len(df)} SNPs")

    # Reference Alignment
    ref_df = pd.read_csv(snplist_path, sep='\t', usecols=["SNP", "A1", "A2"])
    ref_df.rename(columns={"A1": "Ref_A1", "A2": "Ref_A2"}, inplace=True)
    merged = pd.merge(df, ref_df, on="SNP", how="inner")
    
    mask_match = (merged["A1"] == merged["Ref_A1"]) & (merged["A2"] == merged["Ref_A2"])
    mask_flip = (merged["A1"] == merged["Ref_A2"]) & (merged["A2"] == merged["Ref_A1"])
    merged = merged[mask_match | mask_flip].copy()
    
    flip_indices = (merged["A1"] == merged["Ref_A2"]) & (merged["A2"] == merged["Ref_A1"])
    logger.info(f"Flipping {flip_indices.sum()} SNPs (aligning alleles and adjusting Beta/EAF).")
    
    if "BETA" in merged.columns:
        merged.loc[flip_indices, "BETA"] = -merged.loc[flip_indices, "BETA"]
    if "EAF" in merged.columns:
        merged.loc[flip_indices, "EAF"] = 1 - merged.loc[flip_indices, "EAF"]
        
    merged.loc[flip_indices, "A1"] = merged.loc[flip_indices, "Ref_A1"]
    merged.loc[flip_indices, "A2"] = merged.loc[flip_indices, "Ref_A2"]
    merged.drop(columns=["Ref_A1", "Ref_A2"], inplace=True)

    if "EAF" in merged.columns and "MAF" in merged.columns:
        merged.drop(columns=["MAF"], inplace=True)
    
    out_path = analysis_dir / "preprocessed_input.tsv"
    merged.to_csv(out_path, sep='\t', index=False)
    return out_path

# -----------------------------
# Munge integration /  LDSC munge
# -----------------------------
def run_munge_sumstats(
    munge_py: Path, input_path: Path, snplist_path: Path,
    analysis_dir: Path, sample_N: int | None, clean_name: str,
    env: dict, log: logging.Logger,
) -> Path:
    out_prefix = analysis_dir / f"{clean_name}.munged"
    cmd = [
        sys.executable, str(munge_py),
        "--sumstats", str(input_path),
        "--out", str(out_prefix),
        "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
        "--signed-sumstats", "BETA,0",
        "--merge-alleles", str(snplist_path),
        "--chunksize", "500000",
    ]
    if sample_N:
        cmd.extend(["--N", str(sample_N)])
        
    log.info("Running munge_sumstats.py: %s", " ".join(cmd))
    subprocess.run(cmd, env=env, check=True)
    return Path(str(out_prefix) + ".sumstats.gz")

def setup_logger(logs_dir: Path, clean_name: str) -> logging.Logger:
    log_path = logs_dir / f"{clean_name}.log"
    logger = logging.getLogger(clean_name)
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh = logging.FileHandler(log_path, encoding='utf-8')
    fh.setFormatter(fmt)
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.handlers.clear()
    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger

# -----------------------------
# Main / 
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Unified GWAS → LAVA preprocessing template")
    parser.add_argument("--input", required=True, help="Path to input GWAS summary stats")
    parser.add_argument("--munge_py", required=True, help="Path to ldsc munge_sumstats.py")
    parser.add_argument("--snplist", required=True, help="Path to UKB LAVA LD snplist")
    parser.add_argument("--out_base", required=True, help="Output base directory")
    parser.add_argument("--phenotype", required=True, help="Phenotype ID for input.info")
    parser.add_argument("--N", type=int, help="Sample size (N)")
    parser.add_argument("--cases", type=int, help="Number of cases (for binary)")
    parser.add_argument("--controls", type=int, help="Number of controls (for binary)")
    parser.add_argument("--maf", type=float, default=0.01, help="MAF threshold (default: 0.01)")
    args = parser.parse_args()

    input_path = Path(args.input).resolve()
    munge_py = Path(args.munge_py).resolve()
    snplist_path = Path(args.snplist).resolve()
    out_base = Path(args.out_base).resolve()

    clean_name = f"{args.phenotype}_LAVA_Input"
    dirs = init_output_structure(out_base, clean_name)
    logger = setup_logger(dirs["logs"], clean_name)
    env = sanitize_env()

    # Checking input dependencies
    for p in (input_path, munge_py, snplist_path):
        if not p.exists():
            logger.error("Path not found: %s", p)
            sys.exit(1)

    # 1. Detect mapping
    header = read_header(input_path)
    det = detect_mapping(header)
    if det["missing"]:
        logger.error("Missing required standard columns: %s", det["missing"])
        sys.exit(2)
    write_colname_report(dirs["reports"], clean_name, header, det["mapping"])

    # 2. Strict Preprocessing & Munge
    cleaned_input_path = preprocess_and_align(
        input_path, snplist_path, det["mapping"], dirs["analysis"], logger, args.maf
    )
    
    use_n_col = bool(det["mapping"].get("N"))
    munged_path = run_munge_sumstats(
        munge_py, cleaned_input_path, snplist_path, dirs["analysis"],
        args.N if not use_n_col else None, clean_name, env, logger
    )
    
    # 3. Write input.info.txt
    info_out = dirs["analysis"] / "input.info.txt"
    with info_out.open("wt", encoding="utf-8") as w:
        w.write("phenotype\tcases\tcontrols\tprevalence\tfilename\n")
        w.write(f"{args.phenotype}\t{args.cases or 'NA'}\t{args.controls or 'NA'}\t\t{munged_path}\n")

    logger.info("All done. Standardized LAVA input generated successfully.")

if __name__ == "__main__":
    main()
