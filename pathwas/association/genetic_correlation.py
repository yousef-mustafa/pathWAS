## ------------------------------------------------------------------------------------------- ##
## Genetic Correlation Estimation                                                              ##
## ------------------------------------------------------------------------------------------- ##
## @script: genetic_correlation.py                                                             ##
##                                                                                             ##
## @description: Utilities for generating PAS summary statistics and running LDSC              ##
##               regression to estimate genetic correlation between PAS and a trait.           ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

"""Utilities for estimating genetic correlation using LDSC."""

from __future__ import annotations

from typing import Optional
import pandas as pd
import numpy as np
import gzip
import subprocess
import os
import logging


# --------------------------- Summary statistics generation ---------------------------

def generate_pas_summary_stats(
    snp_weights: pd.Series,
    snp_var: Optional[pd.Series] = None,
    output_path: str = "pas_pathway.sumstats.gz",
    trait_name: str = "PAS",
    ld_score_path: Optional[str] = None,
) -> pd.DataFrame:
    """Create LDSC-compatible summary statistics for a PAS.

    Parameters
    ----------
    snp_weights : pandas.Series
        SNP-to-PAS weights indexed by SNP identifier.
    snp_var : pandas.Series, optional
        Variance of the SNPs. If ``None``, variance of 1 is assumed.
    output_path : str, optional
        File path to write ``.sumstats.gz`` output.
    trait_name : str, optional
        Name for the PAS trait.
    ld_score_path : str, optional
        Optional path to LD scores used to compute ``chi^2`` statistics.

    Returns
    -------
    pandas.DataFrame
        DataFrame compatible with LDSC ``--sumstats`` input.
    """

    if snp_var is None:
        snp_var = pd.Series(1.0, index=snp_weights.index)
    else:
        snp_var = snp_var.reindex(snp_weights.index, fill_value=1.0)

    if (snp_var <= 0).any():
        logging.warning("Some SNP variances are zero or negative")

    z_scores = snp_weights * np.sqrt(snp_var.astype(float))

    df = pd.DataFrame({
        "SNP": snp_weights.index,
        "A1": "A",
        "A2": "G",
        "Z": z_scores.astype(float),
        "N": 10000,
        "CHR": 1,
        "BP": range(1, len(snp_weights) + 1),
    })

    if ld_score_path is not None:
        try:
            ld_scores = pd.read_csv(ld_score_path, sep="\t")
            df = df.merge(ld_scores[["SNP", "L2"]], on="SNP", how="left")
            df["CHI2"] = df["Z"] ** 2
        except Exception as exc:  # pragma: no cover - runtime I/O errors
            logging.warning("Could not merge LD scores: %s", exc)

    with gzip.open(output_path, "wt") as gz:
        df.to_csv(gz, sep="\t", index=False)
    logging.info("Wrote PAS summary stats to %s", output_path)
    return df


# --------------------------- Run LDSC regression ---------------------------

def run_ldsc_rg(
    pas_sumstats_path: str,
    trait_sumstats_path: str,
    ref_ld_chr: str,
    w_ld_chr: str,
    output_prefix: str,
) -> None:
    """Run LDSC to estimate genetic correlation between PAS and trait.

    Parameters
    ----------
    pas_sumstats_path : str
        Path to PAS summary statistics produced by :func:`generate_pas_summary_stats`.
    trait_sumstats_path : str
        Path to GWAS summary statistics for the trait of interest.
    ref_ld_chr : str
        Directory containing reference LD scores per chromosome.
    w_ld_chr : str
        Directory containing weight LD scores per chromosome.
    output_prefix : str
        Prefix for LDSC output files.
    """

    for f in (pas_sumstats_path, trait_sumstats_path):
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file not found: {f}")

    cmd = [
        "ldsc.py",
        "--rg",
        f"{pas_sumstats_path},{trait_sumstats_path}",
        "--ref-ld-chr",
        ref_ld_chr,
        "--w-ld-chr",
        w_ld_chr,
        "--out",
        output_prefix,
    ]

    logging.info("Running LDSC: %s", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)

    if proc.returncode != 0:
        logging.error(proc.stderr)
        raise RuntimeError(f"LDSC failed with code {proc.returncode}")

    logging.info(proc.stdout)
    logging.info("LDSC finished; output prefix %s", output_prefix)