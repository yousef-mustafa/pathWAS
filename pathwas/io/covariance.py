## ------------------------------------------------------------------------------------------- ##
## Variant Covariance Utilities                                                                ##
## ------------------------------------------------------------------------------------------- ##
## @script: covariance.py                                                                      ##
##                                                                                             ##
## @description: Functions to load VCFs as genotype matrices and compute variant covariance.   ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

import numpy as np
import pandas as pd
import logging
from pathlib import Path


def compute_covariance(genotype_matrix: pd.DataFrame) -> pd.DataFrame:
    """Compute covariance matrix between variants.

    Parameters
    ----------
    genotype_matrix : pandas.DataFrame
        DataFrame with samples in rows and variants in columns (0/1/2 encoded).

    Returns
    -------
    pandas.DataFrame
        Covariance matrix of variants.
    """
    logging.debug("Computing covariance matrix for %d variants", genotype_matrix.shape[1])
    return genotype_matrix.cov()


def load_vcf_as_matrix(vcf_path: str) -> pd.DataFrame:
    """Load VCF into a genotype matrix using pandas.

    This function expects the VCF to be small enough to fit in memory and only
    parses genotype calls (assuming 0/1/2 encoding)."""
    logging.info("Loading VCF %s", vcf_path)
    import allel
    callset = allel.read_vcf(vcf_path, fields=['samples', 'calldata/GT'])
    gt = allel.GenotypeArray(callset['calldata/GT']).to_n_alt().T
    df = pd.DataFrame(gt, columns=[f'var_{i}' for i in range(gt.shape[1])])
    df.index = callset['samples']
    logging.info("Loaded genotype matrix with %d variants", df.shape[1])
    return df