## ------------------------------------------------------------------------------------------- ##
## Pathway Association Testing Utilities                                                       ##
## ------------------------------------------------------------------------------------------- ##
## @script: pathway_test.py                                                                    ##
##                                                                                             ##
## @description: Aggregate variants by pathway and perform simple correlation-based            ##
##               association tests with pathway activation scores.                             ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

from typing import Dict, Iterable, Optional, Tuple, Union

import numpy as np
import pandas as pd
import logging
from scipy import stats


def aggregate_variants(genotypes: pd.DataFrame, pathways: Dict[str, Iterable[str]]) -> Dict[str, pd.DataFrame]:
    """Group variant genotypes by pathway based on gene membership."""
    grouped = {}
    for pw, genes in pathways.items():
        cols = [c for c in genotypes.columns if c.split(':')[0] in genes]
        if cols:
            grouped[pw] = genotypes[cols]
            logging.debug("Pathway %s has %d variants", pw, len(cols))
    return grouped


def association_test(pas: pd.DataFrame, variants: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Simple association test between pathway activation scores and variant sets.

    For each pathway, compute correlation between PAS and each variant genotype, then
    record minimal p-value.
    """
    results = []
    for pw, geno in variants.items():
        if pw not in pas.columns:
            continue
        scores = pas[pw]
        pvals = []
        for col in geno.columns:
            r, p = stats.pearsonr(scores, geno[col])
            pvals.append(p)
        if pvals:
            min_p = min(pvals)
            logging.debug("%s min p-value %f", pw, min_p)
            results.append({'pathway': pw, 'min_p': min_p})
    return pd.DataFrame(results)


def test_pathway_twas(
    snp_weights: pd.Series,
    z_scores: pd.Series,
    ld_matrix: np.ndarray,
    snp_var: Optional[np.ndarray] = None,
    return_contributions: bool = False,
) -> Union[Tuple[float, float], Tuple[float, float, pd.Series]]:
    """Summary-level Pathway-TWAS association test.

    Parameters
    ----------
    snp_weights : pandas.Series
        SNP-to-pathway weights indexed by SNP identifier.
    z_scores : pandas.Series
        GWAS Z-scores indexed by SNP identifier.
    ld_matrix : numpy.ndarray
        Linkage disequilibrium (LD) matrix matching the order of ``snp_weights``.
    snp_var : numpy.ndarray, optional
        Variance of each SNP. If ``None``, the diagonal of ``ld_matrix`` is used.
    return_contributions : bool, optional
        If ``True``, also return per-SNP contribution terms.

    Returns
    -------
    tuple
        ``(Z_pathway, Var_PAS)`` or ``(Z_pathway, Var_PAS, contributions)`` when
        ``return_contributions`` is ``True``.
    """

    # Ensure matching indices
    if not snp_weights.index.equals(z_scores.index):
        try:
            z_scores = z_scores.loc[snp_weights.index]
        except KeyError as exc:  # pragma: no cover - defensive
            raise ValueError("Indices of snp_weights and z_scores must match") from exc
    logging.debug("Running pathway TWAS with %d SNPs", len(snp_weights))

    weights = snp_weights.astype(float).values
    z = z_scores.astype(float).values

    if ld_matrix.shape[0] != ld_matrix.shape[1] or ld_matrix.shape[0] != len(weights):
        raise ValueError("ld_matrix dimensions must match number of SNPs")

    if snp_var is None:
        snp_var = np.diag(ld_matrix)
    snp_var = np.asarray(snp_var, dtype=float)
    if snp_var.shape[0] != len(weights):
        raise ValueError("snp_var length must match number of SNPs")

    var_pas = float(weights.T @ ld_matrix @ weights)
    logging.debug("Var(PAS) = %f", var_pas)
    if var_pas <= 0:
        logging.warning("Non-positive PAS variance")
        if return_contributions:
            contrib = pd.Series(0.0, index=snp_weights.index)
            return np.nan, 0.0, contrib
        return np.nan, 0.0

    contrib_values = weights * (var_pas / snp_var) * z
    z_pathway = float(contrib_values.sum())

    if return_contributions:
        contrib = pd.Series(contrib_values, index=snp_weights.index)
        logging.info("Computed pathway Z-score %.3f", z_pathway)
        return z_pathway, var_pas, contrib
    logging.info("Computed pathway Z-score %.3f", z_pathway)
    return z_pathway, var_pas


# prevent pytest from collecting this function as a test
test_pathway_twas.__test__ = False
