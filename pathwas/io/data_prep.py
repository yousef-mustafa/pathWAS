## ------------------------------------------------------------------------------------------- ##
## Expression Data Preparation Utilities                                                      ##
## ------------------------------------------------------------------------------------------- ##
## @script: data_prep.py                                                                      ##
##                                                                                             ##
## @description: Convert gene identifiers, download MSigDB gene sets and otherwise prepare     ##
##               expression data for downstream analyses.                                      ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

"""Utilities for preparing gene expression data and pathway definitions."""
from typing import Iterable, Dict, List

import pandas as pd
import logging

try:
    import mygene
except ImportError:  # pragma: no cover - handled in runtime
    mygene = None

try:
    import gseapy as gp
except ImportError:  # pragma: no cover - handled in runtime
    gp = None


def _needs_conversion(genes: Iterable[str]) -> bool:
    """Return True if genes look like HGNC symbols or Entrez IDs."""
    for g in genes:
        if g.startswith("ENS"):
            continue
        # numeric -> likely Entrez
        if g.isdigit():
            return True
        # otherwise assume symbol
        if g.isalpha():
            return True
    return False


def convert_gene_ids(expr: pd.DataFrame, species: str = "human") -> pd.DataFrame:
    """Convert HGNC or Entrez gene IDs in ``expr`` to Ensembl IDs."""
    if mygene is None:
        raise ImportError("mygene is required for gene ID conversion")

    genes = list(expr.columns)
    if not _needs_conversion(genes):
        logging.debug("Gene IDs already look like Ensembl IDs")
        return expr

    mg = mygene.MyGeneInfo()
    res = mg.querymany(genes, scopes=["symbol", "entrezgene"], fields="ensembl.gene", species=species)
    mapping: Dict[str, str] = {}
    for r in res:
        if "notfound" in r and r["notfound"]:
            continue
        query = r.get("query")
        ens = r.get("ensembl")
        if isinstance(ens, list):
            ens = ens[0]
        if isinstance(ens, dict):
            ens = ens.get("gene")
        if isinstance(ens, str):
            mapping[query] = ens

    expr = expr.rename(columns=mapping)
    expr = expr.loc[:, ~expr.columns.duplicated()]
    logging.info("Converted %d gene identifiers", len(mapping))
    return expr


def convert_gene_list(genes: Iterable[str], species: str = "human") -> List[str]:
    """Convert a list of genes to Ensembl IDs."""
    if mygene is None:
        raise ImportError("mygene is required for gene ID conversion")

    genes = list(genes)
    if not _needs_conversion(genes):
        logging.debug("Gene list already in Ensembl format")
        return genes

    mg = mygene.MyGeneInfo()
    res = mg.querymany(genes, scopes=["symbol", "entrezgene"], fields="ensembl.gene", species=species)
    mapping: Dict[str, str] = {}
    for r in res:
        if "notfound" in r and r["notfound"]:
            continue
        query = r.get("query")
        ens = r.get("ensembl")
        if isinstance(ens, list):
            ens = ens[0]
        if isinstance(ens, dict):
            ens = ens.get("gene")
        if isinstance(ens, str):
            mapping[query] = ens
    logging.info("Converted %d gene identifiers in list", len(mapping))
    return [mapping.get(g, g) for g in genes]


def load_msigdb_library(library: str, organism: str = "Human") -> Dict[str, List[str]]:
    """Download a gene set library from MSigDB using ``gseapy``."""
    if gp is None:
        raise ImportError("gseapy is required to load MSigDB gene sets")

    logging.info("Downloading MSigDB library %s", library)
    return gp.get_library(library, organism=organism)