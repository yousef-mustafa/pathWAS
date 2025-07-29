## ------------------------------------------------------------------------------------------- ##
## VCF Processing Utilities                                                                   ##
## ------------------------------------------------------------------------------------------- ##
## @script: vcf_processing.py                                                                 ##
##                                                                                             ##
## @description: Wrapper functions for running PLINK to perform LD pruning of variant files.   ##
##                                                                                             ##
## @author: Yousef Mustafa, Lab of Dr. William Bush.                                           ##
## ------------------------------------------------------------------------------------------- ##

import subprocess
import logging
from pathlib import Path
from typing import Optional


def ld_prune(vcf: str, plink_path: str = 'plink', out_dir: str = 'ld_pruned', r2_threshold: float = 0.1, window: int = 50, step: int = 5) -> Path:
    """Run PLINK to perform LD pruning on a VCF file.

    Parameters
    ----------
    vcf : str
        Path to input VCF file.
    plink_path : str, default 'plink'
        Path to PLINK executable.
    out_dir : str, default 'ld_pruned'
        Directory to store pruned output files.
    r2_threshold : float, default 0.1
        Threshold for pairwise LD (r^2) pruning.
    window : int, default 50
        Window size in kb.
    step : int, default 5
        Step size in variants.

    Returns
    -------
    Path
        Path to the pruned VCF file.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logging.info("LD pruning %s", vcf)
    out_prefix = out_dir / 'pruned'
    cmd = [
        plink_path,
        '--vcf', vcf,
        '--indep-pairwise', str(window), str(step), str(r2_threshold),
        '--out', str(out_prefix)
    ]
    subprocess.run(cmd, check=True)
    logging.debug("First PLINK command finished")
    cmd2 = [
        plink_path,
        '--vcf', vcf,
        '--extract', f'{out_prefix}.prune.in',
        '--recode', 'vcf',
        '--out', str(out_prefix)
    ]
    subprocess.run(cmd2, check=True)
    pruned_vcf = out_prefix.with_suffix('.vcf')
    logging.info("LD pruning complete -> %s", pruned_vcf)
    return pruned_vcf