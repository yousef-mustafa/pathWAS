import os
import pandas as pd
import numpy as np
import subprocess

from pathwas.association.genetic_correlation import generate_pas_summary_stats, run_ldsc_rg


def test_generate_pas_summary_stats(tmp_path):
    weights = pd.Series([0.2, -0.3], index=["rs1", "rs2"])
    vars_ = pd.Series([1.0, 0.5], index=["rs1", "rs2"])
    out = tmp_path / "pas.sumstats.gz"
    df = generate_pas_summary_stats(weights, vars_, output_path=str(out), trait_name="PAS_test")
    assert out.exists()
    assert list(df.columns[:7]) == ["SNP", "A1", "A2", "Z", "N", "CHR", "BP"]
    assert len(df) == 2


def test_run_ldsc_rg(monkeypatch, tmp_path):
    calls = {}

    def fake_run(cmd, capture_output=True, text=True):
        calls["cmd"] = cmd
        class R:
            returncode = 0
            stdout = "ok"
            stderr = ""
        return R()

    monkeypatch.setattr(subprocess, "run", fake_run)
    p1 = tmp_path / "pas.sumstats.gz"
    p2 = tmp_path / "trait.sumstats.gz"
    p1.write_text("test")
    p2.write_text("test")

    run_ldsc_rg(str(p1), str(p2), "ref", "weights", str(tmp_path / "out"))
    assert "--rg" in calls["cmd"]