# 🧬 pathWAS: Analysis Framework for the Pathway-Wide Association Study

**pathWAS** is an open-source Python package for conducting **pathway-level analogs of transcriptome-wide association studies (TWAS)**. It enables researchers to model genetically regulated pathway activity and test its association with complex traits using both **individual-level** and **summary-level** data.

---

## 🔍 Overview

Traditional TWAS focuses on associating genetically predicted **gene expression** with traits of interest. `pathWAS` generalizes this approach to the **pathway level** by:

- Modeling **SNP-to-pathway activation scores (PAS)** using gene expression and eQTL data
- Supporting multiple methods for computing PAS, including a proprietary **activity-weighted sum** method
- Testing for **trait associations** using both individual-level and summary-level data
- Enabling **genetic correlation** testing between pathway activity and phenotypes via LD Score Regression

---

## ⚙️ Features

- ✅ PAS computation using mean, median, sum, or activity-weighted methods  
- ✅ Support for KEGG, GO, MSigDB, and custom gene sets  
- ✅ SNP-to-PAS model building using elastic net or other methods  
- ✅ Individual-level or summary-based association tests (TWAS-style Z-scores)  
- ✅ Genetic correlation estimation via LDSC  
- ✅ Modular, testable architecture  

---

## 📦 Installation

```bash
git clone https://github.com/yourusername/pathWAS.git
cd pathWAS
pip install -e .
```

## 📁 Project Structure

pathWAS/
├── pathwas/                # Core package
│   ├── pas/                # PAS computation methods
│   ├── pathways/           # Gene set loading (KEGG, GO, etc.)
│   ├── model/              # SNP → PAS modeling
│   ├── association/        # Trait association & genetic correlation
│   ├── io/                 # Input/output utilities
│   └── cli.py              # Command-line interface (in development)
├── data/                   # Reference files (pathways, LD, annotations)
├── models/                 # Trained SNP→PAS models
├── scripts/                # Example analysis scripts
├── tests/                  # Unit tests
├── notebooks/              # Development notebooks
└── examples/               # Example runs

## 🧪 Example Usage

### Step 1: Compute PAS

```python  
from pathwas.pas.compute_pas import compute_pas

pas_matrix = compute_pas(expression_df, pathway_definitions, method="activity_weighted")
```

### Step 2: Build SNP → PAS Models

```python  
from pathwas.model.train_model import train_elastic_net_model

model = train_elastic_net_model(genotype_df, pas_matrix["Oxidative_Phosphorylation"])
```

### Step 3: Perform Pathway-TWAS (summary-level)

```python  
from pathwas.association.pathway_test import test_pathway_twas

z, var_pas = test_pathway_twas(snp_weights, gwas_zscores, ld_matrix)
```

### Step 4: Estimate Genetic Correlation

```python  
from pathwas.association.genetic_correlation 

import run_ldsc_rg

run_ldsc_rg("pas.sumstats.gz", "trait.sumstats.gz", "ref_ld_chr/", "w_ld_chr/", "output/pathway_trait_rg")
``` 

## 📚 Documentation

📖 Full documentation and tutorials (coming soon).

## 🤝 Contributing

We welcome contributions! To get started:

1. Fork the repo

2. Create a feature branch (git checkout -b my-feature)

3. Commit and push your changes

4. Open a pull request

## 📄 License

Distributed under the MIT License. See LICENSE for details.

## 👥 Authors

Yousef Mustafa, MS, PhD (lead developer)

Contributions welcome!

## 🧠 Acknowledgments

This project is inspired by methodologies from:

PrediXcan / S-PrediXcan

MAGMA / GSEA

WGCNA

LDSC
