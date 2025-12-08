# SIPdb v1.0 – Stable Isotope Probing Meta-Analysis Pipeline

Welcome to the official pipeline used to create **SIPdb v1.0**, a 22-study, cross-environment Stable Isotope Probing (SIP) database.  
This repository includes the **analysis and database building scripts**, while the **large versioned outputs** are hosted on OSF.

---

## 📁 Repository Structure

```
sipdb-pipeline/
├─ R/
│   ├─ sip_meta_analysis_pipeline_v1.0.R
│   └─ sip_database_builder_v1.0.R
├─ README.md
└─ (Large data stored on OSF, see link below)
```

---

## 🚀 What This Repository Provides

🔥 Reproducible SIP meta-analysis pipeline  
🔥 SQLite database construction script  
🔥 Cross-method consensus of incorporators (DESeq2, ALDEx2, limma, edgeR)  
🔥 Applicable to both new and legacy SIP experiments  

---

## 📦 OSF Repository (Large Files)

All large outputs (>100 MB) required to build **SIPdb v1.0** are hosted on OSF:

🔗 **https://osf.io/cwg6r/**

OSF contains:
- SQLite database (`sip_integrated_v1.sqlite`)
- Relaxed/strict master tables (`sip_analysis_v1.0_*.csv`)
- Unified multi-study R object (`studies_sipdb_v1.0.rds`)
- ASV FASTA archive (`all_sequences_v1.fasta`)
- Full documentation (Wiki + README)

---

## 💻 Quick Start

### 1) Load the study inputs
```r
studies <- readRDS("studies_sipdb_v1.0.rds")
length(studies)
names(studies)[1:5]
```

### 2) Run the SIP statistical meta-analysis
```r
source("sip_meta_analysis_pipeline_v1.0.R")
sip_results <- run_sip_meta_analysis(
  studies               = studies,
  cores                 = 8,
  return_summary_tables = TRUE,
  relaxed_threshold     = TRUE,
  stringent_threshold   = TRUE
)
```

### 3) Build the SIPdb SQLite database
```r
source("sip_database_builder_v1.0.R")
build_res <- build_optimized_sip_database_v7(
  studies            = studies,
  relax_tbl          = sip_results$relaxed,
  stringent_tbl      = sip_results$stringent
)
```

More detailed reproducibility instructions are available on OSF.

---

## 🧠 Why SIPdb Matters

- Creates **consensus incorporator profiles**
- Enables **cross-environment comparisons**
- Supports **meta-analysis-ready FASTA**
- Designed to support:
  - SIP
  - environmental microbiology
  - biodegradation assessments
  - soil/rhizosphere systems
  - wastewater bioreactors

---

## 📝 Citation

If you use SIPdb, please cite:

Trentin, A.B. (2025). *SIPdb v1.0: Stable Isotope Probing Database.*  
OSF Repository: https://osf.io/cwg6r/  
GitHub: https://github.com/alexbtrentin/sipdb-pipeline

---

## 🪪 License
MIT License — free to use, modify, and adapt with attribution.

---

## 📫 Contact

For questions, ideas, or contributions:  
a.trentinx@gmail.com
