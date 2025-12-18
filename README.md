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

✅ Reproducible SIP meta-analysis pipeline  
✅ SQLite database construction script  
✅ Cross-method consensus of incorporators (DESeq2, ALDEx2, limma, edgeR)  
✅ Dual-threshold approach (relaxed + stringent)  
✅ Sliding window analysis for density gradients  
✅ Substrate-specific analysis support  
✅ Applicable to both new and legacy SIP experiments

---

## 📦 OSF Repository (Large Files)

All large outputs (>100 MB) required to build **SIPdb v1.0** are hosted on OSF:

🔗 **https://osf.io/cwg6r/**

OSF contains:
- SQLite database (`sip_integrated_v1.sqlite`)
- Relaxed/stringent master tables (`sip_analysis_v1.0_*.csv`)
- Unified multi-study R object (`studies_sipdb_v1.0.rds`)
- ASV FASTA archive (`all_sequences_v1.fasta`)
- Full documentation (Wiki + README)

---

## 💻 Quick Start

### 1) Load Required Packages
```r
# The pipeline loads these automatically, but ensure they're installed:
required_packages <- c(
  "phyloseq", "tidyverse", "data.table", "DESeq2", "edgeR", 
  "limma", "ALDEx2", "openxlsx", "parallel", "foreach", "doParallel"
)

# Install any missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

### 2) Load the Study Inputs
```r
# Download studies_sipdb_v1.0.rds from OSF first
studies <- readRDS("studies_sipdb_v1.0.rds")

# Verify loaded studies
cat(sprintf("Loaded %d studies:\n", length(studies)))
print(names(studies))
```

### 3) Source the Pipeline
```r
source("R/sip_meta_analysis_pipeline_v1.0.R")

# The pipeline will display initialization message when loaded:
# ████████████████████████████████████████████████████████████
# █     SIP META-ANALYSIS PIPELINE v1.0 - READY             █
# ████████████████████████████████████████████████████████████
```

### 4) Run the Complete SIP Analysis
```r
# Run with dual-threshold approach (recommended)
results <- run_complete_sip_analysis(
  studies = studies,
  output_prefix = "sip_analysis_v1.0",
  heavy_threshold = 1.700,         # Relaxed: 1.700 g/ml (default: 1.700)
  run_dual_threshold = TRUE,       # Run both relaxed AND stringent
  stringent_threshold = 1.725,     # Stringent: 1.725 g/ml (default: 1.725)
  n_cores = 4                      # Adjust based on your system
)

# Analysis will take 15-60 minutes depending on system and number of cores
```

### 5) Access Results
```r
# Master tables (all records)
master_relaxed <- results$master_relaxed      # Relaxed threshold results
master_stringent <- results$master_stringent  # Stringent threshold results

# Consensus incorporators (high-confidence taxa)
consensus_taxa <- results$consensus$incorporators
high_confidence <- results$consensus$consensus_taxa  # Found in BOTH thresholds
relaxed_only <- results$consensus$relaxed_only       # Only in relaxed

# Summary statistics
summaries <- results$summaries

# Files saved to disk
results$files  # Paths to all saved files

# View summary
cat(sprintf("\n📊 RESULTS SUMMARY:\n"))
cat(sprintf("Runtime: %.1f minutes\n", as.numeric(results$runtime)))
cat(sprintf("Relaxed results: %d records, %d significant (p<0.05)\n",
            nrow(master_relaxed), 
            sum(master_relaxed$padj < 0.05, na.rm = TRUE)))
if(!is.null(master_stringent)) {
  cat(sprintf("Stringent results: %d records, %d significant (p<0.05)\n",
              nrow(master_stringent),
              sum(master_stringent$padj < 0.05, na.rm = TRUE)))
}
cat(sprintf("Consensus incorporators: %d high-confidence, %d relaxed-only\n",
            results$consensus$summary$n_consensus,
            results$consensus$summary$n_relaxed_only))
```

### 6) Build the SQLite Database (Optional)
```r
source("R/sip_database_builder_v1.0.R")

# Build database from analysis results
db_results <- build_optimized_sip_database_v7(
  studies = studies,
  relax_tbl = results$master_relaxed,
  stringent_tbl = results$master_stringent,
  output_file = "sip_integrated_v1.sqlite"
)
```

---

## 📊 Output Files

All results are saved to `~/output/output_SIPdb_pipeline_v1.0/`:
```
output_SIPdb_pipeline_v1.0/
├─ tables/
│   ├─ sip_analysis_v1.0_YYYYMMDD.xlsx      # Excel report with summaries
│   └─ sip_analysis_summary_v1.0_*.txt      # Text summary
├─ raw_data/
│   ├─ sip_analysis_v1.0_relaxed_*.csv      # Relaxed results
│   ├─ sip_analysis_v1.0_stringent_*.csv    # Stringent results
│   └─ sip_analysis_v1.0_complete_*.RData   # Complete R workspace
├─ consensus_results/
│   └─ sip_analysis_v1.0_consensus_*.csv    # High-confidence incorporators
├─ debug_logs/                               # Analysis logs per study
├─ qc_reports/                               # Quality control reports
└─ plots/                                    # Generated visualizations
```

---

## 🔧 Advanced Usage

### Single Threshold Analysis (Faster)
```r
# Run only relaxed threshold
results <- run_complete_sip_analysis(
  studies = studies,
  heavy_threshold = 1.700,
  run_dual_threshold = FALSE,  # Skip stringent analysis
  n_cores = 8
)
```

### Custom Density Thresholds
```r
# Adjust for different isotopes or experimental conditions
results <- run_complete_sip_analysis(
  studies = studies,
  heavy_threshold = 1.695,      # More permissive
  stringent_threshold = 1.735,  # More conservative
  run_dual_threshold = TRUE,
  n_cores = 4
)
```

### Check Study Compatibility Before Running
```r
# Quick compatibility check
compatibility <- check_pipeline_compatibility(studies)

# View individual study overviews
for(study_name in names(studies)) {
  get_study_overview(studies[[study_name]], study_name)
}
```

### Monitor Memory Usage
```r
# Check memory footprint during analysis
report_memory_usage()

# Clear caches if needed
clear_cache()
```

---

## 🧠 Why SIPdb Matters

- **Multi-method consensus**: Combines DESeq2, edgeR, limma-voom, and ALDEx2
- **Dual-threshold validation**: High-confidence incorporators found at both thresholds
- **Cross-environment comparisons**: 22 studies across soil, rhizosphere, and wastewater
- **Sliding window analysis**: Robust identification across density gradients
- **Substrate-specific analysis**: Handles shared-control and matched-control designs
- **Meta-analysis ready**: Standardized taxonomy and FASTA sequences

**Applications:**
- Microbial stable isotope probing (SIP)
- Environmental microbiology
- Biodegradation assessments
- Soil/rhizosphere microbiome studies
- Wastewater bioreactor characterization
- Carbon and nitrogen cycling

---

## 📊 Key Analysis Features

### Statistical Methods
- **DESeq2**: Negative binomial modeling with size factor normalization
- **edgeR**: TMM normalization with quasi-likelihood testing
- **limma-voom**: Variance modeling with quality weights
- **ALDEx2**: Compositional data analysis (adaptive, complexity-aware)

### Meta-Analysis Components
- Fisher's method for p-value combination
- Inverse variance weighting for effect sizes
- Heterogeneity assessment across windows
- Direction consistency scoring

### Quality Control
- Automatic sample identification (MISIP-compliant)
- Library size and sparsity checks
- Study type detection (density/binary/multifraction)
- Comprehensive logging and diagnostics

---

## 📝 Citation

If you use SIPdb v1.0 or this pipeline, please cite:

**Trentin, A.B., Wilhelm, R.C., Gonzalez, J.M., & Shabtai, I.A. (2025).**  
*SIPdb v1.0: A comprehensive database for stable isotope probing meta-analysis.*  
GitHub: https://github.com/alexbtrentin/sipdb-pipeline  
OSF: https://osf.io/cwg6r/

---

## 🪪 License

MIT License — free to use, modify, and adapt with attribution.

---

## 🐛 Troubleshooting

### Common Issues

**"No valid gradient data"**
- Ensure `gradient_pos_density` or `gradient_position` columns exist
- Check density values are in biological range (1.650-1.800 g/ml)

**"Insufficient samples"**
- Minimum 2 labeled and 2 unlabeled samples required per comparison
- Check sample identification with `identify_labeled_unlabeled_samples()`

**"Object 'PIPELINE_VERSION' not found"**
- Add `PIPELINE_VERSION <- "1.0"` after `set.seed(42)` in the script (around line 46)

**ALDEx2 skipped**
- Normal for large datasets (>3000 taxa or >50 samples)
- Automatic decision based on computational complexity

**Memory issues**
- Reduce `n_cores` to lower memory footprint
- Run studies individually and combine later
- Use `clear_cache()` between analyses

---

## 📫 Contact

**Alex Batista Trentin**  
PhD Student, Agronomy (Soil Microbiology)  
Purdue University  

📧 batista2@purdue.edu  
📧 a.trentinx@gmail.com  
🔗 https://github.com/alexbtrentin

For questions, bug reports, or collaboration inquiries, please open an issue on GitHub or email directly.

---

**Version:** 1.0  
**Last Updated:** December 2025
