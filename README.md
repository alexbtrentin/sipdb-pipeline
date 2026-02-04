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
✅ Checkpoint/resume functionality for long-running analyses  
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
  out_dir = "~/path/to/your/output_folder",  # <- Set your output directory here
  heavy_threshold = 1.700,           # Relaxed: 1.700 g/ml
  run_dual_threshold = TRUE,         # Run both relaxed AND stringent
  stringent_threshold = 1.725,       # Stringent: 1.725 g/ml
  n_cores = 4,                       # Adjust based on your system
  run_id = "my_analysis_v1"          # Name your run for checkpointing
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

## 🔄 Checkpoint/Resume Functionality

The pipeline automatically saves checkpoints after each major step. If an analysis fails (memory issues, cluster timeout, etc.), you can resume from where it left off.

### How It Works

Checkpoints are saved to `{out_dir}/checkpoints/` with the pattern:
```
{run_id}_step01_compatibility.rds
{run_id}_step03_relaxed.rds
{run_id}_step04_master_relaxed.rds
{run_id}_step05_stringent.rds
{run_id}_step06_master_stringent.rds
{run_id}_step07_consensus.rds
{run_id}_step08_summaries.rds
```

### First Run (Give It a Memorable Name)
```r
results <- run_complete_sip_analysis(
  studies = studies,
  out_dir = main_dir,
  run_id = "sipdb_analysis_v1",  # Name your run!
  n_cores = 4
)
```

### Check Progress (If Run Failed)
```r
# See which steps completed
check_pipeline_progress(run_id = "sipdb_analysis_v1")

# Output:
# ════════════════════════════════════════════════════════════
# CHECKPOINT STATUS
# ════════════════════════════════════════════════════════════
# 
# 🔖 Run ID: sipdb_analysis_v1
# ----------------------------------------
#   ✅ step01: Compatibility check (2025-01-15 14:32)
#   ✅ step03: Relaxed analysis (2025-01-15 14:45)
#   ✅ step04: Master relaxed (2025-01-15 14:46)
#   ⬜ step05: Stringent analysis
#   ⬜ step06: Master stringent
#   ⬜ step07: Consensus
#   ⬜ step08: Summaries
```

### Resume from Last Checkpoint
```r
# Resume with the SAME run_id
results <- run_complete_sip_analysis(
  studies = studies,
  out_dir = main_dir,
  resume = TRUE,                   # ← Enable resume mode
  run_id = "sipdb_analysis_v1",    # ← Same run_id as before
  n_cores = 4
)

# Pipeline will skip completed steps and continue from where it stopped
```

### Checkpoint Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `resume` | `FALSE` | Set `TRUE` to load from existing checkpoints |
| `run_id` | Auto-generated timestamp | Unique identifier for the run |
| `checkpoint_dir` | `{out_dir}/checkpoints` | Where checkpoint files are stored |

---

## 📊 Output Files

All results are saved to `{out_dir}/`:
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
├─ checkpoints/                              # Checkpoint files for resume
│   └─ {run_id}_step*.rds
├─ debug_logs/                               # Analysis logs per study
├─ qc_reports/                               # Quality control reports
└─ plots/                                    # Generated visualizations
```

---

## 🔧 Advanced Usage

### Full Configuration (All Parameters)
```r
results <- run_complete_sip_analysis(
  # Required
  studies = studies,
  
  # Output settings
  output_prefix = "sip_analysis_v1.0",
  out_dir = "~/output/my_analysis",
  
  # Density thresholds
  heavy_threshold = 1.700,
  stringent_threshold = 1.725,
  run_dual_threshold = TRUE,
  
  # Performance
  n_cores = 4,
  
  # Checkpoint/resume
  resume = FALSE,
  run_id = "my_analysis_v1",
  checkpoint_dir = NULL  # Uses out_dir/checkpoints by default
)
```

### Single Threshold Analysis (Faster)
```r
# Run only relaxed threshold
results <- run_complete_sip_analysis(
  studies = studies,
  heavy_threshold = 1.700,
  run_dual_threshold = FALSE,  # Skip stringent analysis
  n_cores = 8,
  run_id = "quick_test"
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
  n_cores = 4,
  run_id = "custom_thresholds"
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
- **Checkpoint/resume**: Long analyses can be resumed after failures
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

*Manuscript in preparation — citation will be updated upon publication.*

GitHub: https://github.com/alexbtrentin/sipdb-pipeline  
OSF: https://osf.io/cwg6r/

---

## 🪪 License

MIT License — free to use, modify, and adapt with attribution.

---

## 🐛 Troubleshooting

### Common Issues

**"Object 'main_dir' not found"**
- Set `main_dir` BEFORE sourcing the pipeline:
```r
  main_dir <- "~/output/my_analysis"
  source("R/sip_meta_analysis_pipeline_v1.0.R")
```
- Or always pass `out_dir` explicitly in the function call

**"No valid gradient data"**
- Ensure `gradient_pos_density` or `gradient_position` columns exist
- Check density values are in biological range (1.650-1.800 g/ml)

**"Insufficient samples"**
- Minimum 2 labeled and 2 unlabeled samples required per comparison
- Check sample identification with `identify_labeled_unlabeled_samples()`

**Analysis failed mid-run**
- Use checkpoint/resume functionality:
```r
  results <- run_complete_sip_analysis(
    studies = studies,
    resume = TRUE,
    run_id = "your_original_run_id"
  )
```

**ALDEx2 skipped**
- Normal for large datasets (>3000 taxa or >50 samples)
- Automatic decision based on computational complexity

**Memory issues**
- Reduce `n_cores` to lower memory footprint
- Use checkpoint/resume to break analysis into chunks
- Run studies individually and combine later
- Use `clear_cache()` between analyses

**HPC/cluster timeout**
- Use `run_id` to name your run
- If job times out, resubmit with `resume = TRUE`
- Checkpoints persist across sessions

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

**Version:** 1.0.1  
**Last Updated:** February 2026
