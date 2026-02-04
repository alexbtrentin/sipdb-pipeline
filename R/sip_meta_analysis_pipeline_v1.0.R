################################################################################
#' @title SIP Meta-Analysis Pipeline v1.0 - COMPREHENSIVE DOCUMENTED VERSION
#'
#' @description
#'   Complete Stable Isotope Probing (SIP) analysis pipeline with:
#'   - Multi-method differential abundance (DESeq2, edgeR, limma-voom, ALDEx2)
#'   - Sliding window analysis across density gradients
#'   - Meta-analytic consensus calling
#'   - Reproducibility scoring and confidence tiers
#'   - Dual threshold support (relaxed + stringent)
#'   - Substrate-specific analysis
#'   - Comprehensive outputs (Excel, CSV, plots)
#'
#' @version 1.0
#' @date 2025
#' @author SIP Pipeline Team
#'
#' @details
#'   This pipeline implements state-of-the-art methods for identifying
#'   isotope-incorporating taxa in SIP experiments. The modular design
#'   allows for flexible analysis with multiple statistical approaches
#'   combined through meta-analytic methods.
#'
#'   **Key Features:**
#'   - Automatic sample identification (labeled vs unlabeled)
#'   - Four complementary statistical methods
#'   - Sliding window analysis for gradient data
#'   - Meta-analytic combination with heterogeneity assessment
#'   - Reproducibility scoring across windows
#'   - Confidence tier classification
#'   - Parallel processing for performance
#'   - Comprehensive quality control
#'
#' @note This is version 1.0 - the production release
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(data.table)
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(ALDEx2)
  library(parallel)
  library(openxlsx)
  library(ggplot2)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(stringr)
  library(mgcv)
  library(MASS)
  library(Biostrings)
  library(pheatmap)
  library(vegan)
  library(foreach)
  library(doParallel)
})

options(stringsAsFactors = FALSE)
`%+%` <- function(a, b) paste0(a, b)

# Version check and reproducibility
PIPELINE_VERSION <- "1.0"
set.seed(42)  # For reproducibility

main_dir <- getwd()

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ensure_dir_for_file <- function(filepath) {
  ensure_dir(dirname(filepath))
  invisible(filepath)
}


# ============================================================================
# BIOLOGICAL PARAMETERS WITH JUSTIFICATION
# ============================================================================

BIOLOGICAL_PARAMS <- list(
  density = list(
    heavy_threshold = 1.725,      # CsCl density where 13C-DNA typically separates (Neufeld et al. 2007)
    light_threshold = 1.710,      # Upper bound for unlabeled DNA (Lueders et al. 2004)
    min_density = 1.650,          # Minimum biological CsCl density
    max_density = 1.800,          # Maximum biological CsCl density
    window_width = 0.020,         # Based on typical fraction resolution
    overlap = 0.010               # 50% overlap for smooth transitions
  ),
  
  filtering = list(
    min_prevalence = 0.001,       # 0.1% of average library size
    min_counts = 3,               # Minimum total counts across samples
    min_samples = 2,              # Minimum samples for detection
    max_sparsity = 0.995,         # Maximum allowed sparsity
    pseudocount = 0.5             # For numerical stability (Weiss et al. 2017)
  ),
  
  statistics = list(
    log2fc_threshold = 1.0,       # 2-fold change minimum
    padj_threshold = 0.05,        # Adjusted p-value threshold
    min_samples_per_group = 2     # Minimum for statistical comparison
  )
)

# Global taxonomy cache
.GlobalEnv$.sip_taxonomy_cache <- new.env()
.GlobalEnv$.sip_validation_cache <- new.env()

# ============================================================================
# ENHANCED LOGGING SYSTEM
# ============================================================================

#' Log Message
#'
#' @description [Function description]
#' @export
log_message <- function(message, log_file = NULL, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_message <- sprintf("[%s] %s: %s\n", timestamp, level, message)
  
  if (!is.null(log_file)) {
    cat(formatted_message, file = log_file, append = TRUE)
  }
  
  # Color coding for console
  if (level == "ERROR") {
    cat(crayon::red(formatted_message))
  } else if (level == "WARNING") {
    cat(crayon::yellow(formatted_message))
  } else if (level == "SUCCESS") {
    cat(crayon::green(formatted_message))
  } else {
    cat(formatted_message)
  }
}

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

#' Validate Density Values
#'
#' @description [Function description]
#' @export
validate_density_values <- function(densities, log_file = NULL) {
  if (is.null(densities) || length(densities) == 0) {
    log_message("No density values provided", log_file, "ERROR")
    return(FALSE)
  }
  
  valid <- !is.na(densities) &
    densities >= BIOLOGICAL_PARAMS$density$min_density &
    densities <= BIOLOGICAL_PARAMS$density$max_density
  
  if (sum(valid) == 0) {
    log_message("No valid density values in biological range", log_file, "ERROR")
    return(FALSE)
  }
  
  log_message(sprintf("Validated %d/%d density values (%.1f%% valid)",
                      sum(valid), length(densities),
                      100 * sum(valid)/length(densities)),
              log_file, "SUCCESS")
  return(TRUE)
}

#' Validate Sample Sizes
#'
#' @description [Function description]
#' @export
validate_sample_sizes <- function(labeled, unlabeled, min_required = 2, log_file = NULL) {
  n_labeled <- sum(labeled)
  n_unlabeled <- sum(unlabeled)
  
  if (n_labeled < min_required) {
    log_message(sprintf("Insufficient labeled samples: %d < %d", n_labeled, min_required),
                log_file, "ERROR")
    return(FALSE)
  }
  
  if (n_unlabeled < min_required) {
    log_message(sprintf("Insufficient unlabeled samples: %d < %d", n_unlabeled, min_required),
                log_file, "ERROR")
    return(FALSE)
  }
  
  log_message(sprintf("Sample sizes OK: %d labeled, %d unlabeled", n_labeled, n_unlabeled),
              log_file, "SUCCESS")
  return(TRUE)
}

# ============================================================================
# CORE HELPER FUNCTIONS
# ============================================================================

#' Standardize Results
#'
#' @description [Function description]
#' @export
standardize_results <- function(results) {
  if (is.null(results) || nrow(results) == 0) return(NULL)
  
  # Standardize column names
  results <- results %>%
    rename_with(~case_when(
      . == "FDR"       ~ "padj",
      . == "adj.P.Val" ~ "padj",
      . == "logFC"     ~ "log2FC",
      TRUE             ~ .
    )) %>%
    mutate(
      padj   = if("padj" %in% colnames(.)) as.numeric(padj) else NA_real_,
      log2FC = if("log2FC" %in% colnames(.)) as.numeric(log2FC) else NA_real_,
      pvalue = if("pvalue" %in% colnames(.)) as.numeric(pvalue) else NA_real_
    ) %>%
    filter(!is.na(log2FC), !is.na(pvalue))  # Remove invalid results
  
  return(results)
}

# Improved p-value combination with method selection
#' Combine P Values
#'
#' @description [Function description]
#' @export
combine_p_values <- function(p_values, method = "fisher", weights = NULL) {
  p_values <- p_values[!is.na(p_values) & is.finite(p_values)]
  p_values <- pmax(p_values, 1e-300)  # Avoid numerical issues
  
  if (length(p_values) == 0) return(1)
  if (length(p_values) == 1) return(p_values[1])
  
  if (method == "fisher") {
    stat <- -2 * sum(log(p_values))
    return(pchisq(stat, df = 2 * length(p_values), lower.tail = FALSE))
    
  } else if (method == "stouffer") {
    if (is.null(weights)) weights <- rep(1, length(p_values))
    z_scores <- qnorm(1 - p_values)
    combined_z <- sum(z_scores * weights) / sqrt(sum(weights^2))
    return(pnorm(combined_z, lower.tail = FALSE))
    
  } else if (method == "tippett") {
    return(1 - (1 - min(p_values))^length(p_values))
    
  } else {
    stop("Unknown p-value combination method")
  }
}

# Density value parsing with validation
#' Parse Density Values
#'
#' @description [Function description]
#' @export
parse_density_values <- function(density_strings) {
  densities <- sapply(density_strings, function(x) {
    if (is.na(x) || x == "" || x == "NA" || x == "UNK") return(NA_real_)
    if (is.numeric(x)) return(as.numeric(x))
    
    x <- as.character(x)
    
    # Handle ranges
    if (grepl("-", x) && !grepl("^-", x)) {
      parts <- strsplit(x, "-")[[1]]
      nums <- as.numeric(parts)
      if (length(nums) == 2 && all(!is.na(nums))) {
        return(mean(nums))
      }
    }
    
    return(as.numeric(gsub("[^0-9.-]", "", x)))
  })
  
  densities <- unname(densities)
  
  # Validate biological range
  invalid <- !is.na(densities) &
    (densities < BIOLOGICAL_PARAMS$density$min_density |
       densities > BIOLOGICAL_PARAMS$density$max_density)
  
  if (any(invalid)) {
    densities[invalid] <- NA_real_
  }
  
  return(densities)
}

# ============================================================================
# ENHANCED SAMPLE IDENTIFICATION
# ============================================================================

#' Identify Labeled Unlabeled Samples
#'
#' @description [Function description]
#' @export
identify_labeled_unlabeled_samples <- function(physeq, log_file = NULL) {
  meta <- as.data.frame(sample_data(physeq))
  
  log_message("========================================", log_file)
  log_message("SAMPLE IDENTIFICATION", log_file)
  log_message("========================================", log_file)
  
  # Initialize detection vectors
  unlabeled_samples <- rep(FALSE, nrow(meta))
  labeled_samples <- rep(FALSE, nrow(meta))
  
  # Track detection methods and confidence
  methods_detected <- character()
  detection_confidence <- "low"
  
  # Priority 1: MISIP Standard (highest confidence)
  if ("isotopolog_label" %in% colnames(meta)) {
    unlabeled_misip <- grepl("natural abundance|^12C$|^14N$|^16O$|^H2O$|unlabeled|unlabelled",
                             meta$isotopolog_label, ignore.case = TRUE)
    labeled_misip <- grepl("isotopically labelled|isotopically labeled|^13C$|^15N$|^18O$|^2H$|^D2O$|labeled|labelled",
                           meta$isotopolog_label, ignore.case = TRUE)
    
    if (any(unlabeled_misip) || any(labeled_misip)) {
      unlabeled_samples <- unlabeled_samples | unlabeled_misip
      labeled_samples <- labeled_samples | labeled_misip
      methods_detected <- c(methods_detected, "MISIP_standard")
      detection_confidence <- "high"
      log_message(sprintf("  ✓ MISIP standard: %d unlabeled, %d labeled",
                          sum(unlabeled_misip), sum(labeled_misip)), log_file, "SUCCESS")
    }
  }
  
  # Priority 2: isotopolog column
  if ("isotopolog" %in% colnames(meta) && detection_confidence != "high") {
    unlabeled_patterns <- c("none", "natural", "unlabeled", "unlabelled", "control",
                            "12C", "14N", "16O", "H2O")
    isotope_patterns <- c("13C", "15N", "18O", "2H", "D2O")
    
    unlabeled_regex <- paste0("(^|[^0-9])(", paste(unlabeled_patterns, collapse = "|"), ")($|[^0-9])")
    unlabeled_iso <- grepl(unlabeled_regex, meta$isotopolog, ignore.case = TRUE)
    
    has_isotope <- grepl(paste(isotope_patterns, collapse = "|"), meta$isotopolog, ignore.case = FALSE)
    not_natural <- !grepl("natural|none|unlabeled|unlabelled|control", meta$isotopolog, ignore.case = TRUE)
    labeled_iso <- has_isotope & not_natural
    
    if (any(unlabeled_iso) || any(labeled_iso)) {
      unlabeled_samples <- unlabeled_samples | unlabeled_iso
      labeled_samples <- labeled_samples | labeled_iso
      methods_detected <- c(methods_detected, "isotopolog")
      if (detection_confidence != "high") detection_confidence <- "medium"
      log_message(sprintf("  ✓ Isotopolog: %d unlabeled, %d labeled",
                          sum(unlabeled_iso), sum(labeled_iso)), log_file, "SUCCESS")
    }
  }
  
  # Conflict resolution
  both <- unlabeled_samples & labeled_samples
  if (any(both)) {
    log_message(sprintf("  ⚠ Resolving %d conflicting classifications", sum(both)),
                log_file, "WARNING")
    
    if ("isotopolog_label" %in% colnames(meta)) {
      unlabeled_samples[both] <- grepl("natural abundance", meta$isotopolog_label[both], ignore.case = TRUE)
      labeled_samples[both] <- !unlabeled_samples[both]
    } else {
      # Conservative: mark as neither
      unlabeled_samples[both] <- FALSE
      labeled_samples[both] <- FALSE
    }
  }
  
  # Validation
  n_unlabeled <- sum(unlabeled_samples)
  n_labeled <- sum(labeled_samples)
  
  log_message(sprintf("SUMMARY: %d unlabeled, %d labeled (confidence: %s)",
                      n_unlabeled, n_labeled, detection_confidence),
              log_file, if(n_unlabeled > 0 && n_labeled > 0) "SUCCESS" else "WARNING")
  
  # Extract isotope information
  isotopes <- character()
  if ("isotope" %in% colnames(meta)) {
    isotopes <- unique(meta$isotope[labeled_samples & !is.na(meta$isotope)])
  }
  
  # Extract substrate information
  substrates <- character()
  substrate_counts <- list()
  
  if ("isotopolog" %in% colnames(meta)) {
    labeled_isotopologs <- meta$isotopolog[labeled_samples]
    labeled_isotopologs <- labeled_isotopologs[!is.na(labeled_isotopologs) & labeled_isotopologs != ""]
    
    # Remove pure isotope entries
    pure_isotopes <- c("13C", "15N", "18O", "2H", "D2O")
    labeled_isotopologs <- labeled_isotopologs[!labeled_isotopologs %in% pure_isotopes]
    
    if (length(labeled_isotopologs) > 0) {
      substrate_table <- table(labeled_isotopologs)
      substrates <- names(substrate_table)
      substrate_counts <- as.list(substrate_table)
      
      log_message(sprintf("  Substrates detected: %s",
                          paste(substr(substrates, 1, 20), collapse = ", ")),
                  log_file)
    }
  }
  
  return(list(
    unlabeled = unlabeled_samples,
    labeled = labeled_samples,
    metadata = meta,
    isotopes = isotopes,
    substrates = substrates,
    substrate_counts = substrate_counts,
    n_unlabeled = n_unlabeled,
    n_labeled = n_labeled,
    methods_detected = methods_detected,
    detection_confidence = detection_confidence,
    has_both = (n_unlabeled > 0 && n_labeled > 0)
  ))
}

# ============================================================================
# ENHANCED STUDY TYPE DETECTION
# ============================================================================

#' Detect Study Type
#'
#' @description [Function description]
#' @export
detect_study_type <- function(physeq, sample_info, log_file = NULL) {
  meta <- sample_info$metadata
  study_type <- list()
  
  log_message("----------------------------------------", log_file)
  log_message("STUDY TYPE DETECTION", log_file)
  
  # Density analysis
  has_density <- FALSE
  density_quality <- "none"
  if ("gradient_pos_density" %in% colnames(meta)) {
    density_values <- parse_density_values(meta$gradient_pos_density)
    sip_mask <- sample_info$labeled | sample_info$unlabeled
    density_values_sip <- density_values[sip_mask]
    
    valid_densities <- !is.na(density_values_sip) &
      density_values_sip >= BIOLOGICAL_PARAMS$density$min_density &
      density_values_sip <= BIOLOGICAL_PARAMS$density$max_density
    
    n_valid <- sum(valid_densities)
    if (n_valid > 0) {
      has_density <- TRUE
      density_range <- max(density_values_sip[valid_densities]) -
        min(density_values_sip[valid_densities])
      
      # Quality assessment
      if (n_valid >= 20 && density_range >= 0.05) {
        density_quality <- "high"
      } else if (n_valid >= 10 && density_range >= 0.03) {
        density_quality <- "medium"
      } else {
        density_quality <- "low"
      }
      
      log_message(sprintf("  Density: %d values, range %.3f g/ml, quality: %s",
                          n_valid, density_range, density_quality),
                  log_file, "SUCCESS")
    }
  }
  
  # Position analysis
  has_positions <- FALSE
  n_positions <- 0
  position_quality <- "none"
  if ("gradient_position" %in% colnames(meta)) {
    sip_samples <- sample_info$labeled | sample_info$unlabeled
    positions <- meta$gradient_position[sip_samples]
    
    valid_positions <- positions[!is.na(positions) & positions > 0 & positions < 100]
    unique_positions <- unique(valid_positions)
    n_positions <- length(unique_positions)
    
    if (n_positions >= 1) {
      has_positions <- TRUE
      
      # Quality assessment
      samples_per_position <- table(valid_positions)
      min_samples <- min(samples_per_position)
      
      if (n_positions >= 5 && min_samples >= 3) {
        position_quality <- "high"
      } else if (n_positions >= 3 && min_samples >= 2) {
        position_quality <- "medium"
      } else {
        position_quality <- "low"
      }
      
      log_message(sprintf("  Positions: %d unique, quality: %s",
                          n_positions, position_quality),
                  log_file, "SUCCESS")
    }
  }
  
  # Decision logic based on quality
  if (has_density && density_quality != "none") {
    study_type$type <- "density"
    study_type$gradient_var <- "gradient_pos_density"
    study_type$quality <- density_quality
    log_message(sprintf("  → Study type: DENSITY (%s quality)", density_quality),
                log_file, "SUCCESS")
    
  } else if (has_positions) {
    if (n_positions <= 2) {
      study_type$type <- "binary"
      study_type$quality <- position_quality
      log_message("  → Study type: BINARY", log_file, "SUCCESS")
    } else {
      study_type$type <- "multifraction"
      study_type$quality <- position_quality
      log_message(sprintf("  → Study type: MULTIFRACTION (%d positions)", n_positions),
                  log_file, "SUCCESS")
    }
    study_type$gradient_var <- "gradient_position"
    
  } else {
    log_message("  → No valid gradient data - SKIPPING", log_file, "ERROR")
    return(NULL)
  }
  
  study_type$has_density <- has_density
  study_type$n_positions <- n_positions
  study_type$density_quality <- density_quality
  study_type$position_quality <- position_quality
  
  return(study_type)
}

# ============================================================================
# SUBSTRATE DESIGN DETECTION
# ============================================================================

#' Detect Substrate Design
#'
#' @description [Function description]
#' @export
detect_substrate_design <- function(physeq, sample_info, log_file = NULL) {
  meta <- sample_info$metadata
  design_info <- list()
  
  if (!"isotopolog" %in% colnames(meta)) {
    design_info$type <- "no_substrates"
    design_info$substrates <- character()
    design_info$control_type <- NA
    design_info$n_substrates <- 0
    return(design_info)
  }
  
  control_isotopologs <- unique(meta$isotopolog[sample_info$unlabeled])
  control_isotopologs <- control_isotopologs[!is.na(control_isotopologs)]
  
  labeled_substrates <- sample_info$substrates
  n_substrates <- length(labeled_substrates)
  
  log_message(sprintf("Substrate design: %d substrates, %d control types",
                      n_substrates, length(control_isotopologs)), log_file)
  
  # Determine design type
  if (n_substrates == 0) {
    design_info$type <- "no_substrates"
    design_info$control_type <- "unknown"
    
  } else if (n_substrates == 1) {
    design_info$type <- "single_substrate"
    design_info$control_type <- "standard"
    
  } else if (any(control_isotopologs %in% c("several", "mixed", "all", "control"))) {
    design_info$type <- "shared_control"
    design_info$control_type <- "pooled"
    log_message("  Detected SHARED CONTROL design", log_file)
    
  } else {
    # Check for matched controls
    substrate_names <- gsub("^[0-9]+[A-Z]+-", "", labeled_substrates)
    overlap <- length(intersect(substrate_names, control_isotopologs))
    
    if (overlap >= length(substrate_names) * 0.5) {
      design_info$type <- "matched_control"
      design_info$control_type <- "substrate_specific"
      log_message("  Detected MATCHED CONTROL design", log_file)
    } else {
      design_info$type <- "mixed"
      design_info$control_type <- "partial"
    }
  }
  
  design_info$substrates <- labeled_substrates
  design_info$n_substrates <- n_substrates
  design_info$control_isotopologs <- control_isotopologs
  
  # Calculate sample sizes per substrate
  if (n_substrates > 0) {
    substrate_sample_sizes <- list()
    for (sub in labeled_substrates) {
      n_samples <- sum(meta$isotopolog == sub & sample_info$labeled, na.rm = TRUE)
      substrate_sample_sizes[[sub]] <- n_samples
    }
    design_info$substrate_sample_sizes <- substrate_sample_sizes
    
    # Validate substrate sample sizes
    min_samples <- min(unlist(substrate_sample_sizes))
    if (min_samples < BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
      log_message(sprintf("  ⚠ Some substrates have few samples (min: %d)", min_samples),
                  log_file, "WARNING")
    }
  }
  
  return(design_info)
}

# ============================================================================
# QUALITY CONTROL CHECKS
# ============================================================================

#' Perform Qc Checks
#'
#' @description [Function description]
#' @export
perform_qc_checks <- function(physeq, sample_info, study_name, log_file = NULL) {
  log_message("----------------------------------------", log_file)
  log_message("QUALITY CONTROL CHECKS", log_file)
  
  qc_results <- list(
    passed = TRUE,
    warnings = character(),
    errors = character()
  )
  
  # Check 1: Sample sizes
  if (sample_info$n_labeled < BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
    msg <- sprintf("Insufficient labeled samples: %d", sample_info$n_labeled)
    qc_results$errors <- c(qc_results$errors, msg)
    qc_results$passed <- FALSE
    log_message(paste("  ✗", msg), log_file, "ERROR")
  }
  
  if (sample_info$n_unlabeled < BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
    msg <- sprintf("Insufficient unlabeled samples: %d", sample_info$n_unlabeled)
    qc_results$errors <- c(qc_results$errors, msg)
    qc_results$passed <- FALSE
    log_message(paste("  ✗", msg), log_file, "ERROR")
  }
  
  # Check 2: Library sizes
  lib_sizes <- sample_sums(physeq)
  low_lib <- sum(lib_sizes < 1000)
  if (low_lib > 0) {
    msg <- sprintf("%d samples with <1000 reads", low_lib)
    qc_results$warnings <- c(qc_results$warnings, msg)
    log_message(paste("  ⚠", msg), log_file, "WARNING")
  }
  
  # Check 3: Sparsity
  otu_mat <- as(otu_table(physeq), "matrix")
  sparsity <- sum(otu_mat == 0) / length(otu_mat)
  if (sparsity > BIOLOGICAL_PARAMS$filtering$max_sparsity) {
    msg <- sprintf("Very sparse data: %.1f%% zeros", sparsity * 100)
    qc_results$warnings <- c(qc_results$warnings, msg)
    log_message(paste("  ⚠", msg), log_file, "WARNING")
  }
  
  # Check 4: Taxa count
  n_taxa <- ntaxa(physeq)
  if (n_taxa < 100) {
    msg <- sprintf("Low taxa count: %d", n_taxa)
    qc_results$warnings <- c(qc_results$warnings, msg)
    log_message(paste("  ⚠", msg), log_file, "WARNING")
  }
  
  # Check 5: Detection confidence
  if (sample_info$detection_confidence == "low") {
    msg <- "Low confidence in sample classification"
    qc_results$warnings <- c(qc_results$warnings, msg)
    log_message(paste("  ⚠", msg), log_file, "WARNING")
  }
  
  # Summary
  if (qc_results$passed && length(qc_results$warnings) == 0) {
    log_message("  ✓ All QC checks passed", log_file, "SUCCESS")
  } else if (qc_results$passed) {
    log_message(sprintf("  ✓ QC passed with %d warnings", length(qc_results$warnings)),
                log_file, "WARNING")
  } else {
    log_message(sprintf("  ✗ QC failed: %d errors, %d warnings",
                        length(qc_results$errors), length(qc_results$warnings)),
                log_file, "ERROR")
  }
  
  return(qc_results)
}

# ============================================================================
# METADATA EXTRACTION
# ============================================================================

#' Extract Comprehensive Metadata
#'
#' @description [Function description]
#' @export
extract_comprehensive_metadata <- function(physeq, sample_info, study_name,
                                           substrate_specific = FALSE,
                                           current_substrate = NULL) {
  meta <- sample_info$metadata
  
  study_metadata <- list(
    study_id = study_name,
    pipeline_version = PIPELINE_VERSION,
    analysis_date = as.character(Sys.Date()),  # Keep as character
    n_labeled_total = sample_info$n_labeled,
    n_unlabeled_total = sample_info$n_unlabeled,
    total_samples = nrow(meta),
    detection_confidence = sample_info$detection_confidence,
    substrate_specific = substrate_specific
  )
  
  # Isotope information
  if (length(sample_info$isotopes) > 0) {
    study_metadata$isotopes_list <- paste(sample_info$isotopes, collapse = ";")
    study_metadata$primary_isotope <- sample_info$isotopes[1]
  } else {
    study_metadata$isotopes_list <- NA
    study_metadata$primary_isotope <- NA
  }
  
  # Substrate information
  if (substrate_specific && !is.null(current_substrate)) {
    study_metadata$isotopolog <- current_substrate
    study_metadata$substrates_list <- current_substrate
    study_metadata$current_substrate <- current_substrate
  } else if (length(sample_info$substrates) > 0) {
    study_metadata$substrates_list <- paste(sample_info$substrates, collapse = ";")
    study_metadata$isotopolog <- sample_info$substrates[1]
    study_metadata$n_substrates <- length(sample_info$substrates)
  } else {
    study_metadata$isotopolog <- NA
    study_metadata$substrates_list <- NA
    study_metadata$n_substrates <- 0
  }
  
  # Environmental metadata (WITHOUT DATE COLUMNS)
  env_mapping <- list(
    env_biome = c("Env.Biome..sample.", "env_biome", "Environment", "Env.Biome"),
    env_feature = c("Env.Feature..sample.", "env_feature", "Env.Feature"),
    env_material = c("Env.Material..sample.", "env_material", "Env.Material"),
    environment_label = c("Environment", "environment_label", "environment", "env_label")
  )
  
  for (meta_col in names(env_mapping)) {
    possible_cols <- env_mapping[[meta_col]]
    found_col <- intersect(possible_cols, colnames(meta))
    
    if (length(found_col) > 0) {
      # Ensure we convert to character to avoid date parsing
      values <- as.character(meta[[found_col[1]]])
      values <- unique(values[!is.na(values) & values != "" & values != "NA"])
      study_metadata[[meta_col]] <- if (length(values) > 0) {
        paste(values[1:min(3, length(values))], collapse = "; ")
      } else NA
    } else {
      study_metadata[[meta_col]] <- NA
    }
  }
  
  # Publication metadata
  doi_cols <- c("DOI.URL.x", "DOI.URL.y", "DOI", "DOI.URL", "DOI_URL", "doi")
  for (col in doi_cols) {
    if (col %in% colnames(meta)) {
      dois <- as.character(unique(meta[[col]]))
      dois <- dois[!is.na(dois) & dois != "" & dois != "NA"]
      if (length(dois) > 0) {
        study_metadata$DOI_URL <- dois[1]
        break
      }
    }
  }
  
  # Bioproject
  bioproject_cols <- c("Bioproject.ID", "Bioproject", "Bioproject_ID", "bioproject_id")
  for (col in bioproject_cols) {
    if (col %in% colnames(meta)) {
      bioprojects <- as.character(unique(meta[[col]]))
      bioprojects <- bioprojects[!is.na(bioprojects) & bioprojects != ""]
      if (length(bioprojects) > 0) {
        study_metadata$Bioproject_ID <- bioprojects[1]
        break
      }
    }
  }
  
  # Gradient information
  if ("gradient_position" %in% colnames(meta)) {
    positions <- unique(meta$gradient_position[meta$gradient_position > 0 &
                                                 !is.na(meta$gradient_position)])
    study_metadata$n_fractions <- length(positions)
  } else {
    study_metadata$n_fractions <- NA
  }
  
  if ("gradient_pos_density" %in% colnames(meta)) {
    densities <- parse_density_values(meta$gradient_pos_density)
    study_metadata$has_density <- any(!is.na(densities) & densities > 0)
    if (study_metadata$has_density) {
      study_metadata$density_range <- paste(
        round(range(densities[!is.na(densities)]), 3),
        collapse = "-"
      )
    }
  } else {
    study_metadata$has_density <- FALSE
  }
  
  return(study_metadata)
}

# ============================================================================
# CACHE MANAGEMENT
# ============================================================================

#' Get Cached Taxonomy
#'
#' @description [Function description]
#' @export
get_cached_taxonomy <- function(physeq, study_id) {
  cache_key <- paste0(study_id, "_tax")
  
  if (exists(cache_key, envir = .GlobalEnv$.sip_taxonomy_cache)) {
    return(get(cache_key, envir = .GlobalEnv$.sip_taxonomy_cache))
  }
  
  tax_info <- tryCatch({
    tax_df <- as.data.frame(tax_table(physeq))
    tax_df$taxa <- rownames(tax_df)
    setDT(tax_df)
    setkey(tax_df, taxa)
    tax_df
  }, error = function(e) {
    data.table(taxa = character(0))
  })
  
  assign(cache_key, tax_info, envir = .GlobalEnv$.sip_taxonomy_cache)
  return(tax_info)
}

#' Clear Cache
#'
#' @description [Function description]
#' @export
clear_cache <- function() {
  rm(list = ls(envir = .GlobalEnv$.sip_taxonomy_cache),
     envir = .GlobalEnv$.sip_taxonomy_cache)
  rm(list = ls(envir = .GlobalEnv$.sip_validation_cache),
     envir = .GlobalEnv$.sip_validation_cache)
  message("✓ All caches cleared")
  invisible(NULL)
}

################################################################################
# SIP Meta-Analysis Pipeline v1.0
################################################################################

#' Apply Smart Filtering
#'
#' @description [Function description]
#' @export
apply_smart_filtering <- function(otu_mat, study_type, log_file = NULL) {
  n_samples <- ncol(otu_mat)
  n_taxa_original <- nrow(otu_mat)
  
  keep_taxa <- rowSums(otu_mat) > 0
  
  otu_mat_filtered <- otu_mat[keep_taxa, ]
  n_taxa_after <- nrow(otu_mat_filtered)
  
  log_message(sprintf("  Taxa: %d → %d (%.1f%% retained)",
                      n_taxa_original, n_taxa_after,
                      100 * n_taxa_after/n_taxa_original),
              log_file)
  
  # Add small pseudocount for numerical stability
  otu_mat_filtered <- otu_mat_filtered + BIOLOGICAL_PARAMS$filtering$pseudocount
  
  return(otu_mat_filtered)
}
# ============================================================================
# DESEQ2 WITH ROBUST ERROR HANDLING
# ============================================================================

#' Run Deseq2 Robust
#'
#' @description [Function description]
#' @export
run_deseq2_robust <- function(otu_mat, condition, study_name, log_file = NULL) {
  log_message("  Running DESeq2...", log_file)
  
  tryCatch({
    counts <- round(otu_mat)
    counts[counts < 0] <- 0
    
    colData <- data.frame(
      row.names = colnames(counts),
      condition = condition
    )
    
    # Create DESeq dataset
    dds <- DESeqDataSetFromMatrix(
      countData = counts,
      colData = colData,
      design = ~ condition
    )
    
    # Robust size factor estimation
    dds <- tryCatch({
      # Estimate normalization factors from count data
      estimateSizeFactors(dds)
    }, error = function(e) {
      log_message("    Standard size factors failed, trying poscounts", log_file, "WARNING")
      tryCatch({
        # Estimate normalization factors from count data
        estimateSizeFactors(dds, type = "poscounts")
      }, error = function(e2) {
        log_message("    Poscounts failed, using median ratio", log_file, "WARNING")
        lib_sizes <- colSums(counts)
        size_factors <- lib_sizes / median(lib_sizes)
        size_factors[size_factors == 0] <- 1
        sizeFactors(dds) <- size_factors
        dds
      })
    })
    
    # Run DESeq2 with appropriate settings
    dds <- tryCatch({
      DESeq(dds, quiet = TRUE, parallel = FALSE)
    }, error = function(e) {
      if (grepl("dispersion", e$message)) {
        log_message("    Using gene-wise dispersion", log_file, "WARNING")
        # Estimate dispersion parameters for negative binomial model
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        # Perform Wald test for differential expression
        dds <- nbinomWaldTest(dds)
      } else if (grepl("convergence", e$message)) {
        log_message("    Using local fit", log_file, "WARNING")
        DESeq(dds, fitType = "local", quiet = TRUE)
      } else {
        stop(e)
      }
    })
    
    # Extract results
    res <- results(dds, contrast = c("condition", "Labeled", "Unlabeled"))
    
    result_df <- data.frame(
      method = "DESeq2",
      taxa = rownames(res),
      log2FC = res$log2FoldChange,
      pvalue = res$pvalue,
      padj = res$padj,
      baseMean = res$baseMean,
      lfcSE = res$lfcSE,
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(padj), !is.na(log2FC), is.finite(log2FC))
    
    log_message(sprintf("    ✓ DESeq2 completed: %d results", nrow(result_df)),
                log_file, "SUCCESS")
    
    return(result_df)
    
  }, error = function(e) {
    log_message(sprintf("    ✗ DESeq2 failed: %s", e$message), log_file, "ERROR")
    return(NULL)
  })
}

# ============================================================================
# EDGER WITH ROBUST SETTINGS
# ============================================================================

#' Run Edger Robust
#'
#' @description [Function description]
#' @export
run_edger_robust <- function(otu_mat, condition, study_name, log_file = NULL) {
  log_message("  Running edgeR...", log_file)
  
  tryCatch({
    y <- DGEList(counts = otu_mat, group = condition)
    
    # Adaptive filtering
    keep <- rowSums(cpm(y) > 0.5) >= 2
    if (sum(keep) < 10) {
      keep <- rowSums(y$counts > 0) >= 2
    }
    
    y <- y[keep, , keep.lib.sizes = FALSE]
    
    # Normalization
    # Calculate TMM normalization factors
    y <- calcNormFactors(y, method = "TMM")
    
    # Design matrix
    design <- model.matrix(~ condition)
    
    # Robust dispersion estimation
    y <- estimateDisp(y, design, robust = TRUE)
    
    # Fit model
    # Fit quasi-likelihood negative binomial model
    fit <- glmQLFit(y, design, robust = TRUE)
    # Perform quasi-likelihood F-test
    qlf <- glmQLFTest(fit, coef = 2)
    
    # Extract results
    # Extract top differential abundant taxa
    res <- topTags(qlf, n = Inf, adjust.method = "BH")$table
    
    result_df <- data.frame(
      method = "edgeR",
      taxa = rownames(res),
      log2FC = res$logFC,
      pvalue = res$PValue,
      padj = res$FDR,
      logCPM = res$logCPM,
      F = res$F,
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(padj), !is.na(log2FC), is.finite(log2FC))
    
    log_message(sprintf("    ✓ edgeR completed: %d results", nrow(result_df)),
                log_file, "SUCCESS")
    
    return(result_df)
    
  }, error = function(e) {
    log_message(sprintf("    ✗ edgeR failed: %s", e$message), log_file, "ERROR")
    return(NULL)
  })
}

# ============================================================================
# LIMMA-VOOM WITH QUALITY WEIGHTS
# ============================================================================

#' Run Limma Robust
#'
#' @description [Function description]
#' @export
run_limma_robust <- function(otu_mat, condition, study_name, log_file = NULL) {
  log_message("  Running limma-voom...", log_file)
  
  tryCatch({
    dge <- DGEList(counts = otu_mat)
    
    # Filter
    keep <- filterByExpr(dge, group = condition)
    if (sum(keep) < 10) {
      keep <- rowSums(dge$counts > 0) >= 2
    }
    
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    
    # Normalize
    # Calculate TMM normalization factors
    dge <- calcNormFactors(dge, method = "TMM")
    
    # Design
    design <- model.matrix(~ condition)
    
    # Voom with quality weights
    v <- tryCatch({
      voomWithQualityWeights(dge, design, plot = FALSE)
    }, error = function(e) {
      log_message("    Standard voom failed, trying without weights", log_file, "WARNING")
      # Transform counts to log2-CPM with precision weights
      voom(dge, design, plot = FALSE)
    })
    
    # Fit model
    fit <- lmFit(v, design)
    fit <- eBayes(fit, robust = TRUE)
    
    # Extract results
    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    
    result_df <- data.frame(
      method = "limma",
      taxa = rownames(res),
      log2FC = res$logFC,
      pvalue = res$P.Value,
      padj = res$adj.P.Val,
      AveExpr = res$AveExpr,
      t = res$t,
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(padj), !is.na(log2FC), is.finite(log2FC))
    
    log_message(sprintf("    ✓ limma-voom completed: %d results", nrow(result_df)),
                log_file, "SUCCESS")
    
    return(result_df)
    
  }, error = function(e) {
    log_message(sprintf("    ✗ limma-voom failed: %s", e$message), log_file, "ERROR")
    return(NULL)
  })
}

# ============================================================================
# ALDEX2 WITH SMART SELECTION
# ============================================================================

#' Should Run Aldex2
#'
#' @description [Function description]
#' @export
should_run_aldex2 <- function(n_taxa, n_samples, study_type) {
  # Decision based on computational complexity
  complexity <- n_taxa * n_samples
  
  if (complexity > 50000) return(FALSE)
  if (n_taxa > 3000) return(FALSE)
  if (n_samples > 50) return(FALSE)
  if (n_samples < 4) return(FALSE)
  if (study_type$type == "binary" && n_samples < 6) return(FALSE)
  
  return(TRUE)
}

#' Run Aldex2 Smart
#'
#' @description [Function description]
#' @export
run_aldex2_smart <- function(otu_mat, labeled_subset, unlabeled_subset,
                             study_name, study_type, log_file = NULL) {
  
  n_taxa <- nrow(otu_mat)
  n_samples <- ncol(otu_mat)
  
  if (!should_run_aldex2(n_taxa, n_samples, study_type)) {
    log_message(sprintf("  Skipping ALDEx2 (complexity: %d × %d = %d)",
                        n_taxa, n_samples, n_taxa * n_samples),
                log_file, "WARNING")
    return(NULL)
  }
  
  log_message("  Running ALDEx2...", log_file)
  
  tryCatch({
    # Prepare data
    otu_mat_int <- round(otu_mat)
    
    # Create conditions vector
    conds <- character(n_samples)
    conds[labeled_subset] <- "Labeled"
    conds[unlabeled_subset] <- "Unlabeled"
    
    # Adaptive MC samples
    mc_samples <- min(16, max(4, floor(128 / sqrt(n_taxa))))
    
    log_message(sprintf("    Using %d MC samples", mc_samples), log_file)
    
    # Run ALDEx2
    # Run ALDEx2 compositional analysis
    x <- aldex(
      reads = otu_mat_int,
      conditions = conds,
      mc.samples = mc_samples,
      test = "t",
      effect = TRUE,
      denom = "iqlr",
      verbose = FALSE
    )
    
    result_df <- data.frame(
      method = "ALDEx2",
      taxa = rownames(x),
      log2FC = x$diff.btw,
      pvalue = x$we.ep,
      padj = x$we.eBH,
      effect = x$effect,
      overlap = x$overlap,
      row.names = NULL,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(padj), !is.na(log2FC), is.finite(log2FC))
    
    log_message(sprintf("    ✓ ALDEx2 completed: %d results", nrow(result_df)),
                log_file, "SUCCESS")
    
    return(result_df)
    
  }, error = function(e) {
    log_message(sprintf("    ✗ ALDEx2 failed: %s", e$message), log_file, "ERROR")
    return(NULL)
  })
}

# ============================================================================
# MAIN DA ANALYSIS COORDINATOR
# ============================================================================

#' Run Da Analysis
#'
#' @description [Function description]
#' @export
run_da_analysis <- function(physeq, labeled_samples, unlabeled_samples,
                            study_name, study_type = NULL) {
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_name, "_da_analysis.log"))
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message(sprintf("DA ANALYSIS: %s", study_name), log_file)
  log_message("========================================", log_file)
  
  # Validate inputs
  if (!validate_sample_sizes(labeled_samples, unlabeled_samples,
                             BIOLOGICAL_PARAMS$statistics$min_samples_per_group,
                             log_file)) {
    return(NULL)
  }
  
  # Prepare data
  otu_mat <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq)) otu_mat <- t(otu_mat)
  
  # Handle NAs
  if (any(is.na(otu_mat))) {
    otu_mat[is.na(otu_mat)] <- 0
    log_message("  Replaced NA values with 0", log_file, "WARNING")
  }
  
  # Subset to relevant samples
  all_samples <- labeled_samples | unlabeled_samples
  otu_mat <- otu_mat[, all_samples]
  labeled_subset <- labeled_samples[all_samples]
  unlabeled_subset <- unlabeled_samples[all_samples]
  
  # Apply smart filtering
  if (is.null(study_type)) {
    study_type <- list(type = "unknown", quality = "medium")
  }
  
  otu_mat <- apply_smart_filtering(otu_mat, study_type, log_file)
  
  if (nrow(otu_mat) < 2) {
    log_message("Insufficient taxa after filtering", log_file, "ERROR")
    return(NULL)
  }
  
  # Create condition factor
  condition <- factor(
    ifelse(labeled_subset, "Labeled", "Unlabeled"),
    levels = c("Unlabeled", "Labeled")
  )
  
  # Initialize results
  results_list <- list()
  
  # Run each method
  log_message("Running differential abundance methods:", log_file)
  
  # DESeq2
  results_list$deseq2 <- run_deseq2_robust(otu_mat, condition, study_name, log_file)
  
  # edgeR
  results_list$edger <- run_edger_robust(otu_mat, condition, study_name, log_file)
  
  # limma-voom
  results_list$limma <- run_limma_robust(otu_mat, condition, study_name, log_file)
  
  # ALDEx2 (conditional)
  results_list$aldex2 <- run_aldex2_smart(otu_mat, labeled_subset, unlabeled_subset,
                                          study_name, study_type, log_file)
  
  # Combine results
  all_results <- lapply(results_list[!sapply(results_list, is.null)], standardize_results)
  results <- bind_rows(all_results)
  
  if (nrow(results) == 0) {
    log_message("No results generated from any method", log_file, "ERROR")
    return(NULL)
  }
  
  # Add metadata
  results$n_labeled <- sum(labeled_subset)
  results$n_unlabeled <- sum(unlabeled_subset)
  results$study_name <- study_name
  
  # Calculate consensus metrics
  consensus <- results %>%
    group_by(taxa) %>%
    summarise(
      n_methods = n(),
      mean_log2FC = mean(log2FC, na.rm = TRUE),
      median_log2FC = median(log2FC, na.rm = TRUE),
      min_padj = min(padj, na.rm = TRUE),
      max_padj = max(padj, na.rm = TRUE),
      consensus_direction = ifelse(all(log2FC > 0), "up",
                                   ifelse(all(log2FC < 0), "down", "mixed")),
      .groups = "drop"
    )
  
  log_message(sprintf("Summary: %d methods, %d unique taxa, %d consensus enriched (padj<0.05)",
                      length(all_results),
                      n_distinct(results$taxa),
                      sum(consensus$min_padj < 0.05, na.rm = TRUE)),
              log_file, "SUCCESS")
  
  # Clean up
  rm(otu_mat)
  gc(verbose = FALSE)
  
  return(results)
}

################################################################################
# SIP Meta-Analysis Pipeline v1.0
################################################################################

# ============================================================================
# WINDOW GENERATION STRATEGIES
# ============================================================================

#' Generate biologically-informed windows for density gradients
generate_density_windows <- function(densities, labeled, unlabeled,
                                     heavy_threshold = NULL,
                                     log_file = NULL) {
  
  if (is.null(heavy_threshold)) {
    heavy_threshold <- BIOLOGICAL_PARAMS$density$heavy_threshold
  }
  
  valid_mask <- !is.na(densities) & densities > 0
  valid_densities <- densities[valid_mask]
  
  if (length(valid_densities) < 4) {
    log_message("Insufficient density values for windowing", log_file, "WARNING")
    return(NULL)
  }
  
  windows <- list()
  
  # ONLY analyze heavy fraction - this is what matters for SIP
  heavy_mask <- densities >= heavy_threshold & valid_mask
  if (sum(heavy_mask & labeled) >= 2 && sum(heavy_mask & unlabeled) >= 2) {
    windows[["heavy_fraction"]] <- list(
      indices = which(heavy_mask),
      center = mean(densities[heavy_mask], na.rm = TRUE),
      range = range(densities[heavy_mask], na.rm = TRUE),
      type = "heavy",
      priority = 1
    )
  }
  
  # Optional: Add a slightly lower threshold window for borderline incorporators
  medium_threshold <- heavy_threshold - 0.01
  borderline_mask <- densities >= medium_threshold & densities < heavy_threshold & valid_mask
  if (sum(borderline_mask & labeled) >= 2 && sum(borderline_mask & unlabeled) >= 2) {
    windows[["borderline"]] <- list(
      indices = which(borderline_mask),
      center = mean(densities[borderline_mask], na.rm = TRUE),
      range = range(densities[borderline_mask], na.rm = TRUE),
      type = "borderline",
      priority = 2
    )
  }
  
  log_message(sprintf("Generated %d density windows for SIP analysis",
                      length(windows)), log_file)
  
  return(windows)
}

#' Generate windows for position-based gradients
generate_position_windows <- function(positions, labeled, unlabeled,
                                      n_positions, log_file = NULL) {
  
  valid_mask <- !is.na(positions) & positions > 0
  unique_positions <- sort(unique(positions[valid_mask]))
  
  if (length(unique_positions) < 2) {
    return(NULL)
  }
  
  windows <- list()
  
  if (n_positions <= 3) {
    # For few positions, use each as a window
    for (i in seq_along(unique_positions)) {
      pos <- unique_positions[i]
      pos_mask <- positions == pos & valid_mask
      
      if (sum(pos_mask & labeled) >= 2 && sum(pos_mask & unlabeled) >= 2) {
        windows[[paste0("position_", pos)]] <- list(
          indices = which(pos_mask),
          center = pos,
          range = c(pos, pos),
          type = "single_position",
          priority = i
        )
      }
    }
    
  } else {
    # For many positions, group into windows
    # Use hierarchical clustering to find natural groups
    pos_dist <- dist(unique_positions)
    
    if (length(unique_positions) >= 5) {
      # Adaptive number of clusters
      n_clusters <- min(ceiling(length(unique_positions) / 3), 5)
      hc <- hclust(pos_dist, method = "ward.D2")
      clusters <- cutree(hc, k = n_clusters)
      
      for (k in 1:n_clusters) {
        cluster_positions <- unique_positions[clusters == k]
        pos_mask <- positions %in% cluster_positions & valid_mask
        
        if (sum(pos_mask & labeled) >= 2 && sum(pos_mask & unlabeled) >= 2) {
          windows[[paste0("cluster_", k)]] <- list(
            indices = which(pos_mask),
            center = mean(cluster_positions),
            range = range(cluster_positions),
            type = "position_cluster",
            priority = k
          )
        }
      }
    } else {
      # Simple adjacent grouping for 4-5 positions
      window_size <- 2
      for (i in 1:(length(unique_positions) - window_size + 1)) {
        window_positions <- unique_positions[i:(i + window_size - 1)]
        pos_mask <- positions %in% window_positions & valid_mask
        
        if (sum(pos_mask & labeled) >= 2 && sum(pos_mask & unlabeled) >= 2) {
          windows[[paste0("window_", i)]] <- list(
            indices = which(pos_mask),
            center = mean(window_positions),
            range = range(window_positions),
            type = "adjacent_positions",
            priority = i
          )
        }
      }
    }
  }
  
  log_message(sprintf("Generated %d position windows", length(windows)), log_file)
  
  return(windows)
}

#' Adaptive KNN window generation (only for specific cases)
generate_knn_windows <- function(values, labeled, unlabeled,
                                 k_per_group = 3,
                                 use_all_centers = FALSE,
                                 log_file = NULL) {
  
  valid_mask <- !is.na(values) & is.finite(values)
  valid_values <- values[valid_mask]
  valid_labeled <- labeled[valid_mask]
  valid_unlabeled <- unlabeled[valid_mask]
  
  n_labeled <- sum(valid_labeled)
  n_unlabeled <- sum(valid_unlabeled)
  
  # Adaptive k based on sample sizes
  k_per_group <- min(k_per_group, floor(min(n_labeled, n_unlabeled) / 3))
  
  if (k_per_group < 2) {
    log_message("Insufficient samples for KNN windows", log_file, "WARNING")
    return(NULL)
  }
  
  # Order by value
  ord <- order(valid_values)
  ordered_values <- valid_values[ord]
  ordered_labeled <- valid_labeled[ord]
  ordered_unlabeled <- valid_unlabeled[ord]
  
  windows <- list()
  
  # Select center points
  if (use_all_centers) {
    centers <- seq_along(ordered_values)
  } else {
    # Use strategic centers (quartiles + extremes)
    n_centers <- min(10, length(ordered_values))
    centers <- unique(round(seq(1, length(ordered_values), length.out = n_centers)))
  }
  
  for (i in seq_along(centers)) {
    center_idx <- centers[i]
    center_val <- ordered_values[center_idx]
    
    # Calculate distances
    distances <- abs(ordered_values - center_val)
    
    # Find k nearest labeled and unlabeled
    labeled_idx <- which(ordered_labeled)
    unlabeled_idx <- which(ordered_unlabeled)
    
    if (length(labeled_idx) < k_per_group || length(unlabeled_idx) < k_per_group) {
      next
    }
    
    # Get k nearest from each group
    labeled_dists <- distances[labeled_idx]
    unlabeled_dists <- distances[unlabeled_idx]
    
    k_labeled <- labeled_idx[order(labeled_dists)[1:k_per_group]]
    k_unlabeled <- unlabeled_idx[order(unlabeled_dists)[1:k_per_group]]
    
    window_idx <- sort(unique(c(k_labeled, k_unlabeled)))
    
    # Map back to original indices
    original_idx <- ord[window_idx]
    
    windows[[paste0("knn_", i)]] <- list(
      indices = original_idx,
      center = center_val,
      range = range(ordered_values[window_idx]),
      type = "knn",
      priority = abs(i - length(centers)/2) + 1,  # Prioritize middle windows
      k = k_per_group
    )
  }
  
  log_message(sprintf("Generated %d KNN windows (k=%d per group)",
                      length(windows), k_per_group), log_file)
  
  return(windows)
}

# ============================================================================
# ADAPTIVE WINDOW SELECTION
# ============================================================================

#' Select Window Strategy
#'
#' @description [Function description]
#' @export
select_window_strategy <- function(physeq, sample_info, study_type, variable,
                                   log_file = NULL) {
  
  log_message("========================================", log_file)
  log_message("WINDOW STRATEGY SELECTION", log_file)
  log_message("========================================", log_file)
  
  meta <- sample_info$metadata
  strategy <- list()
  
  if (!(variable %in% colnames(meta))) {
    log_message(sprintf("Variable '%s' not found", variable), log_file, "ERROR")
    return(NULL)
  }
  
  # Parse values
  if (variable == "gradient_pos_density") {
    values <- parse_density_values(meta[[variable]])
  } else {
    values <- as.numeric(meta[[variable]])
  }
  
  # Get SIP samples only
  sip_mask <- sample_info$labeled | sample_info$unlabeled
  values_sip <- values[sip_mask]
  labeled_sip <- sample_info$labeled[sip_mask]
  unlabeled_sip <- sample_info$unlabeled[sip_mask]
  
  # Count valid values
  valid_mask <- !is.na(values_sip) & is.finite(values_sip) & values_sip > 0
  n_valid <- sum(valid_mask)
  
  log_message(sprintf("Valid gradient values: %d/%d", n_valid, length(values_sip)), log_file)
  
  if (n_valid < 4) {
    log_message("Insufficient valid values for windowing", log_file, "ERROR")
    return(NULL)
  }
  
  # Determine strategy based on study type and data quality
  if (study_type$type == "density") {
    if (study_type$density_quality == "high") {
      strategy$method <- "fixed_density"
      strategy$description <- "Fixed biological windows for high-quality density data"
    } else if (study_type$density_quality == "medium") {
      strategy$method <- "adaptive_density"
      strategy$description <- "Adaptive density windows for medium-quality data"
    } else {
      strategy$method <- "knn_density"
      strategy$description <- "KNN windows for low-quality density data"
    }
    
  } else if (study_type$type == "binary") {
    strategy$method <- "binary_direct"
    strategy$description <- "Direct comparison for binary design"
    
  } else if (study_type$type == "multifraction") {
    if (study_type$n_positions <= 5) {
      strategy$method <- "position_grouped"
      strategy$description <- "Grouped positions for few fractions"
    } else {
      strategy$method <- "position_clustered"
      strategy$description <- "Hierarchical clustering for many fractions"
    }
  } else {
    strategy$method <- "fallback_knn"
    strategy$description <- "Fallback KNN strategy"
  }
  
  strategy$variable <- variable
  strategy$n_valid_samples <- n_valid
  
  log_message(sprintf("Selected strategy: %s", strategy$description), log_file, "SUCCESS")
  
  return(strategy)
}

# ============================================================================
# RUN WINDOWED ANALYSIS
# ============================================================================

#' Run Windowed Da Analysis
#'
#' @description [Function description]
#' @export
run_windowed_da_analysis <- function(physeq, sample_info, study_type,
                                     variable, strategy,
                                     study_name, substrate_specific = FALSE,
                                     current_substrate = NULL,
                                     n_cores = NULL,
                                     log_file = NULL) {
  
  log_message("========================================", log_file)
  log_message("WINDOWED ANALYSIS", log_file)
  log_message("========================================", log_file)
  
  meta <- sample_info$metadata
  
  # Parse gradient values
  if (variable == "gradient_pos_density") {
    values <- parse_density_values(meta[[variable]])
  } else {
    values <- as.numeric(meta[[variable]])
  }
  
  # Filter for relevant samples
  if (substrate_specific && !is.null(current_substrate)) {
    substrate_mask <- (meta$isotopolog == current_substrate & sample_info$labeled) |
      sample_info$unlabeled
    relevant_mask <- substrate_mask & !is.na(values) & values > 0
  } else {
    relevant_mask <- (sample_info$labeled | sample_info$unlabeled) &
      !is.na(values) & values > 0
  }
  
  values_relevant <- values[relevant_mask]
  labeled_relevant <- sample_info$labeled[relevant_mask]
  unlabeled_relevant <- sample_info$unlabeled[relevant_mask]
  sample_ids_relevant <- rownames(meta)[relevant_mask]
  
  # Generate windows based on strategy
  windows <- NULL
  
  if (strategy$method %in% c("fixed_density", "adaptive_density")) {
    windows <- generate_density_windows(
      values_relevant,
      labeled_relevant,
      unlabeled_relevant,
      heavy_threshold = BIOLOGICAL_PARAMS$density$heavy_threshold,
      log_file = log_file
    )
    
  } else if (strategy$method == "binary_direct") {
    # For binary, just use all samples as one window
    windows <- list(
      all = list(
        indices = seq_along(values_relevant),
        center = mean(values_relevant, na.rm = TRUE),
        range = range(values_relevant, na.rm = TRUE),
        type = "binary",
        priority = 1
      )
    )
    
  } else if (strategy$method %in% c("position_grouped", "position_clustered")) {
    windows <- generate_position_windows(
      values_relevant,
      labeled_relevant,
      unlabeled_relevant,
      study_type$n_positions,
      log_file = log_file
    )
    
  } else {
    # Fallback to KNN
    windows <- generate_knn_windows(
      values_relevant,
      labeled_relevant,
      unlabeled_relevant,
      k_per_group = max(3, floor(sum(labeled_relevant) / 5)),
      use_all_centers = FALSE,
      log_file = log_file
    )
  }
  
  if (is.null(windows) || length(windows) == 0) {
    log_message("No valid windows generated", log_file, "ERROR")
    return(NULL)
  }
  
  log_message(sprintf("Processing %d windows", length(windows)), log_file)
  
  # Run DA analysis for each window
  window_results <- list()
  window_metadata <- list()
  
  # Determine if parallel processing is beneficial
  use_parallel <- FALSE
  if (!is.null(n_cores) && n_cores > 1 && length(windows) >= 3) {
    use_parallel <- TRUE
    n_cores <- min(n_cores, length(windows))
  }
  
  if (use_parallel) {
    log_message(sprintf("Using parallel processing (%d cores)", n_cores), log_file)
    
    # Setup parallel backend
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    
    # Export required objects
    clusterExport(cl, c("physeq", "sample_ids_relevant", "run_da_analysis",
                        "validate_sample_sizes", "apply_smart_filtering",
                        "run_deseq2_robust", "run_edger_robust",
                        "run_limma_robust", "run_aldex2_smart", "should_run_aldex2",
                        "standardize_results", "BIOLOGICAL_PARAMS",
                        "main_dir", "log_message"),
                  envir = environment())
    
    # Run in parallel
    results <- foreach(window_name = names(windows),
                       .combine = 'c',
                       .packages = c("phyloseq", "tidyverse", "DESeq2",
                                     "edgeR", "limma", "ALDEx2")) %dopar% {
                                       
                                       window <- windows[[window_name]]
                                       window_samples <- sample_ids_relevant[window$indices]
                                       
                                       # Subset phyloseq
                                       physeq_window <- prune_samples(window_samples, physeq)
                                       physeq_window <- prune_taxa(taxa_sums(physeq_window) > 0, physeq_window)
                                       
                                       # Determine labeled/unlabeled in window
                                       window_meta <- as.data.frame(sample_data(physeq_window))
                                       labeled_window <- labeled_relevant[window$indices][match(rownames(window_meta),
                                                                                                window_samples)]
                                       unlabeled_window <- unlabeled_relevant[window$indices][match(rownames(window_meta),
                                                                                                    window_samples)]
                                       
                                       # Run DA analysis
                                       da_results <- run_da_analysis(
                                         physeq_window,
                                         labeled_window,
                                         unlabeled_window,
                                         paste0(study_name, "_", window_name),
                                         study_type
                                       )
                                       
                                       if (!is.null(da_results) && nrow(da_results) > 0) {
                                         da_results$window_name <- window_name
                                         da_results$window_center <- window$center
                                         da_results$window_range_start <- window$range[1]
                                         da_results$window_range_end <- window$range[2]
                                         da_results$window_type <- window$type
                                         da_results$window_priority <- window$priority
                                         da_results$n_samples_window <- length(window_samples)
                                       }
                                       
                                       list(da_results)
                                     }
    
    stopCluster(cl)
    
    # Collect results
    for (res in results) {
      if (!is.null(res) && nrow(res) > 0) {
        window_name <- unique(res$window_name)[1]
        window_results[[window_name]] <- res
        window_metadata[[window_name]] <- windows[[window_name]]
      }
    }
    
  } else {
    # Sequential processing
    pb <- txtProgressBar(min = 0, max = length(windows), style = 3)
    
    for (i in seq_along(windows)) {
      setTxtProgressBar(pb, i)
      
      window_name <- names(windows)[i]
      window <- windows[[window_name]]
      window_samples <- sample_ids_relevant[window$indices]
      
      # Subset phyloseq
      physeq_window <- prune_samples(window_samples, physeq)
      physeq_window <- prune_taxa(taxa_sums(physeq_window) > 0, physeq_window)
      
      if (ntaxa(physeq_window) < 10) {
        next
      }
      
      # Determine labeled/unlabeled in window
      window_meta <- as.data.frame(sample_data(physeq_window))
      labeled_window <- rownames(window_meta) %in%
        sample_ids_relevant[window$indices][labeled_relevant[window$indices]]
      unlabeled_window <- rownames(window_meta) %in%
        sample_ids_relevant[window$indices][unlabeled_relevant[window$indices]]
      
      # Run DA analysis
      da_results <- run_da_analysis(
        physeq_window,
        labeled_window,
        unlabeled_window,
        paste0(study_name, "_", window_name),
        study_type
      )
      
      if (!is.null(da_results) && nrow(da_results) > 0) {
        da_results$window_name <- window_name
        da_results$window_center <- window$center
        da_results$window_range_start <- window$range[1]
        da_results$window_range_end <- window$range[2]
        da_results$window_type <- window$type
        da_results$window_priority <- window$priority
        da_results$n_samples_window <- length(window_samples)
        
        window_results[[window_name]] <- da_results
        window_metadata[[window_name]] <- window
      }
    }
    
    close(pb)
  }
  
  log_message(sprintf("Completed: %d/%d windows produced results",
                      length(window_results), length(windows)),
              log_file, "SUCCESS")
  
  # Combine results
  if (length(window_results) > 0) {
    combined_results <- bind_rows(window_results)
    
    # Add additional metadata
    combined_results$variable_type <- variable
    combined_results$strategy_method <- strategy$method
    
    if (substrate_specific && !is.null(current_substrate)) {
      combined_results$substrate_analyzed <- current_substrate
    }
    
    return(list(
      combined_results = combined_results,
      window_metadata = window_metadata,
      strategy = strategy,
      n_windows = length(window_results)
    ))
  } else {
    return(NULL)
  }
}

# ============================================================================
# META-ANALYSIS OF WINDOW RESULTS
# ============================================================================

#' Perform Window Meta Analysis
#'
#' @description [Function description]
#' @export
perform_window_meta_analysis <- function(window_results, min_windows = 2,
                                         log_file = NULL) {
  
  if (is.null(window_results) || nrow(window_results) == 0) {
    return(NULL)
  }
  
  log_message("========================================", log_file)
  log_message("WINDOW META-ANALYSIS", log_file)
  
  # Convert to data.table for speed
  dt <- as.data.table(window_results)
  
  # Per-method meta-analysis with fixed column consistency
  method_meta <- dt[, {
    n_windows_local <- .N
    
    # Always return the same structure regardless of conditions
    if (n_windows_local < min_windows) {
      # Single window - return values as-is
      result <- list(
        n_windows = n_windows_local,
        meta_log2FC = log2FC[1],
        meta_pvalue = pvalue[1],
        consistency = 1.0,
        peak_window = if(exists("window_name") && length(window_name) > 0) window_name[1] else NA_character_,
        peak_log2FC = log2FC[1],
        peak_center = if(exists("window_center") && length(window_center) > 0) window_center[1] else NA_real_
      )
    } else {
      # Multiple windows - perform meta-analysis
      # Handle potential NAs in padj
      padj_clean <- padj
      padj_clean[is.na(padj_clean)] <- 1
      
      # Inverse variance weighting
      weights <- 1 / pmax(padj_clean, 0.01)
      weights <- weights / sum(weights, na.rm = TRUE)
      
      # Weighted mean effect
      meta_log2FC_val <- sum(log2FC * weights, na.rm = TRUE)
      
      # Combined p-value (Fisher's method)
      pval_clean <- pvalue
      pval_clean[is.na(pval_clean)] <- 1
      pval_clean <- pmax(pval_clean, 1e-300)
      fisher_stat <- -2 * sum(log(pval_clean))
      meta_pvalue_val <- pchisq(fisher_stat, df = 2 * n_windows_local, lower.tail = FALSE)
      
      # Direction consistency
      consistency_val <- mean(sign(log2FC) == sign(median(log2FC, na.rm = TRUE)), na.rm = TRUE)
      
      # Find peak window (lowest padj)
      peak_idx <- which.min(padj_clean)
      if (length(peak_idx) == 0) peak_idx <- 1
      
      result <- list(
        n_windows = n_windows_local,
        meta_log2FC = meta_log2FC_val,
        meta_pvalue = meta_pvalue_val,
        consistency = consistency_val,
        peak_window = if(exists("window_name") && length(window_name) >= peak_idx) window_name[peak_idx] else NA_character_,
        peak_log2FC = if(length(log2FC) >= peak_idx) log2FC[peak_idx] else NA_real_,
        peak_center = if(exists("window_center") && length(window_center) >= peak_idx) window_center[peak_idx] else NA_real_
      )
    }
    
    result
  }, by = .(taxa, method)]
  
  # Adjust p-values
  # Multiple testing correction using FDR
  method_meta[, meta_padj := p.adjust(meta_pvalue, method = "BH"), by = method]
  
  # Cross-method consensus with fixed column consistency
  consensus <- method_meta[, {
    n_methods_local <- .N
    
    # Handle potential NAs
    meta_pvalue_clean <- meta_pvalue
    meta_pvalue_clean[is.na(meta_pvalue_clean)] <- 1
    
    # Combine p-values across methods using helper function
    if (exists("combine_p_values", mode = "function")) {
      consensus_pvalue_val <- combine_p_values(meta_pvalue_clean, method = "stouffer")
    } else {
      # Fallback if combine_p_values is not available
      z_scores <- qnorm(1 - meta_pvalue_clean)
      z_scores[!is.finite(z_scores)] <- 0
      combined_z <- sum(z_scores) / sqrt(length(z_scores))
      consensus_pvalue_val <- pnorm(combined_z, lower.tail = FALSE)
    }
    
    # Average effect size
    consensus_log2FC_val <- mean(meta_log2FC, na.rm = TRUE)
    
    # Agreement score
    agreement_val <- mean(consistency, na.rm = TRUE)
    
    # Count significant methods
    meta_padj_clean <- meta_padj
    meta_padj_clean[is.na(meta_padj_clean)] <- 1
    n_sig_relaxed_val <- sum(meta_padj_clean < 0.1)
    n_sig_standard_val <- sum(meta_padj_clean < 0.05)
    
    # Return consistent structure
    list(
      n_methods = n_methods_local,
      consensus_log2FC = consensus_log2FC_val,
      consensus_pvalue = consensus_pvalue_val,
      agreement_score = agreement_val,
      n_sig_relaxed = n_sig_relaxed_val,
      n_sig_standard = n_sig_standard_val
    )
  }, by = taxa]
  
  # Adjust consensus p-values
  # Multiple testing correction using FDR
  consensus[, consensus_padj := p.adjust(consensus_pvalue, method = "BH")]
  
  # Classify incorporators
  consensus[, is_incorporator := {
    # Handle NAs properly
    padj_check <- !is.na(consensus_padj) & consensus_padj < 0.05
    log2fc_check <- !is.na(consensus_log2FC) &
      abs(consensus_log2FC) > BIOLOGICAL_PARAMS$statistics$log2fc_threshold
    agreement_check <- !is.na(agreement_score) & agreement_score > 0.5
    
    padj_check & log2fc_check & agreement_check
  }]
  
  # Assign confidence levels using base R to avoid case_when issues
  consensus[, confidence_level := {
    conf <- rep("not_significant", .N)
    
    # High confidence
    high_mask <- !is.na(consensus_padj) & consensus_padj < 0.01 &
      !is.na(agreement_score) & agreement_score > 0.8
    conf[high_mask] <- "high"
    
    # Medium confidence
    med_mask <- !is.na(consensus_padj) & consensus_padj < 0.05 &
      !is.na(agreement_score) & agreement_score > 0.6 &
      conf != "high"
    conf[med_mask] <- "medium"
    
    # Low confidence
    low_mask <- !is.na(consensus_padj) & consensus_padj < 0.1 &
      !is.na(agreement_score) & agreement_score > 0.5 &
      conf != "high" & conf != "medium"
    conf[low_mask] <- "low"
    
    conf
  }]
  
  # Convert back to data.frame
  setDF(method_meta)
  setDF(consensus)
  
  log_message(sprintf("Meta-analysis complete: %d taxa analyzed, %d incorporators identified",
                      nrow(consensus),
                      sum(consensus$is_incorporator, na.rm = TRUE)),
              log_file, "SUCCESS")
  
  return(list(
    method_level = method_meta,
    consensus = consensus
  ))
}


# ============================================================================
# DENSITY STUDY ANALYSIS
# ============================================================================


#' Analyze Density Study
#'
#' @description [Function description]
#' @export
analyze_density_study <- function(physeq, sample_info, study_name, study_metadata,
                                  substrate_specific = FALSE,
                                  current_substrate = NULL,
                                  heavy_threshold = NULL) {
  
  if (is.null(heavy_threshold)) {
    heavy_threshold <- BIOLOGICAL_PARAMS$density$heavy_threshold
  }
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_name, "_density_analysis.log"))
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message("DENSITY GRADIENT ANALYSIS", log_file)
  log_message("========================================", log_file)
  
  if (substrate_specific && !is.null(current_substrate)) {
    log_message(sprintf("Substrate: %s", current_substrate), log_file)
  }
  
  meta <- sample_info$metadata
  densities <- parse_density_values(meta$gradient_pos_density)
  
  # FIX: Handle NAs in density values
  valid_density <- !is.na(densities) & densities > 0
  
  if (sum(valid_density) == 0) {
    log_message("No valid density values", log_file, "ERROR")
    return(NULL)
  }
  
  # Get density range from valid values only
  density_range <- range(densities[valid_density], na.rm = TRUE)
  log_message(sprintf("Density range: %.3f - %.3f g/ml",
                      density_range[1], density_range[2]), log_file)
  log_message(sprintf("Heavy threshold: %.3f g/ml", heavy_threshold), log_file)
  
  results_list <- list()
  
  # Detect study type for this specific analysis
  study_type <- detect_study_type(physeq, sample_info, log_file)
  
  # 1. HEAVY FRACTION ANALYSIS
  # FIX: Properly handle NAs in all logical operations
  
  # Create masks with explicit NA handling
  density_mask <- !is.na(densities) & densities >= heavy_threshold
  labeled_mask <- sample_info$labeled
  unlabeled_mask <- sample_info$unlabeled
  
  # Replace any NAs with FALSE
  density_mask[is.na(density_mask)] <- FALSE
  labeled_mask[is.na(labeled_mask)] <- FALSE
  unlabeled_mask[is.na(unlabeled_mask)] <- FALSE
  
  # Now safely combine
  heavy_labeled <- density_mask & labeled_mask
  heavy_unlabeled <- density_mask & unlabeled_mask
  
  n_heavy_labeled <- sum(heavy_labeled)
  n_heavy_unlabeled <- sum(heavy_unlabeled)
  
  log_message(sprintf("Heavy fraction: %d labeled, %d unlabeled samples",
                      n_heavy_labeled, n_heavy_unlabeled), log_file)
  
  if (n_heavy_labeled >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group &&
      n_heavy_unlabeled >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
    
    log_message("Analyzing heavy fraction...", log_file)
    
    heavy_results <- run_da_analysis(
      physeq,
      heavy_labeled,
      heavy_unlabeled,
      paste0(study_name, "_heavy"),
      study_type
    )
    
    if (!is.null(heavy_results)) {
      heavy_results$comparison_type <- "heavy_fraction"
      heavy_results$density_threshold <- heavy_threshold
      heavy_results$density_range <- paste(density_range, collapse = "-")
      
      if (substrate_specific && !is.null(current_substrate)) {
        heavy_results$substrate_analyzed <- current_substrate
      }
      
      results_list$heavy_fraction <- heavy_results
      
      log_message(sprintf("  Heavy fraction: %d significant taxa (p<0.05)",
                          sum(heavy_results$padj < 0.05, na.rm = TRUE)),
                  log_file, "SUCCESS")
    }
  } else {
    log_message("  Insufficient samples in heavy fraction", log_file, "WARNING")
  }
  
  # 2. SLIDING WINDOW ANALYSIS
  # Count valid density values for windowing
  n_valid_for_windows <- sum(valid_density)
  
  if (n_valid_for_windows >= 10) {
    log_message("Running sliding window analysis...", log_file)
    
    # Select window strategy
    strategy <- select_window_strategy(
      physeq, sample_info, study_type,
      "gradient_pos_density", log_file
    )
    
    if (!is.null(strategy)) {
      window_results <- tryCatch({
        run_windowed_da_analysis(
          physeq, sample_info, study_type,
          "gradient_pos_density", strategy,
          study_name, substrate_specific, current_substrate,
          n_cores = if(strategy$n_valid_samples > 30) 4 else 1,
          log_file
        )
      }, error = function(e) {
        log_message(sprintf("  Window analysis failed: %s", e$message),
                    log_file, "ERROR")
        NULL
      })
      
      if (!is.null(window_results)) {
        # Perform meta-analysis
        meta_results <- tryCatch({
          perform_window_meta_analysis(
            window_results$combined_results,
            min_windows = 2,
            log_file
          )
        }, error = function(e) {
          log_message(sprintf("  Meta-analysis failed: %s", e$message),
                      log_file, "ERROR")
          NULL
        })
        
        if (!is.null(window_results$combined_results)) {
          results_list$sliding_window <- window_results
        }
        
        if (!is.null(meta_results)) {
          results_list$window_meta_analysis <- meta_results
          
          log_message(sprintf("  Window analysis: %d incorporators identified",
                              sum(meta_results$consensus$is_incorporator, na.rm = TRUE)),
                      log_file, "SUCCESS")
        }
      }
    }
  } else {
    log_message("  Insufficient density values for sliding windows", log_file, "WARNING")
  }
  
  
  return(results_list)
}

# ============================================================================
# BINARY STUDY ANALYSIS
# ============================================================================

#' Analyze Binary Study
#'
#' @description [Function description]
#' @export
analyze_binary_study <- function(physeq, sample_info, study_name, study_metadata,
                                 substrate_specific = FALSE,
                                 current_substrate = NULL) {
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_name, "_binary_analysis.log"))
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message("BINARY GRADIENT ANALYSIS", log_file)
  log_message("========================================", log_file)
  
  if (substrate_specific && !is.null(current_substrate)) {
    log_message(sprintf("Substrate: %s", current_substrate), log_file)
  }
  
  meta <- sample_info$metadata
  
  # Get positions
  sip_samples <- sample_info$labeled | sample_info$unlabeled
  positions <- meta$gradient_position[sip_samples & !is.na(meta$gradient_position)]
  unique_positions <- sort(unique(positions[positions > 0]))
  
  results_list <- list()
  
  # Detect study type
  study_type <- detect_study_type(physeq, sample_info, log_file)
  
  if (length(unique_positions) == 0) {
    log_message("No valid positions found", log_file, "ERROR")
    return(NULL)
  }
  
  # Analyze heavy position (typically position 1 or lowest number)
  heavy_position <- unique_positions[1]
  
  at_heavy <- meta$gradient_position == heavy_position & !is.na(meta$gradient_position)
  labeled_heavy <- sample_info$labeled & at_heavy
  unlabeled_heavy <- sample_info$unlabeled & at_heavy
  
  log_message(sprintf("Heavy position %d: %d labeled, %d unlabeled",
                      heavy_position, sum(labeled_heavy), sum(unlabeled_heavy)),
              log_file)
  
  if (sum(labeled_heavy) >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group &&
      sum(unlabeled_heavy) >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
    
    da_results <- run_da_analysis(
      physeq,
      labeled_heavy,
      unlabeled_heavy,
      paste0(study_name, "_binary_heavy"),
      study_type
    )
    
    if (!is.null(da_results)) {
      da_results$comparison_type <- "binary_heavy"
      da_results$gradient_position <- heavy_position
      
      if (substrate_specific && !is.null(current_substrate)) {
        da_results$substrate_analyzed <- current_substrate
      }
      
      results_list$binary_heavy <- da_results
      
      log_message(sprintf("  Binary comparison: %d significant taxa (p<0.05)",
                          sum(da_results$padj < 0.05, na.rm = TRUE)),
                  log_file, "SUCCESS")
    }
  } else {
    log_message("  Insufficient samples for binary comparison", log_file, "WARNING")
  }
  
  # If there are 2 positions, also analyze light
  if (length(unique_positions) == 2) {
    light_position <- unique_positions[2]
    
    at_light <- meta$gradient_position == light_position & !is.na(meta$gradient_position)
    labeled_light <- sample_info$labeled & at_light
    unlabeled_light <- sample_info$unlabeled & at_light
    
    if (sum(labeled_light) >= 2 && sum(unlabeled_light) >= 2) {
      da_light <- run_da_analysis(
        physeq,
        labeled_light,
        unlabeled_light,
        paste0(study_name, "_binary_light"),
        study_type
      )
      
      if (!is.null(da_light)) {
        da_light$comparison_type <- "binary_light"
        da_light$gradient_position <- light_position
        
        if (substrate_specific && !is.null(current_substrate)) {
          da_light$substrate_analyzed <- current_substrate
        }
        
        results_list$binary_light <- da_light
      }
    }
  }
  
  return(results_list)
}

# ============================================================================
# MULTIFRACTION STUDY ANALYSIS
# ============================================================================

#' Analyze Multifraction Study
#'
#' @description [Function description]
#' @export
analyze_multifraction_study <- function(physeq, sample_info, study_name, study_metadata,
                                        substrate_specific = FALSE,
                                        current_substrate = NULL) {
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_name, "_multifraction_analysis.log"))
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message("MULTIFRACTION GRADIENT ANALYSIS", log_file)
  log_message("========================================", log_file)
  
  if (substrate_specific && !is.null(current_substrate)) {
    log_message(sprintf("Substrate: %s", current_substrate), log_file)
  }
  
  meta <- sample_info$metadata
  
  # Get positions
  sip_samples <- sample_info$labeled | sample_info$unlabeled
  positions <- meta$gradient_position[sip_samples & !is.na(meta$gradient_position)]
  unique_positions <- sort(unique(positions[positions > 0]))
  
  n_positions <- length(unique_positions)
  log_message(sprintf("Found %d gradient positions", n_positions), log_file)
  
  results_list <- list()
  
  # Detect study type
  study_type <- detect_study_type(physeq, sample_info, log_file)
  
  # 1. HEAVY POOL ANALYSIS
  # Define heavy as top 33% of positions
  heavy_cutoff <- unique_positions[ceiling(n_positions * 0.33)]
  heavy_positions <- unique_positions[unique_positions <= heavy_cutoff]
  
  log_message(sprintf("Heavy positions: %s",
                      paste(heavy_positions, collapse = ", ")), log_file)
  
  in_heavy <- meta$gradient_position %in% heavy_positions
  heavy_labeled <- sample_info$labeled & in_heavy
  heavy_unlabeled <- sample_info$unlabeled & in_heavy
  
  if (sum(heavy_labeled) >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group &&
      sum(heavy_unlabeled) >= BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
    
    heavy_results <- run_da_analysis(
      physeq,
      heavy_labeled,
      heavy_unlabeled,
      paste0(study_name, "_heavy_pool"),
      study_type
    )
    
    if (!is.null(heavy_results)) {
      heavy_results$comparison_type <- "heavy_pool"
      heavy_results$positions_included <- paste(heavy_positions, collapse = ",")
      
      if (substrate_specific && !is.null(current_substrate)) {
        heavy_results$substrate_analyzed <- current_substrate
      }
      
      results_list$heavy_pool <- heavy_results
      
      log_message(sprintf("  Heavy pool: %d significant taxa (p<0.05)",
                          sum(heavy_results$padj < 0.05, na.rm = TRUE)),
                  log_file, "SUCCESS")
    }
  }
  
  # 2. POSITION-BASED SLIDING WINDOWS
  if (n_positions >= 3) {
    log_message("Running position-based window analysis...", log_file)
    
    # Select window strategy
    strategy <- select_window_strategy(
      physeq, sample_info, study_type,
      "gradient_position", log_file
    )
    
    if (!is.null(strategy)) {
      window_results <- run_windowed_da_analysis(
        physeq, sample_info, study_type,
        "gradient_position", strategy,
        study_name, substrate_specific, current_substrate,
        n_cores = if(n_positions > 10) 4 else 1,
        log_file
      )
      
      if (!is.null(window_results)) {
        # Perform meta-analysis
        meta_results <- perform_window_meta_analysis(
          window_results$combined_results,
          min_windows = 2,
          log_file
        )
        
        results_list$sliding_window <- window_results
        results_list$window_meta_analysis <- meta_results
        
        if (!is.null(meta_results)) {
          log_message(sprintf("  Window analysis: %d incorporators identified",
                              sum(meta_results$consensus$is_incorporator, na.rm = TRUE)),
                      log_file, "SUCCESS")
        }
      }
    }
  }
  
  return(results_list)
}

# ============================================================================
# MAIN STUDY ANALYSIS COORDINATOR
# ============================================================================

#' Analyze Study
#'
#' @description [Function description]
#' @export
analyze_study <- function(physeq, study_id,
                          run_substrate_specific = TRUE,
                          heavy_threshold = NULL) {
  
  if (is.null(heavy_threshold)) {
    heavy_threshold <- BIOLOGICAL_PARAMS$density$heavy_threshold
  }
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_id, "_master.log"))
  ensure_dir_for_file(log_file)
  
  log_message("\n" %+% paste(rep("=", 60), collapse = ""), log_file)
  log_message(sprintf("ANALYZING STUDY: %s", study_id), log_file)
  log_message(paste(rep("=", 60), collapse = ""), log_file)
  
  # Sample identification
  sample_info <- identify_labeled_unlabeled_samples(physeq, log_file)
  
  if (!sample_info$has_both) {
    log_message("Study lacks both labeled and unlabeled samples", log_file, "ERROR")
    return(NULL)
  }
  
  # QC checks
  qc_results <- perform_qc_checks(physeq, sample_info, study_id, log_file)
  
  if (!qc_results$passed) {
    log_message("Study failed QC checks", log_file, "ERROR")
    return(NULL)
  }
  
  # Detect substrate design
  substrate_design <- detect_substrate_design(physeq, sample_info, log_file)
  
  # Detect study type
  study_type <- detect_study_type(physeq, sample_info, log_file)
  
  if (is.null(study_type)) {
    log_message("No valid gradient information - skipping", log_file, "ERROR")
    return(NULL)
  }
  
  # Decide on analysis approach
  if (run_substrate_specific &&
      substrate_design$type == "shared_control" &&
      substrate_design$n_substrates > 1) {
    
    log_message("Running SUBSTRATE-SPECIFIC analysis", log_file)
    results_list <- analyze_study_by_substrate(
      physeq, study_id, sample_info, substrate_design,
      study_type, heavy_threshold
    )
    
  } else {
    log_message("Running STANDARD analysis", log_file)
    
    # Extract metadata
    study_metadata <- extract_comprehensive_metadata(physeq, sample_info, study_id)
    study_metadata$study_type <- study_type$type
    study_metadata$qc_warnings <- length(qc_results$warnings)
    
    # Run appropriate analysis
    results_list <- NULL
    
    if (study_type$type == "density") {
      results_list <- analyze_density_study(
        physeq, sample_info, study_id, study_metadata,
        heavy_threshold = heavy_threshold
      )
    } else if (study_type$type == "binary") {
      results_list <- analyze_binary_study(
        physeq, sample_info, study_id, study_metadata
      )
    } else if (study_type$type == "multifraction") {
      results_list <- analyze_multifraction_study(
        physeq, sample_info, study_id, study_metadata
      )
    }
    
    # Add metadata to results
    if (!is.null(results_list) && length(results_list) > 0) {
      results_list <- add_metadata_to_results(
        results_list, study_metadata, physeq, study_id
      )
    }
  }
  
  # Generate QC report
  generate_qc_report(physeq, sample_info, results_list, study_id)
  
  log_message(sprintf("Analysis complete for %s", study_id), log_file, "SUCCESS")
  
  return(results_list)
}

# ============================================================================
# SUBSTRATE-SPECIFIC ANALYSIS
# ============================================================================

#' Analyze Study By Substrate
#'
#' @description [Function description]
#' @export
analyze_study_by_substrate <- function(physeq, study_id, sample_info,
                                       substrate_design, study_type,
                                       heavy_threshold) {
  
  log_file <- file.path(main_dir, "debug_logs",
                        paste0(study_id, "_substrate_analysis.log"))
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message("SUBSTRATE-SPECIFIC ANALYSIS", log_file)
  log_message("========================================", log_file)
  
  meta <- sample_info$metadata
  substrates <- substrate_design$substrates
  
  all_substrate_results <- list()
  
  for (substrate in substrates) {
    log_message(sprintf("\nAnalyzing substrate: %s", substrate), log_file)
    log_message(paste(rep("-", 40), collapse = ""), log_file)
    
    # Filter for this substrate
    if ("isotopolog" %in% colnames(meta)) {
      substrate_samples <- (meta$isotopolog == substrate & sample_info$labeled) |
        sample_info$unlabeled
    } else {
      log_message("No isotopolog column for filtering", log_file, "ERROR")
      next
    }
    
    n_substrate_labeled <- sum(meta$isotopolog == substrate & sample_info$labeled, na.rm = TRUE)
    n_controls <- sum(sample_info$unlabeled)
    
    log_message(sprintf("  Samples: %d labeled, %d controls",
                        n_substrate_labeled, n_controls), log_file)
    
    if (n_substrate_labeled < BIOLOGICAL_PARAMS$statistics$min_samples_per_group) {
      log_message("  Skipping: insufficient labeled samples", log_file, "WARNING")
      next
    }
    
    # Subset phyloseq
    physeq_substrate <- prune_samples(substrate_samples, physeq)
    
    # Create substrate-specific sample info
    substrate_sample_info <- list(
      labeled = sample_info$labeled[substrate_samples],
      unlabeled = sample_info$unlabeled[substrate_samples],
      metadata = as.data.frame(sample_data(physeq_substrate)),
      n_labeled = n_substrate_labeled,
      n_unlabeled = n_controls,
      isotopes = sample_info$isotopes,
      substrates = substrate,
      has_both = TRUE,
      detection_confidence = sample_info$detection_confidence
    )
    
    # Re-detect study type for this subset
    substrate_study_type <- detect_study_type(physeq_substrate, substrate_sample_info, log_file)
    
    if (is.null(substrate_study_type)) {
      log_message("  Skipping: no valid gradient information", log_file, "WARNING")
      next
    }
    
    # Extract metadata
    substrate_study_name <- paste0(study_id, "_", gsub("[^A-Za-z0-9]", "", substrate))
    study_metadata <- extract_comprehensive_metadata(
      physeq_substrate, substrate_sample_info, substrate_study_name,
      substrate_specific = TRUE, current_substrate = substrate
    )
    study_metadata$study_type <- substrate_study_type$type
    
    # Run appropriate analysis
    substrate_results <- NULL
    
    if (substrate_study_type$type == "density") {
      substrate_results <- analyze_density_study(
        physeq_substrate, substrate_sample_info,
        substrate_study_name, study_metadata,
        substrate_specific = TRUE,
        current_substrate = substrate,
        heavy_threshold = heavy_threshold
      )
    } else if (substrate_study_type$type == "binary") {
      substrate_results <- analyze_binary_study(
        physeq_substrate, substrate_sample_info,
        substrate_study_name, study_metadata,
        substrate_specific = TRUE,
        current_substrate = substrate
      )
    } else if (substrate_study_type$type == "multifraction") {
      substrate_results <- analyze_multifraction_study(
        physeq_substrate, substrate_sample_info,
        substrate_study_name, study_metadata,
        substrate_specific = TRUE,
        current_substrate = substrate
      )
    }
    
    if (!is.null(substrate_results) && length(substrate_results) > 0) {
      # Add metadata
      substrate_results <- add_metadata_to_results(
        substrate_results, study_metadata, physeq_substrate, substrate_study_name
      )
      
      # Mark as substrate-specific
      for (analysis_type in names(substrate_results)) {
        if (is.data.frame(substrate_results[[analysis_type]])) {
          substrate_results[[analysis_type]]$substrate_analyzed <- substrate
        }
      }
      
      all_substrate_results[[substrate]] <- substrate_results
    }
  }
  
  # Flatten results
  if (length(all_substrate_results) > 0) {
    flattened_results <- flatten_substrate_results(all_substrate_results, study_id)
    return(flattened_results)
  }
  
  return(NULL)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Add Metadata To Results
#'
#' @description [Function description]
#' @export
add_metadata_to_results <- function(results_list, study_metadata, physeq, study_id) {
  
  # Get taxonomy information
  tax_info <- get_cached_taxonomy(physeq, study_id)
  
  add_metadata <- function(df) {
    if (is.data.frame(df)) {
      # Add taxonomy
      if (nrow(tax_info) > 0 && "taxa" %in% colnames(df)) {
        df <- merge(df, tax_info, by = "taxa", all.x = TRUE)
      }
      
      # Add study metadata
      for (meta_name in names(study_metadata)) {
        if (!meta_name %in% colnames(df)) {
          df[[meta_name]] <- study_metadata[[meta_name]]
        }
      }
    }
    return(df)
  }
  
  # Apply to all results
  enhanced_results <- list()
  
  for (name in names(results_list)) {
    x <- results_list[[name]]
    
    if (is.data.frame(x)) {
      enhanced_results[[name]] <- add_metadata(x)
    } else if (is.list(x)) {
      # Handle nested results
      if ("combined_results" %in% names(x) && is.data.frame(x$combined_results)) {
        x$combined_results <- add_metadata(x$combined_results)
      }
      if ("consensus" %in% names(x) && is.data.frame(x$consensus)) {
        x$consensus <- add_metadata(x$consensus)
      }
      enhanced_results[[name]] <- x
    } else {
      enhanced_results[[name]] <- x
    }
  }
  
  return(enhanced_results)
}

#' Flatten Substrate Results
#'
#' @description [Function description]
#' @export
flatten_substrate_results <- function(substrate_results_list, study_id) {
  flattened <- list()
  
  for (substrate in names(substrate_results_list)) {
    substrate_results <- substrate_results_list[[substrate]]
    substrate_clean <- gsub("[^A-Za-z0-9]", "", substrate)
    
    for (analysis_type in names(substrate_results)) {
      if (analysis_type == "sliding_window") {
        analysis_name <- paste0(substrate_clean, "_sliding_window")
      } else if (analysis_type == "window_meta_analysis") {
        analysis_name <- paste0(substrate_clean, "_meta_analysis")
      } else {
        analysis_name <- paste0(substrate_clean, "_", analysis_type)
      }
      
      flattened[[analysis_name]] <- substrate_results[[analysis_type]]
    }
  }
  
  return(flattened)
}

#' Generate Qc Report
#'
#' @description [Function description]
#' @export
generate_qc_report <- function(physeq, sample_info, results, study_id) {
  # Create a simple QC report (placeholder - can be expanded)
  qc_file <- file.path(main_dir, "qc_reports", paste0(study_id, "_qc_report.txt"))
  ensure_dir_for_file(qc_file)
  
  sink(qc_file)
  cat("QC REPORT:", study_id, "\n")
  cat("Date:", Sys.Date(), "\n\n")
  cat("Samples:", sample_info$n_labeled, "labeled,", sample_info$n_unlabeled, "unlabeled\n")
  cat("Taxa:", ntaxa(physeq), "\n")
  cat("Library size range:", range(sample_sums(physeq)), "\n")
  
  if (!is.null(results)) {
    cat("\nAnalysis Results:\n")
    for (name in names(results)) {
      if (is.data.frame(results[[name]])) {
        cat("  ", name, ":", nrow(results[[name]]), "results\n")
      }
    }
  }
  
  sink()
}

# ============================================================================
# SCHEMA ENFORCEMENT AND VALIDATION
# ============================================================================

#' Enforce Master Schema
#'
#' @description [Function description]
#' @export
enforce_master_schema <- function(df, studies) {
  
  # Define required columns
  required_cols <- c(
    # Core results
    "taxa", "method", "log2FC", "pvalue", "padj",
    
    # Study metadata
    "study_id", "analysis_type", "comparison_type",
    "n_labeled", "n_unlabeled",
    
    # Gradient information
    "density_threshold", "gradient_position", "position_label",
    "window_name", "window_center", "window_type",
    
    # Substrate information
    "substrate_analyzed", "isotopolog", "substrates_list",
    
    # Environmental metadata
    "DOI_URL", "Bioproject_ID", "env_biome", "env_feature",
    "env_material", "environment_label",
    
    # Study characteristics
    "study_type", "n_fractions", "has_density", "primary_isotope",
    
    # Taxonomy
    "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species",
    
    # Quality metrics
    "detection_confidence", "qc_warnings", "pipeline_version"
  )
  
  # Add missing columns with appropriate defaults
  for (col in required_cols) {
    if (!col %in% names(df)) {
      if (col %in% c("n_labeled", "n_unlabeled", "n_fractions", "qc_warnings")) {
        df[[col]] <- NA_integer_
      } else if (col %in% c("log2FC", "pvalue", "padj", "density_threshold",
                            "window_center", "gradient_position")) {
        df[[col]] <- NA_real_
      } else if (col %in% c("has_density")) {
        df[[col]] <- NA
      } else {
        df[[col]] <- NA_character_
      }
    }
  }
  
  # Clean duplicate columns
  df <- df %>%
    select(-ends_with(".y")) %>%
    rename_with(~str_remove(.x, "\\.x$"), ends_with(".x"))
  
  # Ensure critical columns have proper types
  df <- df %>%
    mutate(
      log2FC = as.numeric(log2FC),
      pvalue = as.numeric(pvalue),
      padj = as.numeric(padj),
      n_labeled = suppressWarnings(as.integer(n_labeled)),
      n_unlabeled = suppressWarnings(as.integer(n_unlabeled))
    )
  
  return(df)
}

# ============================================================================
# MASTER RESULTS COMBINATION
# ============================================================================

#' Combine All Results
#'
#' @description [Function description]
#' @export
combine_all_results <- function(results_list, studies) {
  
  log_file <- file.path(main_dir, "debug_logs", "results_combination.log")
  ensure_dir_for_file(log_file)
  
  log_message("========================================", log_file)
  log_message("COMBINING ALL RESULTS", log_file)
  log_message("========================================", log_file)
  
  master_list <- list()
  
  for (study_name in names(results_list)) {
    study_results <- results_list[[study_name]]
    
    if (is.null(study_results) || length(study_results) == 0) {
      log_message(sprintf("Skipping %s: no results", study_name), log_file, "WARNING")
      next
    }
    
    # DON'T access phyloseq objects - this triggers date parsing
    # Just process the results directly
    
    flattened_results <- list()
    
    for (analysis_type in names(study_results)) {
      if (is.data.frame(study_results[[analysis_type]])) {
        df <- study_results[[analysis_type]]
        df$study_id <- study_name
        df$analysis_type <- analysis_type
        flattened_results[[analysis_type]] <- df
        
      } else if (is.list(study_results[[analysis_type]])) {
        if ("combined_results" %in% names(study_results[[analysis_type]])) {
          df <- study_results[[analysis_type]]$combined_results
          if (is.data.frame(df)) {
            df$study_id <- study_name
            df$analysis_type <- analysis_type
            flattened_results[[paste0(analysis_type, "_results")]] <- df
          }
        }
        if ("consensus" %in% names(study_results[[analysis_type]])) {
          df <- study_results[[analysis_type]]$consensus
          if (is.data.frame(df)) {
            df$study_id <- study_name
            df$analysis_type <- analysis_type
            if ("consensus_log2FC" %in% names(df)) df$log2FC <- df$consensus_log2FC
            if ("consensus_pvalue" %in% names(df)) df$pvalue <- df$consensus_pvalue
            if ("consensus_padj" %in% names(df)) df$padj <- df$consensus_padj
            flattened_results[[paste0(analysis_type, "_consensus")]] <- df
          }
        }
      }
    }
    
    if (length(flattened_results) > 0) {
      study_combined <- bind_rows(flattened_results)
      master_list[[study_name]] <- study_combined
      log_message(sprintf("  %s: %d records combined", study_name, nrow(study_combined)),
                  log_file, "SUCCESS")
    }
  }
  
  if (length(master_list) == 0) {
    log_message("No valid results to combine", log_file, "ERROR")
    return(data.frame())
  }
  
  master_df <- bind_rows(master_list)
  
  log_message(sprintf("Master table created: %d total records", nrow(master_df)),
              log_file, "SUCCESS")
  
  return(master_df)
}

# ============================================================================
# SUMMARY GENERATION
# ============================================================================

#' Generate Analysis Summaries
#'
#' @description [Function description]
#' @export
generate_analysis_summaries <- function(master_table, output_dir) {
  
  if (nrow(master_table) == 0) {
    warning("Empty master table - cannot generate summaries")
    return(list())
  }
  
  summaries <- list()
  
  # Overall summary
  summaries$overall <- master_table %>%
    summarise(
      n_studies = n_distinct(study_id, na.rm = TRUE),
      n_taxa_total = n_distinct(taxa, na.rm = TRUE),
      n_methods = n_distinct(method, na.rm = TRUE),
      n_significant_005 = sum(padj < 0.05, na.rm = TRUE),
      n_significant_01 = sum(padj < 0.1, na.rm = TRUE),
      mean_abs_log2FC = mean(abs(log2FC[padj < 0.05]), na.rm = TRUE)
    )
  
  # Study summary - check for column existence
  summary_cols <- c("study_id", "taxa", "padj", "log2FC", "method")
  available_cols <- intersect(summary_cols, names(master_table))
  
  if ("study_id" %in% names(master_table)) {
    summaries$by_study <- master_table %>%
      group_by(study_id) %>%
      summarise(
        n_records = n(),
        n_analyses = if("analysis_type" %in% names(.)) n_distinct(analysis_type, na.rm = TRUE) else 1,
        n_taxa = if("taxa" %in% names(.)) n_distinct(taxa, na.rm = TRUE) else NA,
        n_sig_005 = if("padj" %in% names(.)) sum(padj < 0.05, na.rm = TRUE) else NA,
        n_sig_01 = if("padj" %in% names(.)) sum(padj < 0.1, na.rm = TRUE) else NA,
        mean_effect = if(all(c("padj", "log2FC") %in% names(.))) {
          mean(abs(log2FC[padj < 0.05]), na.rm = TRUE)
        } else NA,
        study_type = if("study_type" %in% names(.)) dplyr::first(study_type) else NA,
        has_density = if("has_density" %in% names(.)) dplyr::first(has_density) else NA,
        n_substrates = if("substrate_analyzed" %in% names(.)) {
          n_distinct(substrate_analyzed[substrate_analyzed != "all_substrates"])
        } else NA,
        .groups = "drop"
      )
  }
  
  return(summaries)
}

# ============================================================================
# EXCEL OUTPUT GENERATION
# ============================================================================

#' Save To Excel
#'
#' @description [Function description]
#' @export
save_to_excel <- function(master_table, summaries, output_dir) {
  
  wb <- createWorkbook()
  
  # Master table
  addWorksheet(wb, "Master_Results")
  writeData(wb, "Master_Results", master_table)
  
  # Summary sheets
  if ("overall" %in% names(summaries)) {
    addWorksheet(wb, "Overall_Summary")
    writeData(wb, "Overall_Summary", summaries$overall)
  }
  
  if ("by_study" %in% names(summaries)) {
    addWorksheet(wb, "Study_Summary")
    writeData(wb, "Study_Summary", summaries$by_study)
  }
  
  if ("by_method" %in% names(summaries)) {
    addWorksheet(wb, "Method_Performance")
    writeData(wb, "Method_Performance", summaries$by_method)
  }
  
  if ("by_substrate" %in% names(summaries)) {
    addWorksheet(wb, "Substrate_Analysis")
    writeData(wb, "Substrate_Analysis", summaries$by_substrate)
  }
  
  if ("top_incorporators" %in% names(summaries)) {
    addWorksheet(wb, "Top_Incorporators")
    writeData(wb, "Top_Incorporators", summaries$top_incorporators)
  }
  
  if ("cross_study" %in% names(summaries)) {
    addWorksheet(wb, "Cross_Study_Taxa")
    writeData(wb, "Cross_Study_Taxa", summaries$cross_study)
  }
  
  # Taxonomic summaries
  if ("taxonomic" %in% names(summaries)) {
    for (level in names(summaries$taxonomic)) {
      sheet_name <- paste0(level, "_Summary")
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, summaries$taxonomic[[level]])
    }
  }
  
  # Save workbook
  filename <- file.path(output_dir, "tables",
                        paste0("sip_analysis_v1.0_",
                               format(Sys.Date(), "%Y%m%d"), ".xlsx"))
  
  ensure_dir_for_file(filename)
  saveWorkbook(wb, filename, overwrite = TRUE)
  
  message(sprintf("✓ Excel report saved: %s", basename(filename)))
  
  return(filename)
}

# ============================================================================
# TEXT SUMMARY GENERATION
# ============================================================================

#' Generate Text Summary
#'
#' @description [Function description]
#' @export
generate_text_summary <- function(master_table, summaries, output_dir) {
  
  summary_file <- file.path(output_dir, "tables",
                            paste0("sip_analysis_summary_v1.0_",
                                   format(Sys.Date(), "%Y%m%d"), ".txt"))
  
  ensure_dir_for_file(summary_file)
  sink(summary_file)
  
  cat("================================================================================\n")
  cat("SIP META-ANALYSIS PIPELINE v1.0 - SUMMARY REPORT\n")
  cat("================================================================================\n\n")
  
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Pipeline version 1.0 (Optimized)\n\n")
  
  # Overall statistics
  if ("overall" %in% names(summaries)) {
    cat("OVERALL STATISTICS\n")
    cat("------------------\n")
    cat(sprintf("Studies analyzed: %d\n", summaries$overall$n_studies))
    cat(sprintf("Total unique taxa: %d\n", summaries$overall$n_taxa_total))
    cat(sprintf("Substrates analyzed: %d\n", summaries$overall$n_substrates))
    cat(sprintf("Significant taxa (p<0.05): %d\n", summaries$overall$n_significant_005))
    cat(sprintf("Mean absolute log2FC: %.2f\n", summaries$overall$mean_abs_log2FC))
    cat("\n")
  }
  
  # Study breakdown
  if ("by_study" %in% names(summaries)) {
    cat("STUDY BREAKDOWN\n")
    cat("---------------\n")
    print(summaries$by_study, n = 20)
    cat("\n")
  }
  
  # Top incorporators
  if ("top_incorporators" %in% names(summaries) &&
      nrow(summaries$top_incorporators) > 0) {
    cat("TOP 20 INCORPORATORS\n")
    cat("--------------------\n")
    print(head(summaries$top_incorporators, 20), n = 20)
    cat("\n")
  }
  
  # Cross-study consensus
  if ("cross_study" %in% names(summaries) &&
      nrow(summaries$cross_study) > 0) {
    cat("CROSS-STUDY CONSENSUS TAXA\n")
    cat("--------------------------\n")
    cat(sprintf("Taxa found in multiple studies: %d\n",
                nrow(summaries$cross_study)))
    print(head(summaries$cross_study, 10), n = 10)
    cat("\n")
  }
  
  # Method comparison
  if ("by_method" %in% names(summaries)) {
    cat("METHOD COMPARISON\n")
    cat("-----------------\n")
    print(summaries$by_method)
    cat("\n")
  }
  
  cat("================================================================================\n")
  cat("END OF REPORT\n")
  cat("================================================================================\n")
  
  sink()
  
  message(sprintf("✓ Text summary saved: %s", basename(summary_file)))
  
  return(summary_file)
}

# ============================================================================
# DUAL THRESHOLD ANALYSIS
# ============================================================================

#' Call Consensus Incorporators
#'
#' @description [Function description]
#' @export
call_consensus_incorporators <- function(relaxed_results, stringent_results,
                                         log2fc_threshold = NULL,
                                         padj_threshold = NULL) {
  
  if (is.null(log2fc_threshold)) {
    log2fc_threshold <- BIOLOGICAL_PARAMS$statistics$log2fc_threshold
  }
  
  if (is.null(padj_threshold)) {
    padj_threshold <- BIOLOGICAL_PARAMS$statistics$padj_threshold
  }
  
  # Filter relaxed results
  relaxed_sig <- relaxed_results %>%
    dplyr::filter(
      !is.na(log2FC),
      !is.na(padj),
      padj < padj_threshold,
      abs(log2FC) > log2fc_threshold
    )
  
  if (nrow(relaxed_sig) == 0) {
    return(list(
      incorporators = data.frame(),
      consensus_taxa = data.frame(),
      relaxed_only = data.frame(),
      summary = list(
        n_consensus = 0,
        n_relaxed_only = 0,
        consensus_rate = 0
      )
    ))
  }
  
  # Filter stringent results
  stringent_sig <- stringent_results %>%
    dplyr::filter(
      !is.na(log2FC),
      !is.na(padj),
      padj < padj_threshold,
      abs(log2FC) > log2fc_threshold
    )
  
  if (nrow(stringent_sig) == 0) {
    # Only relaxed results exist
    relaxed_sig$confidence <- "relaxed_only"
    
    return(list(
      incorporators = relaxed_sig,
      consensus_taxa = data.frame(),
      relaxed_only = relaxed_sig,
      summary = list(
        n_consensus = 0,
        n_relaxed_only = nrow(relaxed_sig),
        consensus_rate = 0
      )
    ))
  }
  
  # Find consensus - use explicit dplyr namespace to avoid MASS conflicts
  join_cols <- c("taxa", "study_id", "method", "substrate_analyzed")
  join_cols <- intersect(join_cols, names(relaxed_sig))
  join_cols <- intersect(join_cols, names(stringent_sig))
  
  # Prepare relaxed data for join
  relaxed_for_join <- relaxed_sig %>%
    dplyr::select(dplyr::all_of(c(join_cols, "log2FC", "padj"))) %>%
    dplyr::rename(log2FC_relaxed = log2FC, padj_relaxed = padj)
  
  # Prepare stringent data for join
  stringent_for_join <- stringent_sig %>%
    dplyr::select(dplyr::all_of(c(join_cols, "log2FC", "padj"))) %>%
    dplyr::rename(log2FC_stringent = log2FC, padj_stringent = padj)
  
  # Perform the join
  consensus <- dplyr::inner_join(
    relaxed_for_join,
    stringent_for_join,
    by = join_cols
  )
  
  if (nrow(consensus) > 0) {
    consensus$mean_log2FC <- (consensus$log2FC_relaxed + consensus$log2FC_stringent) / 2
    consensus$min_padj <- pmin(consensus$padj_relaxed, consensus$padj_stringent)
    consensus$confidence <- "high_consensus"
  }
  
  # Find relaxed-only
  relaxed_only <- dplyr::anti_join(
    relaxed_sig,
    consensus %>% dplyr::select(dplyr::all_of(join_cols)),
    by = join_cols
  )
  
  if (nrow(relaxed_only) > 0) {
    relaxed_only$confidence <- "relaxed_only"
    relaxed_only$mean_log2FC <- relaxed_only$log2FC
    relaxed_only$min_padj <- relaxed_only$padj
  }
  
  # Combine all
  all_incorporators <- dplyr::bind_rows(consensus, relaxed_only)
  
  # Summary
  summary_stats <- list(
    n_consensus = nrow(consensus),
    n_relaxed_only = nrow(relaxed_only),
    total_incorporators = nrow(all_incorporators),
    consensus_rate = if (nrow(all_incorporators) > 0) {
      nrow(consensus) / nrow(all_incorporators)
    } else 0,
    consensus_taxa = if(nrow(consensus) > 0) unique(consensus$taxa) else character()
  )
  
  return(list(
    incorporators = all_incorporators,
    consensus_taxa = consensus,
    relaxed_only = relaxed_only,
    summary = summary_stats
  ))
}

# Define output directory

resolve_out_dir <- function(out_dir = NULL, default_dir = main_dir) {
  if (!is.null(out_dir) && nzchar(out_dir)) {
    return(normalizePath(out_dir, mustWork = FALSE))
  }
  normalizePath(default_dir, mustWork = FALSE)
}


################################################################################
# SIP Meta-Analysis Pipeline v1.0
################################################################################

# ============================================================================
# MAIN PIPELINE FUNCTION
# ============================================================================

#' Run Complete Sip Analysis
#'
#' @description [Function description]
#' @export
run_complete_sip_analysis <- function(studies,
                                      output_prefix = "sip_analysis_v1.0",
                                      out_dir = NULL,
                                      heavy_threshold = NULL,
                                      run_dual_threshold = TRUE,
                                      stringent_threshold = NULL,
                                      n_cores = NULL,
                                      # NEW CHECKPOINT PARAMETERS:
                                      resume = FALSE,
                                      run_id = NULL,
                                      checkpoint_dir = NULL) {
  
  cat("\n")
  cat("████████████████████████████████████████████████████████████\n")
  cat("█  SIP META-ANALYSIS PIPELINE v1.0 (OPTIMIZED)            █\n")
  cat("████████████████████████████████████████████████████████████\n")
  cat("\n")
  
  main_dir <- resolve_out_dir(out_dir, default_dir = main_dir)
  dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
  
  # === CHECKPOINT SETUP ===
  if (is.null(run_id)) {
    run_id <- format(Sys.time(), "run_%Y%m%d_%H%M%S")
  }
  if (is.null(checkpoint_dir)) {
    checkpoint_dir <- file.path(main_dir, "checkpoints")
  }
  ensure_dir(checkpoint_dir)
  
  cat(sprintf("🔖 Run ID: %s\n", run_id))
  cat(sprintf("📁 Checkpoints: %s\n", checkpoint_dir))
  if (resume) {
    cat("🔄 RESUME MODE: Will load from existing checkpoints\n")
  }
  cat("\n")
  
  start_time <- Sys.time()
  
  # Set default thresholds
  if (is.null(heavy_threshold)) {
    heavy_threshold <- BIOLOGICAL_PARAMS$density$heavy_threshold - 0.025  # 1.700
  }
  
  if (is.null(stringent_threshold)) {
    stringent_threshold <- BIOLOGICAL_PARAMS$density$heavy_threshold  # 1.725
  }
  
  # Determine core usage
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  cat(sprintf("Configuration:\n"))
  cat(sprintf("  - Heavy threshold: %.3f g/ml\n", heavy_threshold))
  if (run_dual_threshold) {
    cat(sprintf("  - Stringent threshold: %.3f g/ml\n", stringent_threshold))
  }
  cat(sprintf("  - Parallel cores: %d\n", n_cores))
  cat(sprintf("  - Studies to analyze: %d\n", length(studies)))
  cat("\n")
  
  # FIX: Clean problematic metadata columns
  for (study_name in names(studies)) {
    meta <- as.data.frame(sample_data(studies[[study_name]]))
    problem_cols <- grep("Collection\\.Date|Date\\.|date$|_date", colnames(meta), value = TRUE)
    problem_cols <- setdiff(problem_cols, "isotopolog_incu_time")
    if (length(problem_cols) > 0) {
      for (col in problem_cols) meta[[col]] <- NULL
      sample_data(studies[[study_name]]) <- sample_data(meta)
    }
  }
  
  # ==========================================================================
  # STEP 1: Compatibility check
  # ==========================================================================
  cat("Step 1: Checking study compatibility...\n")
  
  step1_file <- file.path(checkpoint_dir, paste0(run_id, "_step01_compatibility.rds"))
  
  if (resume && file.exists(step1_file)) {
    cat("  ⏭️  Loading from checkpoint...\n")
    compatibility <- readRDS(step1_file)
  } else {
    compatibility <- check_pipeline_compatibility(studies)
    saveRDS(compatibility, step1_file)
    cat("  💾 Checkpoint saved\n")
  }
  
  if (compatibility$overall_rate == 0) {
    stop("No compatible studies found. Check data format and requirements.")
  }
  
  cat(sprintf("  → %.0f%% studies compatible\n", compatibility$overall_rate * 100))
  
  # ==========================================================================
  # STEP 2: Study overviews (no checkpoint - fast step)
  # ==========================================================================
  cat("\nStep 2: Generating study overviews...\n")
  for (study_name in names(studies)) {
    overview <- get_study_overview(studies[[study_name]], study_name, silent = TRUE)
    if (!is.null(overview)) {
      cat(sprintf("  → %s: %s (%d taxa, %d samples)\n",
                  study_name, overview$study_type, overview$n_taxa, overview$n_samples))
    }
  }
  
  # ==========================================================================
  # STEP 3: Relaxed threshold analysis
  # ==========================================================================
  cat(sprintf("\nStep 3: Running relaxed analysis (%.3f g/ml)...\n", heavy_threshold))
  
  step3_file <- file.path(checkpoint_dir, paste0(run_id, "_step03_relaxed.rds"))
  
  if (resume && file.exists(step3_file)) {
    cat("  ⏭️  Loading from checkpoint...\n")
    results_relaxed <- readRDS(step3_file)
    cat(sprintf("  ✓ Loaded %d study results\n", length(results_relaxed)))
    
  } else {
    results_relaxed <- list()
    pb <- txtProgressBar(min = 0, max = length(studies), style = 3)
    
    for (i in seq_along(studies)) {
      study_name <- names(studies)[i]
      setTxtProgressBar(pb, i)
      
      study_results <- tryCatch({
        analyze_study(
          studies[[study_name]],
          study_name,
          run_substrate_specific = TRUE,
          heavy_threshold = heavy_threshold
        )
      }, error = function(e) {
        message(sprintf("\n  Error in %s: %s", study_name, e$message))
        NULL
      })
      
      if (!is.null(study_results)) {
        results_relaxed[[study_name]] <- study_results
      }
    }
    
    close(pb)
    
    saveRDS(results_relaxed, step3_file)
    cat("\n  💾 Checkpoint saved\n")
  }
  
  cat(sprintf("  → Completed: %d/%d studies produced results\n",
              length(results_relaxed), length(studies)))
  
  # ==========================================================================
  # STEP 4: Combine relaxed results
  # ==========================================================================
  cat("\nStep 4: Combining relaxed results...\n")
  
  step4_file <- file.path(checkpoint_dir, paste0(run_id, "_step04_master_relaxed.rds"))
  
  if (resume && file.exists(step4_file)) {
    cat("  ⏭️  Loading from checkpoint...\n")
    master_relaxed <- readRDS(step4_file)
  } else {
    master_relaxed <- combine_all_results(results_relaxed, studies)
    
    if (nrow(master_relaxed) == 0) {
      stop("No data in relaxed master table")
    }
    
    saveRDS(master_relaxed, step4_file)
    cat("  💾 Checkpoint saved\n")
  }
  
  cat(sprintf("  → Combined: %d total records\n", nrow(master_relaxed)))
  
  # ==========================================================================
  # STEP 5-7: Stringent analysis (optional)
  # ==========================================================================
  results_stringent <- NULL
  master_stringent <- NULL
  consensus_results <- NULL
  
  if (run_dual_threshold) {
    
    # ========================================================================
    # STEP 5: Stringent analysis
    # ========================================================================
    cat(sprintf("\nStep 5: Running stringent analysis (%.3f g/ml)...\n",
                stringent_threshold))
    
    step5_file <- file.path(checkpoint_dir, paste0(run_id, "_step05_stringent.rds"))
    
    if (resume && file.exists(step5_file)) {
      cat("  ⏭️  Loading from checkpoint...\n")
      results_stringent <- readRDS(step5_file)
      cat(sprintf("  ✓ Loaded %d study results\n", length(results_stringent)))
      
    } else {
      results_stringent <- list()
      pb <- txtProgressBar(min = 0, max = length(studies), style = 3)
      
      for (i in seq_along(studies)) {
        study_name <- names(studies)[i]
        setTxtProgressBar(pb, i)
        
        study_results <- tryCatch({
          analyze_study(
            studies[[study_name]],
            study_name,
            run_substrate_specific = TRUE,
            heavy_threshold = stringent_threshold
          )
        }, error = function(e) {
          NULL
        })
        
        if (!is.null(study_results)) {
          results_stringent[[study_name]] <- study_results
        }
      }
      
      close(pb)
      
      saveRDS(results_stringent, step5_file)
      cat("\n  💾 Checkpoint saved\n")
    }
    
    cat(sprintf("  → Completed: %d/%d studies produced results\n",
                length(results_stringent), length(studies)))
    
    if (length(results_stringent) > 0) {
      
      # ======================================================================
      # STEP 6: Combine stringent results
      # ======================================================================
      cat("\nStep 6: Combining stringent results...\n")
      
      step6_file <- file.path(checkpoint_dir, paste0(run_id, "_step06_master_stringent.rds"))
      
      if (resume && file.exists(step6_file)) {
        cat("  ⏭️  Loading from checkpoint...\n")
        master_stringent <- readRDS(step6_file)
      } else {
        master_stringent <- combine_all_results(results_stringent, studies)
        saveRDS(master_stringent, step6_file)
        cat("  💾 Checkpoint saved\n")
      }
      
      cat(sprintf("  → Combined: %d total records\n", nrow(master_stringent)))
      
      # ======================================================================
      # STEP 7: Consensus calling
      # ======================================================================
      cat("\nStep 7: Identifying consensus incorporators...\n")
      
      step7_file <- file.path(checkpoint_dir, paste0(run_id, "_step07_consensus.rds"))
      
      if (resume && file.exists(step7_file)) {
        cat("  ⏭️  Loading from checkpoint...\n")
        consensus_results <- readRDS(step7_file)
      } else {
        consensus_results <- call_consensus_incorporators(
          master_relaxed,
          master_stringent
        )
        saveRDS(consensus_results, step7_file)
        cat("  💾 Checkpoint saved\n")
      }
      
      cat(sprintf("  → Consensus: %d high-confidence, %d relaxed-only\n",
                  consensus_results$summary$n_consensus,
                  consensus_results$summary$n_relaxed_only))
    }
  }
  
  # ==========================================================================
  # STEP 8: Generate summaries
  # ==========================================================================
  final_step <- if (run_dual_threshold) 8 else 5
  cat(sprintf("\nStep %d: Generating summaries...\n", final_step))
  
  step8_file <- file.path(checkpoint_dir, paste0(run_id, "_step08_summaries.rds"))
  
  if (resume && file.exists(step8_file)) {
    cat("  ⏭️  Loading from checkpoint...\n")
    summaries_relaxed <- readRDS(step8_file)
  } else {
    summaries_relaxed <- generate_analysis_summaries(master_relaxed, main_dir)
    saveRDS(summaries_relaxed, step8_file)
    cat("  💾 Checkpoint saved\n")
  }
  
  # ==========================================================================
  # STEP 9: Save results (final output - no checkpoint needed)
  # ==========================================================================
  cat(sprintf("\nStep %d: Saving results...\n", final_step + 1))
  
  saved_files <- save_results(
    results_relaxed = results_relaxed,
    master_relaxed = master_relaxed,
    results_stringent = results_stringent,
    master_stringent = master_stringent,
    consensus_results = consensus_results,
    summaries = summaries_relaxed,
    output_prefix = output_prefix,
    output_dir = main_dir
  )
  
  # ==========================================================================
  # FINAL SUMMARY
  # ==========================================================================
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════\n")
  cat("ANALYSIS COMPLETE\n")
  cat("════════════════════════════════════════════════════════════\n")
  
  cat(sprintf("Runtime: %.1f minutes\n", as.numeric(runtime)))
  cat(sprintf("Studies processed: %d\n", length(results_relaxed)))
  
  # Results summary
  cat("\nRELAXED RESULTS:\n")
  cat(sprintf("  - Total records: %d\n", nrow(master_relaxed)))
  cat(sprintf("  - Significant taxa (p<0.05): %d\n",
              sum(master_relaxed$padj < 0.05, na.rm = TRUE)))
  
  if (!is.null(master_stringent)) {
    cat("\nSTRINGENT RESULTS:\n")
    cat(sprintf("  - Total records: %d\n", nrow(master_stringent)))
    cat(sprintf("  - Significant taxa (p<0.05): %d\n",
                sum(master_stringent$padj < 0.05, na.rm = TRUE)))
  }
  
  if (!is.null(consensus_results)) {
    cat("\nCONSENSUS INCORPORATORS:\n")
    cat(sprintf("  - High-confidence: %d\n", consensus_results$summary$n_consensus))
    cat(sprintf("  - Relaxed-only: %d\n", consensus_results$summary$n_relaxed_only))
    cat(sprintf("  - Consensus rate: %.1f%%\n",
                consensus_results$summary$consensus_rate * 100))
  }
  
  cat(sprintf("\n📁 Results saved to: %s\n", main_dir))
  cat("\nOutput files:\n")
  for (file in names(saved_files)) {
    if (!is.null(saved_files[[file]])) {
      cat(sprintf("  • %s: %s\n", file, basename(saved_files[[file]])))
    }
  }
  
  # Clean up checkpoint message
  cat(sprintf("\n🔖 Checkpoints saved in: %s\n", checkpoint_dir))
  cat("   To resume a failed run, use: resume = TRUE, run_id = \"", run_id, "\"\n", sep = "")
  
  return(list(
    results_relaxed = results_relaxed,
    master_relaxed = master_relaxed,
    results_stringent = results_stringent,
    master_stringent = master_stringent,
    consensus = consensus_results,
    files = saved_files,
    runtime = runtime,
    summaries = summaries_relaxed,
    run_id = run_id  # Return run_id so user can reference it
  ))
}


################################################################################
# Helper function to check checkpoint status
################################################################################

#' Check Pipeline Progress
#'
#' @description Shows what checkpoints exist for a given run
#' @param checkpoint_dir Directory where checkpoints are stored
#' @param run_id The run ID to check (if NULL, shows all runs)
#' @export
check_pipeline_progress <- function(checkpoint_dir = NULL, run_id = NULL) {
  
  if (is.null(checkpoint_dir)) {
    checkpoint_dir <- file.path(main_dir, "checkpoints")
  }
  
  if (!dir.exists(checkpoint_dir)) {
    cat("❌ No checkpoint directory found at:", checkpoint_dir, "\n")
    return(invisible(NULL))
  }
  
  # Find all checkpoint files
  all_files <- list.files(checkpoint_dir, pattern = "\\.rds$", full.names = FALSE)
  
  if (length(all_files) == 0) {
    cat("❌ No checkpoints found\n")
    return(invisible(NULL))
  }
  
  # Extract run IDs
  run_ids <- unique(gsub("_step.*", "", all_files))
  
  cat("════════════════════════════════════════════════════════════\n")
  cat("CHECKPOINT STATUS\n")
  cat("════════════════════════════════════════════════════════════\n\n")
  
  # If specific run_id requested
  if (!is.null(run_id)) {
    run_ids <- run_ids[run_ids == run_id]
  }
  
  step_names <- c(
    "step01" = "Compatibility check",
    "step03" = "Relaxed analysis",
    "step04" = "Master relaxed",
    "step05" = "Stringent analysis",
    "step06" = "Master stringent",
    "step07" = "Consensus",
    "step08" = "Summaries"
  )
  
  for (rid in run_ids) {
    cat(sprintf("🔖 Run ID: %s\n", rid))
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    run_files <- all_files[grepl(paste0("^", rid), all_files)]
    
    for (step in names(step_names)) {
      step_file <- run_files[grepl(step, run_files)]
      if (length(step_file) > 0) {
        file_path <- file.path(checkpoint_dir, step_file)
        file_time <- file.info(file_path)$mtime
        cat(sprintf("  ✅ %s: %s (%s)\n", 
                    step, step_names[step], 
                    format(file_time, "%Y-%m-%d %H:%M")))
      } else {
        cat(sprintf("  ⬜ %s: %s\n", step, step_names[step]))
      }
    }
    
    cat("\n")
  }
  
  cat("To resume a run:\n")
  cat("  results <- run_complete_sip_analysis(\n")
  cat("    studies = studies,\n")
  cat("    resume = TRUE,\n")
  cat("    run_id = \"<run_id_here>\"\n")
  cat("  )\n")
  
  invisible(run_ids)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Check Pipeline Compatibility
#'
#' @description [Function description]
#' @export
check_pipeline_compatibility <- function(studies) {
  
  cat("Checking pipeline compatibility...\n")
  
  compatible <- 0
  issues <- list()
  
  for (study_name in names(studies)) {
    physeq <- studies[[study_name]]
    study_issues <- character()
    
    # Check basic requirements
    sample_info <- tryCatch({
      identify_labeled_unlabeled_samples(physeq)
    }, error = function(e) NULL)
    
    if (is.null(sample_info)) {
      study_issues <- c(study_issues, "sample identification failed")
    } else {
      if (sample_info$n_labeled < 2) {
        study_issues <- c(study_issues, "insufficient labeled samples")
      }
      if (sample_info$n_unlabeled < 2) {
        study_issues <- c(study_issues, "insufficient unlabeled samples")
      }
    }
    
    # Check gradient data
    study_type <- tryCatch({
      detect_study_type(physeq, sample_info)
    }, error = function(e) NULL)
    
    if (is.null(study_type)) {
      study_issues <- c(study_issues, "no valid gradient data")
    }
    
    # Check data quality
    if (ntaxa(physeq) < 10) {
      study_issues <- c(study_issues, "too few taxa")
    }
    
    lib_sizes <- sample_sums(physeq)
    if (median(lib_sizes) < 1000) {
      study_issues <- c(study_issues, "low sequencing depth")
    }
    
    # Record results
    if (length(study_issues) == 0 ||
        !any(grepl("no valid gradient|sample identification failed", study_issues))) {
      compatible <- compatible + 1
      cat(sprintf("  ✓ %s: compatible\n", study_name))
    } else {
      issues[[study_name]] <- study_issues
      cat(sprintf("  ✗ %s: %s\n", study_name, paste(study_issues, collapse = ", ")))
    }
  }
  
  overall_rate <- compatible / length(studies)
  
  cat(sprintf("\nCompatibility: %d/%d studies (%.0f%%)\n",
              compatible, length(studies), overall_rate * 100))
  
  return(list(
    compatible = compatible,
    total = length(studies),
    overall_rate = overall_rate,
    issues = issues
  ))
}

#' Get Study Overview
#'
#' @description [Function description]
#' @export
get_study_overview <- function(physeq, study_name = "study", silent = FALSE) {
  
  sample_info <- tryCatch({
    identify_labeled_unlabeled_samples(physeq)
  }, error = function(e) NULL)
  
  if (is.null(sample_info)) {
    return(NULL)
  }
  
  study_type <- tryCatch({
    detect_study_type(physeq, sample_info)
  }, error = function(e) NULL)
  
  overview <- list(
    study_name = study_name,
    n_samples = nsamples(physeq),
    n_taxa = ntaxa(physeq),
    n_labeled = sample_info$n_labeled,
    n_unlabeled = sample_info$n_unlabeled,
    study_type = if(!is.null(study_type)) study_type$type else "unknown",
    n_substrates = length(sample_info$substrates),
    median_lib_size = median(sample_sums(physeq))
  )
  
  if (!silent) {
    cat(sprintf("\n=== %s ===\n", study_name))
    cat(sprintf("Samples: %d (%d labeled, %d unlabeled)\n",
                overview$n_samples, overview$n_labeled, overview$n_unlabeled))
    cat(sprintf("Taxa: %d\n", overview$n_taxa))
    cat(sprintf("Study type: %s\n", overview$study_type))
    if (overview$n_substrates > 0) {
      cat(sprintf("Substrates: %d\n", overview$n_substrates))
    }
    cat(sprintf("Median library size: %d\n", overview$median_lib_size))
  }
  
  return(overview)
}

#' Save Results
#'
#' @description [Function description]
#' @export
save_results <- function(results_relaxed, master_relaxed,
                         results_stringent = NULL, master_stringent = NULL,
                         consensus_results = NULL, summaries = NULL,
                         output_prefix, output_dir) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saved_files <- list()
  
  # Save relaxed results
  relaxed_file <- file.path(output_dir, "raw_data",
                            paste0(output_prefix, "_relaxed_", timestamp, ".csv"))
  ensure_dir_for_file(relaxed_file)
  write_csv(master_relaxed, relaxed_file)
  saved_files$relaxed_csv <- relaxed_file
  
  # Save stringent results
  if (!is.null(master_stringent)) {
    stringent_file <- file.path(output_dir, "raw_data",
                                paste0(output_prefix, "_stringent_", timestamp, ".csv"))
    ensure_dir_for_file(stringent_file)
    write_csv(master_stringent, stringent_file)
    saved_files$stringent_csv <- stringent_file
  }
  
  # Save consensus
  if (!is.null(consensus_results) && !is.null(consensus_results$incorporators)) {
    consensus_file <- file.path(output_dir, "consensus_results",
                                paste0(output_prefix, "_consensus_", timestamp, ".csv"))
    ensure_dir_for_file(consensus_file)
    write_csv(consensus_results$incorporators, consensus_file)
    saved_files$consensus_csv <- consensus_file
  }
  
  # Save Excel report
  if (!is.null(summaries)) {
    excel_file <- save_to_excel(master_relaxed, summaries, output_dir)
    saved_files$excel_report <- excel_file
    
    text_file <- generate_text_summary(master_relaxed, summaries, output_dir)
    saved_files$text_summary <- text_file
  }
  
  # Save RData
  rdata_file <- file.path(output_dir, "raw_data",
                          paste0(output_prefix, "_complete_", timestamp, ".RData"))
  ensure_dir_for_file(rdata_file)
  save(results_relaxed, master_relaxed,
       results_stringent, master_stringent,
       consensus_results, summaries,
       file = rdata_file)
  saved_files$rdata <- rdata_file
  
  return(saved_files)
}

# ============================================================================
# MEMORY MANAGEMENT
# ============================================================================

#' Report Memory Usage
#'
#' @description [Function description]
#' @export
report_memory_usage <- function() {
  cat("\n📊 MEMORY USAGE REPORT\n")
  cat("======================\n")
  
  # Get all objects
  obj_sizes <- sapply(ls(envir = .GlobalEnv), function(x) {
    tryCatch({
      object.size(get(x, envir = .GlobalEnv))
    }, error = function(e) 0)
  })
  
  # Sort by size
  obj_sizes <- sort(obj_sizes[obj_sizes > 0], decreasing = TRUE)
  
  if (length(obj_sizes) > 0) {
    # Convert to MB
    obj_mb <- obj_sizes / 1024^2
    
    # Top 10
    cat("\nTop objects by size:\n")
    for (i in 1:min(10, length(obj_mb))) {
      cat(sprintf("  %2d. %-30s: %8.2f MB\n", i, names(obj_mb)[i], obj_mb[i]))
    }
    
    cat(sprintf("\nTotal memory used: %.2f MB\n", sum(obj_mb)))
  }
  
  # Cache sizes
  tax_cache_size <- tryCatch({
    object.size(.GlobalEnv$.sip_taxonomy_cache) / 1024^2
  }, error = function(e) 0)
  
  if (tax_cache_size > 0) {
    cat(sprintf("Taxonomy cache: %.2f MB\n", tax_cache_size))
  }
  
  invisible(NULL)
}

# ============================================================================
# INITIALIZATION MESSAGE
# ============================================================================

cat("\n")
cat("████████████████████████████████████████████████████████████\n")
cat("█                                                          █\n")
cat("█     SIP META-ANALYSIS PIPELINE v1.0 - READY             █\n")
cat("█                                                          █\n")
cat("████████████████████████████████████████████████████████████\n")
cat("\n")


cat("QUICK START:\n")
cat("============\n")
cat("# Run complete analysis:\n")
cat("results <- run_complete_sip_analysis(\n")
cat("  studies = studies,\n")
cat("  heavy_threshold = 1.700,       # Relaxed threshold\n")
cat("  run_dual_threshold = TRUE,     # Also run stringent\n")
cat("  stringent_threshold = 1.725,   # Stringent threshold\n")
cat("  n_cores = 4                    # Parallel cores\n")
cat(")\n\n")

cat("# Access results:\n")
cat("master_table <- results$master_relaxed\n")
cat("consensus <- results$consensus$incorporators\n")
cat("summaries <- results$summaries\n\n")

cat("BIOLOGICAL PARAMETERS:\n")
cat("=====================\n")
cat(sprintf("Heavy threshold: %.3f g/ml\n", BIOLOGICAL_PARAMS$density$heavy_threshold))
cat(sprintf("Light threshold: %.3f g/ml\n", BIOLOGICAL_PARAMS$density$light_threshold))
cat(sprintf("Min samples per group: %d\n", BIOLOGICAL_PARAMS$statistics$min_samples_per_group))
cat(sprintf("Log2FC threshold: %.1f\n", BIOLOGICAL_PARAMS$statistics$log2fc_threshold))
cat(sprintf("P-value threshold: %.2f\n", BIOLOGICAL_PARAMS$statistics$padj_threshold))
cat("\n")

cat("📁 Output directory:", main_dir, "\n")
cat("📊 For memory usage: report_memory_usage()\n")
cat("🗑️  To clear cache: clear_cache()\n\n")

cat("Ready to analyze! 🚀\n\n")
