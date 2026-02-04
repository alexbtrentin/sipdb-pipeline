################################################################################
#' @title SIP Database Builder v1.0 - COMPREHENSIVE DOCUMENTED VERSION
#' 
#' @description
#'   Complete database builder for SIP Meta-Analysis Pipeline results.
#'   Creates SQLite databases, FASTA files, and BLAST databases from
#'   pipeline outputs.
#' 
#' @version 1.0
#' @date 2025
#' @author SIP Database Team
#' 
#' @details
#'   This builder creates a comprehensive database from SIP pipeline outputs:
#'   
#'   **Features:**
#'   - Normalized database structure with relationships
#'   - Support for dual thresholds (relaxed + stringent)
#'   - Consensus incorporator integration
#'   - Automatic FASTA file generation
#'   - BLAST database creation
#'   - Optimized indexes for fast queries
#'   - Multiple pre-built views
#'   
#'   **Database Schema:**
#'   - studies: Study-level metadata
#'   - taxonomy: Taxonomic classifications
#'   - asv_sequences: DNA sequences with metrics
#'   - sip_results_relaxed: Results at relaxed threshold
#'   - sip_results_stringent: Results at stringent threshold
#'   - consensus_incorporators: High-confidence incorporators
#'   
#'   **Output Files:**
#'   - SQLite database (.sqlite)
#'   - FASTA file (.fasta)
#'   - BLAST database files
#'   - Build report (.txt)
#' 
#' @note Compatible with SIP Pipeline v1.0+
################################################################################

################################################################################
# Unified Database Builder v1.0 - FIXED FOR STRINGENT VIEW
# Pure Base R - Studies deduplication fixed + Stringent view added
################################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(dplyr)
  library(DBI)
  library(RSQLite)
  library(readr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

options(stringsAsFactors = FALSE)

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================


################################################################################
# MAIN DATABASE BUILDER FUNCTION
# Creates comprehensive SQLite database with optional FASTA and BLAST outputs
################################################################################

#' Build Optimized SIP Database
#' 
#' Main function to build comprehensive SQLite database from SIP pipeline results
#' 
#' @param studies [Description]
#' @param pipeline_dir [Description]
#' @param create_blast [Description]
#' @param optimize_structure [Description]
#' @param only_master_asvs [Description]
#' @param include_consensus [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Creates a normalized database structure with:
#'     - Separate tables for results, taxonomy, sequences, and studies
#'     - Optimized indexes for fast queries
#'     - Multiple views for common queries
#'     - Optional FASTA and BLAST database creation
#'     
#'     Supports:
#'     - Relaxed threshold results (required)
#'     - Stringent threshold results (optional)
#'     - Consensus incorporators (optional)
#' 
#' @export
build_optimized_sip_database_v7 <- function(
    studies,
    pipeline_dir = "~/Documents/sipnew/files/substrate_specific_pipeline_v2.0",
    create_blast = TRUE,
    optimize_structure = TRUE,
    only_master_asvs = TRUE,
    include_consensus = TRUE
) {

  cat("\n=====================================\n")
  cat("🚀 UNIFIED DATABASE BUILDER v1.0\n")
  cat("   For Pipeline v1.0 with Dual Thresholds\n")
  cat("   FIXED: Stringent view now created\n")
  cat("=====================================\n\n")

  start_time <- Sys.time()

  raw_data_dir <- file.path(pipeline_dir, "raw_data")
  consensus_dir <- file.path(pipeline_dir, "consensus_results")
  output_dir <- file.path(pipeline_dir, "database")
  sequence_db_dir <- file.path(pipeline_dir, "sequence_database_comprehensive")

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  # Create output directory
  dir.create(sequence_db_dir, showWarnings = FALSE, recursive = TRUE)

  db_path <- file.path(output_dir, "sip_integrated_v1.sqlite")
  fasta_path <- file.path(sequence_db_dir, "all_sequences_v1.fasta")
  blast_db_path <- file.path(sequence_db_dir, "sip_all_sequences_db_v1")
  report_path <- file.path(output_dir, "database_build_report_v1.txt")

  tryCatch({

    # ========== STEP 1: READ MASTER TABLES ==========
    cat("📁 STEP 1: Loading master tables\n")
    cat("---------------------------------\n")

    relaxed_file <- find_latest_master_table(raw_data_dir, pattern = "relaxed")
    if (is.null(relaxed_file)) {
      relaxed_file <- find_latest_master_table(raw_data_dir, pattern = "master_table")
    }

    master_data <- read_master_table_safe(relaxed_file) %>%
      clean_and_prepare_master_data_v7()

    cat("✅ Loaded relaxed results: ", nrow(master_data), " records\n")
  # Get unique values
    cat("   Unique ASVs: ", n_distinct(master_data$ASV), "\n")

    stringent_data <- NULL
    stringent_file <- find_latest_master_table(raw_data_dir, pattern = "stringent",
                                               must_exist = FALSE)
    if (!is.null(stringent_file)) {
      stringent_data <- read_master_table_safe(stringent_file) %>%
        clean_and_prepare_master_data_v7() %>%
  # Transform data columns
        mutate(threshold_type = "stringent")
      cat("✅ Loaded stringent results: ", nrow(stringent_data), " records\n")
    }

    consensus_data <- NULL
    if (include_consensus && dir.exists(consensus_dir)) {
      consensus_files <- list.files(consensus_dir, pattern = "_consensus_.*\\.csv$",
                                    full.names = TRUE)
      if (length(consensus_files) > 0) {
        consensus_file <- consensus_files[order(file.mtime(consensus_files),
                                                decreasing = TRUE)][1]
        consensus_data <- read_master_table_safe(consensus_file)
        cat("✅ Loaded consensus incorporators: ", nrow(consensus_data), " taxa\n")
      }
    }

    cat("\n")

    # ========== STEP 2: EXTRACT SEQUENCES ==========
    cat("🧬 STEP 2: Extracting sequences\n")
    cat("---------------------------------\n")

    sequences <- extract_master_sequences_only(studies, master_data)
    cat("✅ Extracted ", length(sequences$sequences), " sequences\n\n")

    # ========== STEP 3: BUILD DATABASE ==========
    cat("💾 STEP 3: Building database v1.0\n")
    cat("---------------------------------\n")

  # Remove existing file
    if (file.exists(db_path)) file.remove(db_path)
  # Establish database connection
    con <- dbConnect(SQLite(), db_path)

    if (optimize_structure) {
      build_normalized_database_v7(con, master_data, stringent_data,
                                   consensus_data, sequences)
    } else {
      build_simple_database_v7(con, master_data, sequences)
    }

    # ========== STEP 4: CREATE FASTA & BLAST ==========
    if (length(sequences$sequences) > 0) {
      cat("\n📝 STEP 4: Writing FASTA file\n")
      cat("---------------------------------\n")
      write_fasta_file(sequences$sequences, fasta_path)

      if (create_blast && check_blast_available()) {
        cat("\n🔥 STEP 5: Creating BLAST database\n")
        cat("---------------------------------\n")
        create_blast_database(fasta_path, blast_db_path)
      }
    }

    # ========== WRITE REPORT ==========
    write_build_report_v7(master_data, stringent_data, consensus_data,
                          sequences, report_path)

    # ========== SUMMARY ==========
    end_time <- Sys.time()
    runtime <- difftime(end_time, start_time, units = "secs")

    stats <- get_database_stats_v7(con, optimize_structure)

  # Close database connection
    dbDisconnect(con)

    cat("\n=====================================\n")
    cat("🎉 DATABASE BUILD COMPLETE!\n")
    cat("=====================================\n")
    cat("Runtime: ", round(runtime, 2), " seconds\n")
    cat("Database: ", db_path, "\n")
    cat("  - Sequences: ", stats$n_sequences, "\n")
    cat("  - Relaxed results: ", stats$n_relaxed, "\n")
    if (!is.null(stringent_data)) {
      cat("  - Stringent results: ", stats$n_stringent, "\n")
    }
    if (!is.null(consensus_data)) {
      cat("  - Consensus incorporators: ", stats$n_consensus, "\n")
    }
    cat("  - Studies: ", stats$n_studies, "\n")
    cat("  - Size: ", round(file.size(db_path) / 1024^2, 1), " MB\n")

  # Get unique values
    coverage <- round(100 * stats$n_sequences / n_distinct(master_data$ASV), 1)
    cat("  - Coverage: ", coverage, "% of ASVs have sequences\n")

    if (file.exists(fasta_path)) {
      cat("FASTA: ", fasta_path, "\n")
    }
    if (file.exists(paste0(blast_db_path, ".nhr"))) {
      cat("BLAST DB: ", blast_db_path, "\n")
    }
    cat("Report: ", report_path, "\n")

    cat("\n✨ Database optimized for pipeline v1.0!\n")

    return(list(
      success = TRUE,
      db_path = db_path,
      fasta_path = fasta_path,
      blast_db = blast_db_path,
      n_sequences = stats$n_sequences,
      n_relaxed = stats$n_relaxed,
      n_stringent = if(!is.null(stringent_data)) stats$n_stringent else 0,
      n_consensus = if(!is.null(consensus_data)) stats$n_consensus else 0,
      coverage = coverage,
      runtime = runtime
    ))

  }, error = function(e) {
    cat("\n❌ ERROR: ", e$message, "\n")
    cat("Traceback:\n")
    print(traceback())
  # Close database connection
    if (exists("con")) try(dbDisconnect(con), silent = TRUE)
    return(list(success = FALSE, error = e$message))
  })
}

# ==============================================================================
# FILE READING FUNCTIONS
# ==============================================================================

#' Find Latest Master Table
#' 
#' Locates the most recent master results file based on timestamp
#' 
#' @param raw_dir [Description]
#' @param pattern [Description]
#' @param must_exist [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Searches for files matching pattern and returns the newest by modification time
#' 
#' @export
find_latest_master_table <- function(raw_dir, pattern = "master_table",
                                     must_exist = TRUE) {
  if (!dir.exists(raw_dir)) {
    if (must_exist) stop("raw_data dir not found: ", raw_dir)
    return(NULL)
  }

  search_patterns <- c(
    paste0(".*", pattern, ".*_\\d{8}_\\d{6}\\.csv$"),
    paste0(".*", pattern, ".*\\.csv$")
  )

  cands <- character()
  for (pat in search_patterns) {
    found_files <- list.files(raw_dir, pattern = pat, full.names = TRUE)
    cands <- c(cands, found_files)
  }

  if (!length(cands)) {
    if (must_exist) {
      cat("  Searched for patterns:\n")
      for (pat in search_patterns) {
        cat("    ", pat, "\n")
      }
      stop("No master tables found in: ", raw_dir)
    }
    return(NULL)
  }

  info <- file.info(cands)
  latest <- cands[order(info$mtime, decreasing = TRUE)][1]
  cat("  Using: ", basename(latest), "\n")
  latest
}

#' Safely Read Master Table
#' 
#' Reads CSV file with error handling and fallback strategies
#' 
#' @param filepath [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Tries standard read first, falls back to all-character columns if needed
#' 
#' @export
read_master_table_safe <- function(filepath) {
  cat("  Reading: ", basename(filepath), "\n")

  data <- tryCatch({
  # Read CSV file
    read_csv(filepath, show_col_types = FALSE,
             col_types = cols(.default = col_guess()),
             guess_max = 10000)
  }, error = function(e) {
    cat("    Standard read failed, trying with all character columns...\n")
  # Read CSV file
    read_csv(filepath, show_col_types = FALSE,
             col_types = cols(.default = col_character()))
  })

  probs <- problems(data)
  if (nrow(probs) > 0) {
    cat("    ⚠ ", nrow(probs), " parsing issues\n")
  }

  numeric_cols <- c("log2FC", "pvalue", "padj", "baseMean", "window_center",
                    "density_threshold", "gradient_position", "n_labeled",
                    "n_unlabeled", "window_width", "k_per_group",
                    "confidence_score", "mean_log2FC", "min_padj",
                    "window_weight", "seq_length", "gc_content")

  for (col in numeric_cols) {
    if (col %in% names(data) && !is.numeric(data[[col]])) {
      data[[col]] <- suppressWarnings(as.numeric(data[[col]]))
    }
  }

  cat("    ✓ Loaded ", nrow(data), " rows, ", ncol(data), " columns\n")

  return(data)
}

#' Clean and Prepare Master Data
#' 
#' Standardizes column names and data types for database insertion
#' 
#' @param df [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Handles missing columns, renames for consistency, adds derived fields
#' 
#' @export
clean_and_prepare_master_data_v7 <- function(df) {
  if ("taxa" %in% names(df) && !"ASV" %in% names(df)) {
    df$ASV <- df$taxa
  }

  tax_cols <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (col in tax_cols) {
    if (col %in% names(df)) {
      df[[col]] <- clean_taxonomy_string(df[[col]])
    }
  }

  if ("Genus" %in% names(df)) {
    df <- df %>%
  # Filter data rows
      filter(is.na(Genus) |
               (!grepl("Chloroplast|Mitochondria", Genus, ignore.case = TRUE)))
  }

  df <- df %>%
  # Filter data rows
    filter(
      !is.na(ASV),
      !is.na(log2FC),
      !is.na(pvalue),
      pvalue >= 0 & pvalue <= 1,
      abs(log2FC) <= 50
    )

  if ("Domain" %in% names(df)) {
    df <- df %>%
      filter(is.na(Domain) | 
               !Domain %in% c("Archaea", "Eukarya", "Eukaryota"))
  }

  df <- df %>%
  # Transform data columns
    mutate(
      is_enriched = (pvalue < 0.05) & (abs(log2FC) > 1.0),
      threshold_type = if("threshold_type" %in% names(.)) threshold_type else "relaxed"
    )

  df
}

#' Clean Taxonomy String
#' 
#' @description [Add description]
#' @export
clean_taxonomy_string <- function(x) {
  x %>%
    str_remove("^[a-z]__") %>%
    str_remove_all("\\[.*?\\]") %>%
    str_trim() %>%
    na_if("") %>%
    na_if("Unknown") %>%
    na_if("Unclassified") %>%
    na_if("uncultured") %>%
    na_if("unidentified")
}

# ==============================================================================
# SEQUENCE EXTRACTION
# ==============================================================================

#' Extract ASV Sequences
#' 
#' Extracts DNA sequences for ASVs from phyloseq objects
#' 
#' @param studies [Description]
#' @param master_data [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Searches through study list to find sequences, tracks missing ASVs
#' 
#' @export
extract_master_sequences_only <- function(studies, master_data) {
  master_asvs <- unique(master_data$ASV[!is.na(master_data$ASV)])
  cat("  Target ASVs: ", length(master_asvs), "\n")

  all_sequences <- DNAStringSet()
  source_info <- list()
  sequences_found <- character()

  for (study_id in names(studies)) {
    cat("  Processing ", study_id, "...")
    physeq <- studies[[study_id]]

    seqs <- tryCatch(refseq(physeq), error = function(e) NULL)
    if (is.null(seqs)) {
      cat(" no refseq\n")
      next
    }

    if (!inherits(seqs, "DNAStringSet")) {
      seqs <- DNAStringSet(seqs)
    }

    seq_names <- names(seqs)
    matching_asvs <- intersect(master_asvs, seq_names)

    if (length(matching_asvs) > 0) {
      new_asvs <- setdiff(matching_asvs, sequences_found)
      if (length(new_asvs) > 0) {
        all_sequences <- c(all_sequences, seqs[new_asvs])
        sequences_found <- c(sequences_found, new_asvs)
      }

      source_info[[study_id]] <- data.frame(
        asv_id = matching_asvs,
        study_id = study_id,
        stringsAsFactors = FALSE
      )

      cat(" found ", length(matching_asvs), "\n")
    } else {
      cat(" 0 found\n")
    }
  }

  missing_asvs <- setdiff(master_asvs, sequences_found)
  cat("\n  Summary:\n")
  cat("    - With sequences: ", length(sequences_found), "\n")
  cat("    - Without sequences: ", length(missing_asvs), "\n")
  cat("    - Coverage: ", round(100 * length(sequences_found) / length(master_asvs), 1), "%\n")

  list(
    sequences = all_sequences,
    sources = bind_rows(source_info),
    missing_asvs = missing_asvs
  )
}

# ==============================================================================
# NORMALIZED DATABASE
# ==============================================================================


################################################################################
# DATABASE CONSTRUCTION FUNCTIONS
# Build optimized database structure with tables and indexes
################################################################################

#' Build Normalized Database Structure
#' 
#' Creates optimized database with separate tables and relationships
#' 
#' @param con [Description]
#' @param master_data [Description]
#' @param stringent_data [Description]
#' @param consensus_data [Description]
#' @param sequences [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Database schema:
#'     - studies: Study metadata
#'     - taxonomy: Taxonomic classifications
#'     - asv_sequences: DNA sequences with GC content
#'     - sip_results_relaxed: Results at relaxed threshold
#'     - sip_results_stringent: Results at stringent threshold (optional)
#'     - consensus_incorporators: High-confidence taxa (optional)
#'     Plus multiple views for common queries
#' 
#' @export
build_normalized_database_v7 <- function(con, master_data, stringent_data,
                                         consensus_data, sequences) {
  cat("  Building normalized structure v1.0...\n")

  dbBegin(con)

  # ========== 1. STUDIES TABLE - FIXED: ONE ROW PER STUDY ==========
  cat("  Creating studies table...\n")

  # Ensure taxa -> ASV
  if ("taxa" %in% names(master_data) && !"ASV" %in% names(master_data)) {
    master_data$ASV <- master_data$taxa
  }

  study_cols <- c("study_id", "Bioproject_ID", "DOI_URL",
                  "env_biome", "env_feature", "env_material",
                  "environment_label")

  available_study_cols <- intersect(study_cols, names(master_data))

  if (length(available_study_cols) == 0) {
    studies <- unique(data.frame(
      study_id = master_data$study_id,
      stringsAsFactors = FALSE
    ))
  } else {
    # Extract columns with base R
    studies <- master_data[, available_study_cols, drop = FALSE]

    # Keep only FIRST occurrence of each study_id
    studies <- studies[!duplicated(studies$study_id), , drop = FALSE]
  }

  # Add sequential ID
  studies$id <- seq_len(nrow(studies))

  # Rename columns with base R
  colnames_studies <- names(studies)
  colnames_studies[colnames_studies == "Bioproject_ID"] <- "Bioproject"
  colnames_studies[colnames_studies == "DOI_URL"] <- "DOI"
  names(studies) <- colnames_studies

  # Reset row names
  rownames(studies) <- NULL

  # Write data to table
  dbWriteTable(con, "studies", studies, overwrite = TRUE)
  cat("    ✓ ", nrow(studies), " unique studies\n")

  # ========== 2. TAXONOMY TABLE - PURE BASE R ==========
  cat("  Creating taxonomy table...\n")

  tax_cols <- c("ASV", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  available_tax_cols <- intersect(tax_cols, names(master_data))

  if (length(available_tax_cols) == 0 || !"ASV" %in% available_tax_cols) {
    stop("ERROR: No ASV column found")
  }

  # Extract with base R
  taxonomy <- master_data[, available_tax_cols, drop = FALSE]
  # Remove duplicates
  taxonomy <- taxonomy[!duplicated(taxonomy), , drop = FALSE]
  # Add tax_id
  taxonomy$tax_id <- seq_len(nrow(taxonomy))

  # Rename Order -> TaxOrder
  colnames_tax <- names(taxonomy)
  colnames_tax[colnames_tax == "Order"] <- "TaxOrder"
  names(taxonomy) <- colnames_tax

  rownames(taxonomy) <- NULL

  # Write data to table
  dbWriteTable(con, "taxonomy", taxonomy, overwrite = TRUE)
  cat("    ✓ ", nrow(taxonomy), " unique taxa\n")

  # ========== 3. SEQUENCES TABLE ==========
  if (length(sequences$sequences) > 0) {
    cat("  Creating sequences table...\n")

    seq_df <- data.frame(
      asv_id = names(sequences$sequences),
      sequence = as.character(sequences$sequences),
      seq_length = width(sequences$sequences),
      stringsAsFactors = FALSE
    )

    seq_df$gc_content <- round(sapply(seq_df$sequence, function(s) {
      chars <- strsplit(toupper(s), "")[[1]]
      gc_count <- sum(chars %in% c("G", "C"))
      if (length(chars) > 0) gc_count / length(chars) * 100 else 0
    }), 2)

  # Write data to table
    dbWriteTable(con, "asv_sequences", seq_df, overwrite = TRUE)
    cat("    ✓ ", nrow(seq_df), " sequences\n")
  } else {
    seq_df <- data.frame(asv_id = character(), sequence = character(),
                         seq_length = integer(), gc_content = numeric(),
                         stringsAsFactors = FALSE)
  # Write data to table
    dbWriteTable(con, "asv_sequences", seq_df, overwrite = TRUE)
    cat("    ⚠ No sequences\n")
  }

  # ========== 4. RELAXED RESULTS - PURE BASE R ==========
  cat("  Creating relaxed results table...\n")

  results_cols <- c("study_pk", "tax_id", "ASV", "method", "log2FC",
                    "pvalue", "padj", "isotopolog", "is_enriched",
                    "analysis_type", "substrate_analyzed",
                    "gradient_position", "window_name", "window_mode",
                    "k_per_group", "window_center", "density_threshold",
                    "comparison_type", "baseMean", "window_weight")

  # Create lookup for study_pk
  study_lookup <- setNames(studies$id, studies$study_id)
  master_data$study_pk <- study_lookup[master_data$study_id]

  # Create lookup for tax_id
  tax_lookup <- setNames(taxonomy$tax_id, taxonomy$ASV)
  master_data$tax_id <- tax_lookup[master_data$ASV]

  # Select available columns
  available_results_cols <- intersect(results_cols, names(master_data))
  results_relaxed <- master_data[, available_results_cols, drop = FALSE]

  # Add metadata
  results_relaxed$result_id <- seq_len(nrow(results_relaxed))
  results_relaxed$threshold_type <- "relaxed"

  rownames(results_relaxed) <- NULL

  # Write data to table
  dbWriteTable(con, "sip_results_relaxed", results_relaxed, overwrite = TRUE)
  cat("    ✓ ", nrow(results_relaxed), " results\n")

  # ========== 5. STRINGENT RESULTS - PURE BASE R ==========
  if (!is.null(stringent_data)) {
    cat("  Creating stringent results table...\n")

    # Ensure ASV exists
    if ("taxa" %in% names(stringent_data) && !"ASV" %in% names(stringent_data)) {
      stringent_data$ASV <- stringent_data$taxa
    }

    # Add lookups
    stringent_data$study_pk <- study_lookup[stringent_data$study_id]
    stringent_data$tax_id <- tax_lookup[stringent_data$ASV]

    # Select columns
    available_stringent_cols <- intersect(results_cols, names(stringent_data))
    results_stringent <- stringent_data[, available_stringent_cols, drop = FALSE]

    # Add metadata
    results_stringent$result_id <- seq_len(nrow(results_stringent))
    results_stringent$threshold_type <- "stringent"

    rownames(results_stringent) <- NULL

  # Write data to table
    dbWriteTable(con, "sip_results_stringent", results_stringent, overwrite = TRUE)
    cat("    ✓ ", nrow(results_stringent), " results\n")
  }

  # ========== 6. CONSENSUS - PURE BASE R ==========
  if (!is.null(consensus_data)) {
    cat("  Creating consensus table...\n")

    if ("taxa" %in% names(consensus_data) && !"ASV" %in% names(consensus_data)) {
      consensus_data$ASV <- consensus_data$taxa
    }

    consensus_cols <- c("study_pk", "tax_id", "ASV", "method",
                        "substrate_analyzed", "mean_log2FC", "min_padj",
                        "confidence", "n_studies", "methods_sig")

    # Add lookups if study_id exists
    if ("study_id" %in% names(consensus_data)) {
      consensus_data$study_pk <- study_lookup[consensus_data$study_id]
    }
    consensus_data$tax_id <- tax_lookup[consensus_data$ASV]

    # Select columns
    available_consensus_cols <- intersect(consensus_cols, names(consensus_data))
    consensus_table <- consensus_data[, available_consensus_cols, drop = FALSE]

    # Add ID
    consensus_table$consensus_id <- seq_len(nrow(consensus_table))

    rownames(consensus_table) <- NULL

  # Write data to table
    dbWriteTable(con, "consensus_incorporators", consensus_table, overwrite = TRUE)
    cat("    ✓ ", nrow(consensus_table), " incorporators\n")
  }

  # ========== 7. INDEXES ==========
  cat("  Creating indexes...\n")
  create_indexes_v7(con, !is.null(stringent_data), !is.null(consensus_data))

  # ========== 8. VIEWS ==========
  cat("  Creating views...\n")
  create_views_v7(con, !is.null(stringent_data), !is.null(consensus_data))

  dbCommit(con)

  cat("  Optimizing...\n")
  # Optimize database size
  dbExecute(con, "VACUUM")
  dbExecute(con, "ANALYZE")

  cat("  ✅ Database created\n")
}

# ==============================================================================
# INDEXES
# ==============================================================================

#' Create Indexes V7
#' 
#' @description [Add description]
#' @export
create_indexes_v7 <- function(con, has_stringent, has_consensus) {
  indexes <- c(
    "CREATE INDEX idx_studies_id ON studies(study_id)",
    "CREATE INDEX idx_tax_asv ON taxonomy(ASV)",
    "CREATE INDEX idx_tax_genus ON taxonomy(Genus)",
    "CREATE INDEX idx_tax_phylum ON taxonomy(Phylum)",
    "CREATE INDEX idx_seq_asv ON asv_sequences(asv_id)",
    "CREATE INDEX idx_relax_study ON sip_results_relaxed(study_pk)",
    "CREATE INDEX idx_relax_asv ON sip_results_relaxed(ASV)",
    "CREATE INDEX idx_relax_pval ON sip_results_relaxed(pvalue)",
    "CREATE INDEX idx_relax_enriched ON sip_results_relaxed(is_enriched)"
  )

  if (has_stringent) {
    indexes <- c(indexes,
                 "CREATE INDEX idx_string_asv ON sip_results_stringent(ASV)",
                 "CREATE INDEX idx_string_pval ON sip_results_stringent(pvalue)"
    )
  }

  if (has_consensus) {
    indexes <- c(indexes,
                 "CREATE INDEX idx_cons_asv ON consensus_incorporators(ASV)",
                 "CREATE INDEX idx_cons_conf ON consensus_incorporators(confidence)"
    )
  }

  for (idx in indexes) {
    tryCatch(dbExecute(con, idx), error = function(e) NULL)
  }
}


################################################################################
# DATABASE VIEW CREATION
# Generate SQL views for common query patterns
################################################################################

#' Create Database Views
#' 
#' Generates SQL views for common query patterns
#' 
#' @param con [Description]
#' @param has_stringent [Description]
#' @param has_consensus [Description]
#' 
#' @return [Return value description]
#' 
#' @details
#'   Views include:
#'     - enriched_taxa_relaxed: Significantly enriched at relaxed threshold
#'     - enriched_sequences: Sequences of enriched taxa
#'     - consensus_both_thresholds: Taxa enriched at both thresholds
#'     - high_confidence_incorporators: High-quality consensus taxa
#' 
#' @export
create_views_v7 <- function(con, has_stringent, has_consensus) {

  # RELAXED VIEW
  dbExecute(con, "
    CREATE VIEW IF NOT EXISTS analysis_view_relaxed AS
    SELECT
      r.*,
      s.study_id, s.Bioproject, s.DOI, s.environment_label,
      t.Domain, t.Phylum, t.Class, t.TaxOrder, t.Family, t.Genus, t.Species,
      seq.sequence, seq.seq_length, seq.gc_content
    FROM sip_results_relaxed r
    LEFT JOIN studies s ON r.study_pk = s.id
    LEFT JOIN taxonomy t ON r.tax_id = t.tax_id
    LEFT JOIN asv_sequences seq ON r.ASV = seq.asv_id
  ")

  # ========== CREATE STRINGENT VIEW ==========
  if (has_stringent) {
    cat("  Creating analysis_view_stringent...\n")
    dbExecute(con, "
      CREATE VIEW IF NOT EXISTS analysis_view_stringent AS
      SELECT
        r.*,
        s.study_id, s.Bioproject, s.DOI, s.environment_label,
        t.Domain, t.Phylum, t.Class, t.TaxOrder, t.Family, t.Genus, t.Species,
        seq.sequence, seq.seq_length, seq.gc_content
      FROM sip_results_stringent r
      LEFT JOIN studies s ON r.study_pk = s.id
      LEFT JOIN taxonomy t ON r.tax_id = t.tax_id
      LEFT JOIN asv_sequences seq ON r.ASV = seq.asv_id
    ")
  }

  # ENRICHED SEQUENCES VIEW
  dbExecute(con, "
    CREATE VIEW IF NOT EXISTS enriched_sequences AS
    SELECT DISTINCT
      r.ASV,
      t.Genus, t.Phylum,
      seq.sequence,
      MIN(r.pvalue) as best_pvalue,
      MAX(abs(r.log2FC)) as max_log2fc
    FROM sip_results_relaxed r
    JOIN taxonomy t ON r.tax_id = t.tax_id
    LEFT JOIN asv_sequences seq ON r.ASV = seq.asv_id
    WHERE r.is_enriched = 1
    GROUP BY r.ASV
  ")

  # CONSENSUS BOTH THRESHOLDS VIEW
  if (has_stringent) {
    dbExecute(con, "
      CREATE VIEW IF NOT EXISTS consensus_both_thresholds AS
      SELECT
        rel.ASV,
        t.Genus, t.Phylum,
        AVG(rel.log2FC) as mean_fc_relaxed,
        AVG(str.log2FC) as mean_fc_stringent
      FROM sip_results_relaxed rel
      INNER JOIN sip_results_stringent str ON rel.ASV = str.ASV
      JOIN taxonomy t ON rel.tax_id = t.tax_id
      WHERE rel.is_enriched = 1 AND str.is_enriched = 1
      GROUP BY rel.ASV
    ")
  }

  # HIGH CONFIDENCE INCORPORATORS VIEW
  if (has_consensus) {
    dbExecute(con, "
      CREATE VIEW IF NOT EXISTS high_confidence_incorporators AS
      SELECT
        c.*,
        t.Genus, t.Phylum,
        seq.sequence
      FROM consensus_incorporators c
      JOIN taxonomy t ON c.tax_id = t.tax_id
      LEFT JOIN asv_sequences seq ON c.ASV = seq.asv_id
      WHERE c.confidence = 'high_consensus'
    ")
  }
}

# ==============================================================================
# SIMPLE DATABASE
# ==============================================================================

#' Build Simple Database V7
#' 
#' @description [Add description]
#' @export
build_simple_database_v7 <- function(con, master_data, sequences) {
  cat("  Building simple structure...\n")

  dbBegin(con)

  if (!"ASV" %in% names(master_data) && "taxa" %in% names(master_data)) {
    master_data$ASV <- master_data$taxa
  }

  if ("Order" %in% names(master_data)) {
    colnames_master <- names(master_data)
    colnames_master[colnames_master == "Order"] <- "TaxOrder"
    names(master_data) <- colnames_master
  }

  master_data$id <- seq_len(nrow(master_data))
  # Write data to table
  dbWriteTable(con, "sip_results", master_data, overwrite = TRUE)

  if (length(sequences$sequences) > 0) {
    seq_df <- data.frame(
      asv_id = names(sequences$sequences),
      sequence = as.character(sequences$sequences),
      seq_length = width(sequences$sequences),
      stringsAsFactors = FALSE
    )

    seq_df$gc_content <- round(sapply(seq_df$sequence, function(s) {
      chars <- strsplit(toupper(s), "")[[1]]
      gc_count <- sum(chars %in% c("G", "C"))
      if (length(chars) > 0) gc_count / length(chars) * 100 else 0
    }), 2)

  # Write data to table
    dbWriteTable(con, "asv_sequences", seq_df, overwrite = TRUE)
  }

  indexes <- c(
    "CREATE INDEX idx_asv ON sip_results(ASV)",
    "CREATE INDEX idx_study ON sip_results(study_id)",
    "CREATE INDEX idx_pvalue ON sip_results(pvalue)"
  )

  for (idx in indexes) {
    tryCatch(dbExecute(con, idx), error = function(e) NULL)
  }

  dbCommit(con)
  # Optimize database size
  dbExecute(con, "VACUUM")

  cat("  ✅ Simple database created\n")
}

# ==============================================================================
# UTILITIES
# ==============================================================================

#' Get Database Stats V7
#' 
#' @description [Add description]
#' @export
get_database_stats_v7 <- function(con, optimized) {
  if (optimized) {
    n_sequences <- tryCatch({
      res <- dbGetQuery(con, "SELECT COUNT(DISTINCT asv_id) FROM asv_sequences")
      if (is.data.frame(res)) res[[1]] else 0
    }, error = function(e) 0)

    n_relaxed <- tryCatch({
      res <- dbGetQuery(con, "SELECT COUNT(*) FROM sip_results_relaxed")
      if (is.data.frame(res)) res[[1]] else 0
    }, error = function(e) 0)

    n_stringent <- tryCatch({
      tables <- dbListTables(con)
      if ("sip_results_stringent" %in% tables) {
        res <- dbGetQuery(con, "SELECT COUNT(*) FROM sip_results_stringent")
        if (is.data.frame(res)) res[[1]] else 0
      } else {
        0
      }
    }, error = function(e) 0)

    n_consensus <- tryCatch({
      tables <- dbListTables(con)
      if ("consensus_incorporators" %in% tables) {
        res <- dbGetQuery(con, "SELECT COUNT(*) FROM consensus_incorporators")
        if (is.data.frame(res)) res[[1]] else 0
      } else {
        0
      }
    }, error = function(e) 0)

    n_studies <- tryCatch({
      res <- dbGetQuery(con, "SELECT COUNT(DISTINCT study_id) FROM studies")
      if (is.data.frame(res)) res[[1]] else 0
    }, error = function(e) 0)

    return(list(
      n_sequences = n_sequences,
      n_relaxed = n_relaxed,
      n_stringent = n_stringent,
      n_consensus = n_consensus,
      n_studies = n_studies
    ))

  } else {
    list(
      n_sequences = dbGetQuery(con, "SELECT COUNT(DISTINCT asv_id) FROM asv_sequences")[[1]],
      n_relaxed = dbGetQuery(con, "SELECT COUNT(*) FROM sip_results")[[1]],
      n_stringent = 0,
      n_consensus = 0,
      n_studies = dbGetQuery(con, "SELECT COUNT(DISTINCT study_id) FROM sip_results")[[1]]
    )
  }
}


################################################################################
# UTILITY FUNCTIONS
# Helper functions for FASTA, BLAST, and reporting
################################################################################

#' Check Blast Available
#' 
#' @description [Add description]
#' @export
check_blast_available <- function() {
  tryCatch({
  # Execute BLAST database creation
    system2("makeblastdb", "-version", stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) FALSE)
}


################################################################################
# UTILITY FUNCTIONS
# Helper functions for FASTA, BLAST, and reporting
################################################################################

#' Write Fasta File
#' 
#' @description [Add description]
#' @export
write_fasta_file <- function(sequences, fasta_path) {
  # Write sequences to FASTA
  writeXStringSet(sequences, filepath = fasta_path)
  cat("  Wrote ", length(sequences), " sequences\n")
  TRUE
}

#' Create Blast Database
#' 
#' @description [Add description]
#' @export
create_blast_database <- function(fasta_path, blast_prefix) {
  args <- c("-in", fasta_path, "-dbtype", "nucl", "-out", blast_prefix,
            "-title", "SIP_v1.0", "-parse_seqids")

  # Execute BLAST database creation
  system2("makeblastdb", args = args, stdout = TRUE, stderr = TRUE)

  if (file.exists(paste0(blast_prefix, ".nhr"))) {
    cat("  ✅ BLAST database created\n")
    TRUE
  } else {
    cat("  ❌ BLAST failed\n")
    FALSE
  }
}

#' Write Build Report V7
#' 
#' @description [Add description]
#' @export
write_build_report_v7 <- function(master_data, stringent_data, consensus_data,
                                  sequences, report_path) {
  sink(report_path)

  cat("DATABASE BUILD REPORT v1.0\n")
  cat("==========================\n")
  cat("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  cat("RELAXED RESULTS\n")
  cat("  Records: ", nrow(master_data), "\n")
  # Get unique values
  cat("  Unique ASVs: ", n_distinct(master_data$ASV), "\n")
  # Get unique values
  cat("  Studies: ", n_distinct(master_data$study_id), "\n")
  cat("  Enriched: ", sum(master_data$is_enriched, na.rm = TRUE), "\n\n")

  if (!is.null(stringent_data)) {
    cat("STRINGENT RESULTS\n")
    cat("  Records: ", nrow(stringent_data), "\n")
    cat("  Enriched: ", sum(stringent_data$is_enriched, na.rm = TRUE), "\n\n")
  }

  if (!is.null(consensus_data)) {
    cat("CONSENSUS\n")
    cat("  Total: ", nrow(consensus_data), "\n\n")
  }

  cat("SEQUENCES\n")
  cat("  Found: ", length(sequences$sequences), "\n")
  cat("  Missing: ", length(sequences$missing_asvs), "\n")
  # Get unique values
  cat("  Coverage: ", round(100 * length(sequences$sequences) / n_distinct(master_data$ASV), 1), "%\n")

  sink()
}

# ==============================================================================
# INIT
# ==============================================================================

if (interactive()) {
  cat("\n📦 Database Builder v1.0 loaded!\n")
  cat("   Pipeline v1.0 compatible\n")
  cat("   Stringent view created!\n\n")
  cat('Usage: result <- build_optimized_sip_database_v7(ps)\n\n')
}
