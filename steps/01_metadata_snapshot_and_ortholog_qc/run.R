#!/usr/bin/env Rscript

# =============================================================================
# PESTFLY — Step 1–2: Metadata snapshot, label consistency, and ortholog QC
# =============================================================================
#
# Purpose
# -------
# This script prepares the curated metadata and ortholog-level QC objects used by
# the downstream SNP-based origin-tracing workflow.
#
# It performs three main tasks:
#
#   1. Read, clean, and harmonise metadata from:
#        data/000_input_data/metadata.xlsx
#
#   2. Check label consistency between metadata sample IDs and raw PHYLIP labels:
#        data/000_input_data/phy/OG*.fa
#
#   3. Run ortholog-level QC on FASTA alignments from Step 0:
#        results/00_fasta/OG*.fasta
#
# Metadata harmonisation rules
# ----------------------------
#   - Remove hidden spaces, non-breaking spaces, and zero-width characters.
#
#   - Bdors and Blati are expected OMA/read2tree outgroups. They may occur in
#     alignments but are not expected in metadata.
#
#   - Congo reference samples are assigned to:
#        macroregion_3 = Africa
#        subregion     = C_Africa
#
#   - Reunion / Réunion / Reu* samples are harmonised to:
#        country       = Reunion
#        population    = Reunion
#        macroregion_3 = Asia
#        subregion     = SE_Asia
#
#     This is a genetic/source-lineage grouping for the PESTFLY origin-tracing
#     workflow, consistent with the Deschepper interpretation, not a strict
#     geographic classification.
#
# Outputs
# -------
# Written to:
#   results/01_qc/
#
# Main outputs:
#   metadata_clean.tsv / .rds
#   og_qc.tsv / .rds
#   ogs_pass_qc.txt / .rds
#   ogs_fail_qc.txt / .rds
#   sample_missingness.tsv / .rds
#   reference_nonmissing_by_og.rds
#   qc_params.tsv / .rds
#   step1_2_run_info.rds
#
# Label-consistency outputs:
#   label_consistency_by_og_phy.tsv / .rds
#   label_consistency_sample_summary_phy.tsv / .rds
#   label_consistency_labels_not_in_metadata_phy.tsv
#   label_consistency_labels_not_in_metadata_unexpected_phy.tsv
#   label_consistency_metadata_not_seen_phy.tsv
#   label_consistency_report_phy.txt
#
# Run
# ---
#   Rscript steps/step1_2_metadata_snapshot_orthologqc/run.R
#
# Test
# ----
#   Rscript steps/step1_2_metadata_snapshot_orthologqc/run.R --debug_n 100 --cores 4
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(readxl)
  library(Biostrings)
  library(parallel)
})

# =============================================================================
# Repo-root detection
# =============================================================================

find_repo_root <- function(max_up = 10L) {
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  cand <- wd
  
  for (i in seq_len(max_up + 1L)) {
    if (dir.exists(file.path(cand, "data")) && dir.exists(file.path(cand, "steps"))) {
      return(cand)
    }
    
    parent <- normalizePath(file.path(cand, ".."), winslash = "/", mustWork = FALSE)
    if (identical(parent, cand)) break
    cand <- parent
  }
  
  wd
}

repo_root <- find_repo_root()

is_abs_path <- function(p) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/)", p)
}

resolve_path <- function(p, root = repo_root) {
  if (is.null(p) || !nzchar(p)) return(p)
  if (is_abs_path(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  normalizePath(file.path(root, p), winslash = "/", mustWork = FALSE)
}

# =============================================================================
# CLI
# =============================================================================

option_list <- list(
  make_option(
    "--align_dir",
    type = "character",
    default = "results/00_fasta",
    help = "Directory containing Step 0 FASTA alignments [default %default]"
  ),
  make_option(
    "--file_glob",
    type = "character",
    default = "OG*.fasta",
    help = "FASTA file pattern [default %default]"
  ),
  make_option(
    "--metadata_xlsx",
    type = "character",
    default = "data/000_input_data/metadata.xlsx",
    help = "Metadata Excel file [default %default]"
  ),
  make_option(
    "--metadata_sheet",
    type = "character",
    default = "Selection",
    help = "Metadata Excel sheet name [default %default]"
  ),
  make_option(
    "--phy_dir",
    type = "character",
    default = "data/000_input_data/phy",
    help = "Raw PHYLIP directory used for label consistency checks [default %default]"
  ),
  make_option(
    "--phy_glob",
    type = "character",
    default = "OG*.fa",
    help = "Raw PHYLIP file pattern for label consistency checks [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/01_qc",
    help = "Output directory [default %default]"
  ),
  
  # QC thresholds
  make_option(
    "--min_len",
    type = "integer",
    default = 300,
    help = "Minimum alignment length to keep OG [default %default]"
  ),
  make_option(
    "--min_nonmissing_frac",
    type = "double",
    default = 0.50,
    help = "Per-sample minimum non-missing fraction used in occupancy checks [default %default]"
  ),
  make_option(
    "--min_pop_occupancy",
    type = "double",
    default = 0.70,
    help = "Population occupancy threshold: fraction of samples in population represented in OG [default %default]"
  ),
  make_option(
    "--min_pop_occupancy_frac_pops",
    type = "double",
    default = 0.80,
    help = "Required fraction of populations passing min_pop_occupancy [default %default]"
  ),
  make_option(
    "--max_mean_ambig",
    type = "double",
    default = 0.05,
    help = "Maximum mean ambiguous fraction in OG [default %default]"
  ),
  
  # Runtime
  make_option(
    "--exclude_regex",
    type = "character",
    default = "^$",
    help = "Regex of taxa to exclude before QC, if needed [default %default]"
  ),
  make_option(
    "--allowed_extra_labels",
    type = "character",
    default = "Bdors,Blati",
    help = "Comma-separated labels allowed in PHYLIP/FASTA but absent from metadata [default %default]"
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 0,
    help = "Parallel workers. 0 = all detected cores minus one [default %default]"
  ),
  make_option(
    "--debug_n",
    type = "integer",
    default = 0,
    help = "Process only first N FASTA files; 0 = all [default %default]"
  ),
  make_option(
    "--no_sample_missingness",
    action = "store_true",
    default = FALSE,
    help = "Skip writing sample_missingness.tsv/rds [default %default]"
  ),
  
  # Label consistency
  make_option(
    "--label_check",
    type = "character",
    default = "all",
    help = "Label consistency check against raw PHYLIP files: all, first, none [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$align_dir <- resolve_path(opt$align_dir)
opt$metadata_xlsx <- resolve_path(opt$metadata_xlsx)
opt$phy_dir <- resolve_path(opt$phy_dir)
opt$out_dir <- resolve_path(opt$out_dir)

if (opt$cores <= 0) {
  opt$cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
}

allowed_extra_labels <- trimws(unlist(strsplit(opt$allowed_extra_labels, ",")))
allowed_extra_labels <- allowed_extra_labels[nzchar(allowed_extra_labels)]

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Repo root:             ", repo_root)
message("Metadata:              ", opt$metadata_xlsx)
message("Metadata sheet:        ", opt$metadata_sheet)
message("FASTA dir:             ", opt$align_dir)
message("PHYLIP dir:            ", opt$phy_dir)
message("Output dir:            ", opt$out_dir)
message("Cores:                 ", opt$cores)
message("Allowed extra labels:  ", paste(allowed_extra_labels, collapse = ", "))

# =============================================================================
# General helpers
# =============================================================================

standardize_colnames <- function(x) {
  x <- trimws(x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[\u200B-\u200D\uFEFF]", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

as_clean_character <- clean_text

as_clean_logical <- function(x) {
  if (is.logical(x)) return(x)
  
  y <- clean_text(x)
  y <- tolower(y)
  
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE
  
  out
}

# =============================================================================
# Step 1: Read, clean, and harmonise metadata
# =============================================================================

if (!file.exists(opt$metadata_xlsx)) {
  stop("Metadata file not found: ", opt$metadata_xlsx)
}

metadata_raw <- as.data.table(readxl::read_excel(
  path = opt$metadata_xlsx,
  sheet = opt$metadata_sheet
))

setnames(metadata_raw, standardize_colnames(names(metadata_raw)))

required_cols <- c("sample_id", "is_reference", "population", "country")
missing_required <- setdiff(required_cols, names(metadata_raw))

if (length(missing_required) > 0) {
  stop(
    "Metadata is missing required column(s): ",
    paste(missing_required, collapse = ", "),
    "\nFound columns: ",
    paste(names(metadata_raw), collapse = ", ")
  )
}

metadata_clean <- copy(metadata_raw)

for (nm in names(metadata_clean)) {
  if (!inherits(metadata_clean[[nm]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    metadata_clean[[nm]] <- as_clean_character(metadata_clean[[nm]])
  }
}

metadata_clean[, sample_id := as_clean_character(sample_id)]
metadata_clean[, is_reference := as_clean_logical(is_reference)]
metadata_clean[, population := as_clean_character(population)]
metadata_clean[, country := as_clean_character(country)]

if ("macroregion_3" %in% names(metadata_clean)) {
  metadata_clean[, macroregion_3 := as_clean_character(macroregion_3)]
} else {
  metadata_clean[, macroregion_3 := NA_character_]
}

if ("subregion" %in% names(metadata_clean)) {
  metadata_clean[, subregion := as_clean_character(subregion)]
} else {
  metadata_clean[, subregion := NA_character_]
}

# -----------------------------------------------------------------------------
# Reproducible regional harmonisation rules
# -----------------------------------------------------------------------------

# Congo references: Central Africa
metadata_clean[
  (!is.na(country) & country == "Congo") |
    (!is.na(population) & population == "Congo"),
  `:=`(
    macroregion_3 = "Africa",
    subregion = "C_Africa"
  )
]

# Reunion / Réunion / Reu* samples: Asia-linked SE_Asia source lineage.
metadata_clean[
  (!is.na(country) & country %in% c("Reunion", "Réunion")) |
    (!is.na(population) & population %in% c("Reunion", "Réunion")) |
    (!is.na(sample_id) & grepl("^Reu", sample_id)),
  `:=`(
    country = "Reunion",
    population = "Reunion",
    macroregion_3 = "Asia",
    subregion = "SE_Asia"
  )
]

if (any(is.na(metadata_clean$sample_id))) {
  stop("Metadata contains missing sample_id values.")
}

duplicated_ids <- metadata_clean[duplicated(sample_id), unique(sample_id)]

if (length(duplicated_ids) > 0) {
  stop(
    "Metadata contains duplicated sample_id values after cleaning: ",
    paste(head(duplicated_ids, 20), collapse = ", "),
    if (length(duplicated_ids) > 20) " ..." else ""
  )
}

missing_region <- metadata_clean[is.na(macroregion_3) | is.na(subregion)]

if (nrow(missing_region) > 0) {
  warning(
    "Metadata contains ", nrow(missing_region),
    " rows with missing macroregion_3 or subregion. See metadata_missing_region.tsv."
  )
  
  fwrite(
    missing_region,
    file.path(opt$out_dir, "metadata_missing_region.tsv"),
    sep = "\t"
  )
}

fwrite(metadata_clean, file.path(opt$out_dir, "metadata_clean.tsv"), sep = "\t")
saveRDS(metadata_clean, file.path(opt$out_dir, "metadata_clean.rds"))

metadata_region_counts <- metadata_clean[
  ,
  .N,
  by = .(is_reference, macroregion_3, subregion)
][order(is_reference, macroregion_3, subregion)]

fwrite(
  metadata_region_counts,
  file.path(opt$out_dir, "metadata_region_counts.tsv"),
  sep = "\t"
)

saveRDS(
  metadata_region_counts,
  file.path(opt$out_dir, "metadata_region_counts.rds")
)

message("Metadata rows:    ", nrow(metadata_clean))
message("Metadata columns: ", paste(names(metadata_clean), collapse = ", "))

# =============================================================================
# Label consistency check against raw PHYLIP files
# =============================================================================

parse_phy_header <- function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(parts) < 2) stop("Malformed PHYLIP header: ", line)
  
  ntax <- suppressWarnings(as.integer(parts[1]))
  nsites <- suppressWarnings(as.integer(parts[2]))
  
  if (is.na(ntax) || is.na(nsites)) stop("Invalid PHYLIP header: ", line)
  
  list(ntax = ntax, nsites = nsites)
}

get_phy_first_block_labels <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- sub("\r$", "", x)
  
  nonempty <- which(nzchar(trimws(x)))
  if (length(nonempty) == 0) stop("Empty PHYLIP file")
  
  header_idx <- nonempty[1]
  header <- parse_phy_header(trimws(sub("^\ufeff", "", x[header_idx])))
  
  body <- x[(header_idx + 1L):length(x)]
  body <- body[nzchar(trimws(body))]
  
  if (length(body) < header$ntax) {
    stop(
      basename(path),
      ": fewer body lines than ntax. ntax=",
      header$ntax,
      "; body lines=",
      length(body)
    )
  }
  
  first_block <- body[seq_len(header$ntax)]
  labels <- sub("^\\s*(\\S+).*$", "\\1", first_block)
  labels <- clean_text(labels)
  
  list(
    labels = labels,
    ntax = header$ntax,
    nsites = header$nsites,
    body_lines = length(body)
  )
}

run_phy_label_check <- function(phy_files, metadata_ids, out_dir, allowed_extra_labels) {
  message("Running label consistency check on ", length(phy_files), " PHYLIP file(s)...")
  
  one <- function(f) {
    res <- tryCatch(get_phy_first_block_labels(f), error = function(e) e)
    
    og <- tools::file_path_sans_ext(basename(f))
    
    if (inherits(res, "error")) {
      return(list(
        by_og = data.table(
          og = og,
          file = f,
          ntax_header = NA_integer_,
          nsites_header = NA_integer_,
          n_labels = NA_integer_,
          n_labels_in_metadata = NA_integer_,
          n_labels_not_in_metadata = NA_integer_,
          n_labels_not_in_metadata_unexpected = NA_integer_,
          n_metadata_missing_from_og = NA_integer_,
          status = "ERROR",
          error = conditionMessage(res)
        ),
        labels_not_in_metadata = data.table(
          og = og,
          label = character(),
          allowed_extra = logical()
        ),
        metadata_missing = data.table(
          og = og,
          sample_id = character()
        )
      ))
    }
    
    labels <- unique(clean_text(res$labels))
    labels <- labels[!is.na(labels)]
    
    labels_not_meta <- setdiff(labels, metadata_ids)
    
    labels_not_meta_dt <- data.table(
      og = og,
      label = labels_not_meta
    )
    
    if (nrow(labels_not_meta_dt) > 0) {
      labels_not_meta_dt[, allowed_extra := label %in% allowed_extra_labels]
    } else {
      labels_not_meta_dt[, allowed_extra := logical()]
    }
    
    unexpected_not_meta <- labels_not_meta_dt[allowed_extra != TRUE, label]
    meta_not_labels <- setdiff(metadata_ids, labels)
    
    list(
      by_og = data.table(
        og = og,
        file = f,
        ntax_header = res$ntax,
        nsites_header = res$nsites,
        n_labels = length(labels),
        n_labels_in_metadata = sum(labels %in% metadata_ids),
        n_labels_not_in_metadata = length(labels_not_meta),
        n_labels_not_in_metadata_unexpected = length(unexpected_not_meta),
        n_metadata_missing_from_og = length(meta_not_labels),
        status = "OK",
        error = NA_character_
      ),
      labels_not_in_metadata = labels_not_meta_dt,
      metadata_missing = data.table(
        og = og,
        sample_id = meta_not_labels
      )
    )
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(opt$cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    parallel::clusterEvalQ(cl, {
      library(data.table)
    })
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "metadata_ids",
        "allowed_extra_labels",
        "clean_text",
        "parse_phy_header",
        "get_phy_first_block_labels",
        "one"
      ),
      envir = environment()
    )
    
    res <- parallel::parLapply(cl, phy_files, one)
    
  } else {
    res <- parallel::mclapply(phy_files, one, mc.cores = opt$cores)
  }
  
  by_og <- rbindlist(lapply(res, `[[`, "by_og"), fill = TRUE)
  labels_not_meta <- rbindlist(lapply(res, `[[`, "labels_not_in_metadata"), fill = TRUE)
  metadata_missing <- rbindlist(lapply(res, `[[`, "metadata_missing"), fill = TRUE)
  
  if (!("allowed_extra" %in% names(labels_not_meta))) {
    labels_not_meta[, allowed_extra := logical()]
  }
  
  labels_not_meta_unexpected <- labels_not_meta[allowed_extra != TRUE]
  
  n_ok_ogs <- by_og[status == "OK", .N]
  
  if (nrow(metadata_missing) > 0) {
    sample_summary <- data.table(sample_id = metadata_ids)
    sample_summary[, n_ogs_total := n_ok_ogs]
    sample_summary[, n_ogs_missing := vapply(sample_id, function(id) {
      sum(metadata_missing$sample_id == id)
    }, integer(1))]
    sample_summary[, n_ogs_seen := n_ogs_total - n_ogs_missing]
    sample_summary[, frac_ogs_seen := fifelse(n_ogs_total > 0, n_ogs_seen / n_ogs_total, NA_real_)]
  } else {
    sample_summary <- data.table(
      sample_id = metadata_ids,
      n_ogs_total = n_ok_ogs,
      n_ogs_missing = 0L,
      n_ogs_seen = n_ok_ogs,
      frac_ogs_seen = 1
    )
  }
  
  fwrite(by_og, file.path(out_dir, "label_consistency_by_og_phy.tsv"), sep = "\t")
  saveRDS(by_og, file.path(out_dir, "label_consistency_by_og_phy.rds"))
  
  fwrite(sample_summary, file.path(out_dir, "label_consistency_sample_summary_phy.tsv"), sep = "\t")
  saveRDS(sample_summary, file.path(out_dir, "label_consistency_sample_summary_phy.rds"))
  
  fwrite(labels_not_meta, file.path(out_dir, "label_consistency_labels_not_in_metadata_phy.tsv"), sep = "\t")
  fwrite(labels_not_meta_unexpected, file.path(out_dir, "label_consistency_labels_not_in_metadata_unexpected_phy.tsv"), sep = "\t")
  fwrite(metadata_missing, file.path(out_dir, "label_consistency_metadata_not_seen_phy.tsv"), sep = "\t")
  
  report_lines <- c(
    "PESTFLY Step 1-2 label consistency report",
    "==========================================",
    "",
    paste0("PHYLIP files checked: ", nrow(by_og)),
    paste0("Files with errors: ", by_og[status == "ERROR", .N]),
    paste0("Allowed extra labels: ", paste(allowed_extra_labels, collapse = ", ")),
    "",
    paste0("Unique labels not in metadata, including allowed extras: ", uniqueN(labels_not_meta$label)),
    paste0("Unique labels not in metadata, unexpected only: ", uniqueN(labels_not_meta_unexpected$label)),
    paste0("Metadata samples never seen in checked PHYLIP files: ", sample_summary[n_ogs_seen == 0, .N]),
    "",
    "Interpretation:",
    "- Allowed extra labels are expected labels present in alignments but absent from metadata.",
    "- In this dataset, Bdors and Blati are expected OMA/read2tree outgroups.",
    "- Unexpected labels not in metadata should be resolved before Step 3.",
    "- Metadata samples never seen in PHYLIP should also be resolved before Step 3."
  )
  
  writeLines(report_lines, file.path(out_dir, "label_consistency_report_phy.txt"))
  
  invisible(list(
    by_og = by_og,
    sample_summary = sample_summary,
    labels_not_in_metadata = labels_not_meta,
    labels_not_in_metadata_unexpected = labels_not_meta_unexpected,
    metadata_missing = metadata_missing
  ))
}

label_check <- tolower(opt$label_check)

if (!(label_check %in% c("all", "first", "none"))) {
  stop("--label_check must be one of: all, first, none")
}

if (label_check != "none") {
  phy_files <- sort(Sys.glob(file.path(opt$phy_dir, opt$phy_glob)))
  
  if (length(phy_files) == 0) {
    warning("No PHYLIP files found for label check in: ", opt$phy_dir)
  } else {
    if (label_check == "first") {
      phy_files <- phy_files[1]
    }
    
    run_phy_label_check(
      phy_files = phy_files,
      metadata_ids = metadata_clean$sample_id,
      out_dir = opt$out_dir,
      allowed_extra_labels = allowed_extra_labels
    )
  }
}

# =============================================================================
# Step 2: Ortholog FASTA QC
# =============================================================================

fasta_files <- sort(Sys.glob(file.path(opt$align_dir, opt$file_glob)))

if (length(fasta_files) == 0) {
  stop("No FASTA files found in: ", opt$align_dir)
}

if (opt$debug_n > 0L) {
  fasta_files <- head(fasta_files, opt$debug_n)
}

message("FASTA files for QC: ", length(fasta_files))

qc_one_og <- function(fasta_file, metadata_clean, opt, allowed_extra_labels) {
  og <- tools::file_path_sans_ext(basename(fasta_file))
  
  dna <- Biostrings::readDNAStringSet(fasta_file)
  taxa <- clean_text(names(dna))
  seqs <- as.character(dna)
  
  if (length(seqs) == 0L) {
    return(list(
      og_qc = data.table(
        og = og,
        file = fasta_file,
        status = "FAIL",
        fail_reason = "empty_fasta",
        n_taxa = 0L,
        n_metadata_taxa = 0L,
        aln_len = NA_integer_,
        mean_nonmissing_frac = NA_real_,
        mean_ambig_frac = NA_real_,
        n_pops_total = NA_integer_,
        n_pops_passing_occupancy = NA_integer_,
        frac_pops_passing_occupancy = NA_real_
      ),
      sample_missingness = data.table()
    ))
  }
  
  if (nzchar(opt$exclude_regex) && opt$exclude_regex != "^$") {
    keep <- !grepl(opt$exclude_regex, taxa)
    taxa <- taxa[keep]
    seqs <- seqs[keep]
  }
  
  if (length(allowed_extra_labels) > 0) {
    keep_extra <- !(taxa %in% allowed_extra_labels)
    taxa <- taxa[keep_extra]
    seqs <- seqs[keep_extra]
  }
  
  metadata_ids <- metadata_clean$sample_id
  keep_meta <- taxa %in% metadata_ids
  
  taxa <- taxa[keep_meta]
  seqs <- seqs[keep_meta]
  
  if (length(seqs) == 0L) {
    return(list(
      og_qc = data.table(
        og = og,
        file = fasta_file,
        status = "FAIL",
        fail_reason = "no_metadata_taxa",
        n_taxa = length(dna),
        n_metadata_taxa = 0L,
        aln_len = NA_integer_,
        mean_nonmissing_frac = NA_real_,
        mean_ambig_frac = NA_real_,
        n_pops_total = NA_integer_,
        n_pops_passing_occupancy = NA_integer_,
        frac_pops_passing_occupancy = NA_real_
      ),
      sample_missingness = data.table()
    ))
  }
  
  lens <- nchar(seqs)
  aln_len <- max(lens)
  
  if (length(unique(lens)) != 1L) {
    max_len <- max(lens)
    seqs <- ifelse(
      lens < max_len,
      paste0(seqs, strrep("-", max_len - lens)),
      substr(seqs, 1L, max_len)
    )
    aln_len <- max_len
  }
  
  chars_list <- strsplit(seqs, "", fixed = TRUE)
  
  nonmissing_count <- vapply(chars_list, function(z) {
    sum(z %in% c("A", "C", "G", "T"))
  }, integer(1))
  
  ambig_count <- vapply(chars_list, function(z) {
    sum(!(z %in% c("A", "C", "G", "T", "-", "N", "?")))
  }, integer(1))
  
  nonmissing_frac <- nonmissing_count / aln_len
  ambig_frac <- ambig_count / aln_len
  
  sm <- data.table(
    og = og,
    sample_id = taxa,
    aln_len = aln_len,
    nonmissing_count = nonmissing_count,
    nonmissing_frac = nonmissing_frac,
    ambig_count = ambig_count,
    ambig_frac = ambig_frac
  )
  
  meta_cols <- intersect(
    c("sample_id", "is_reference", "population", "country", "macroregion_3", "subregion"),
    names(metadata_clean)
  )
  
  sm <- merge(
    sm,
    metadata_clean[, ..meta_cols],
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  pop_occ <- sm[
    !is.na(population),
    .(
      n_samples_in_og = .N,
      n_samples_nonmissing_ok = sum(nonmissing_frac >= opt$min_nonmissing_frac, na.rm = TRUE),
      occupancy = mean(nonmissing_frac >= opt$min_nonmissing_frac, na.rm = TRUE)
    ),
    by = population
  ]
  
  n_pops_total <- nrow(pop_occ)
  n_pops_passing <- pop_occ[occupancy >= opt$min_pop_occupancy, .N]
  frac_pops_passing <- if (n_pops_total > 0) n_pops_passing / n_pops_total else NA_real_
  
  mean_nonmissing_frac <- mean(sm$nonmissing_frac, na.rm = TRUE)
  mean_ambig_frac <- mean(sm$ambig_frac, na.rm = TRUE)
  
  fail_reasons <- character()
  
  if (aln_len < opt$min_len) {
    fail_reasons <- c(fail_reasons, "too_short")
  }
  
  if (is.na(mean_ambig_frac) || mean_ambig_frac > opt$max_mean_ambig) {
    fail_reasons <- c(fail_reasons, "high_mean_ambig")
  }
  
  if (is.na(frac_pops_passing) || frac_pops_passing < opt$min_pop_occupancy_frac_pops) {
    fail_reasons <- c(fail_reasons, "low_population_occupancy")
  }
  
  status <- if (length(fail_reasons) == 0L) "PASS" else "FAIL"
  fail_reason <- if (length(fail_reasons) == 0L) NA_character_ else paste(fail_reasons, collapse = ";")
  
  og_qc <- data.table(
    og = og,
    file = fasta_file,
    status = status,
    fail_reason = fail_reason,
    n_taxa = length(dna),
    n_metadata_taxa = length(taxa),
    aln_len = aln_len,
    mean_nonmissing_frac = mean_nonmissing_frac,
    mean_ambig_frac = mean_ambig_frac,
    n_pops_total = n_pops_total,
    n_pops_passing_occupancy = n_pops_passing,
    frac_pops_passing_occupancy = frac_pops_passing
  )
  
  list(
    og_qc = og_qc,
    sample_missingness = sm
  )
}

qc_wrapper <- function(f) {
  tryCatch(
    qc_one_og(
      fasta_file = f,
      metadata_clean = metadata_clean,
      opt = opt,
      allowed_extra_labels = allowed_extra_labels
    ),
    error = function(e) {
      og <- tools::file_path_sans_ext(basename(f))
      
      list(
        og_qc = data.table(
          og = og,
          file = f,
          status = "FAIL",
          fail_reason = paste0("ERROR: ", conditionMessage(e)),
          n_taxa = NA_integer_,
          n_metadata_taxa = NA_integer_,
          aln_len = NA_integer_,
          mean_nonmissing_frac = NA_real_,
          mean_ambig_frac = NA_real_,
          n_pops_total = NA_integer_,
          n_pops_passing_occupancy = NA_integer_,
          frac_pops_passing_occupancy = NA_real_
        ),
        sample_missingness = data.table()
      )
    }
  )
}

message("Running OG-level QC in parallel...")

if (.Platform$OS.type == "windows") {
  cl <- parallel::makeCluster(opt$cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(data.table)
    library(Biostrings)
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "metadata_clean",
      "opt",
      "allowed_extra_labels",
      "clean_text",
      "qc_one_og",
      "qc_wrapper"
    ),
    envir = environment()
  )
  
  qc_results <- parallel::parLapply(cl, fasta_files, qc_wrapper)
  
} else {
  qc_results <- parallel::mclapply(
    fasta_files,
    qc_wrapper,
    mc.cores = opt$cores
  )
}

og_qc <- rbindlist(lapply(qc_results, `[[`, "og_qc"), fill = TRUE)

setorder(og_qc, og)

fwrite(og_qc, file.path(opt$out_dir, "og_qc.tsv"), sep = "\t")
saveRDS(og_qc, file.path(opt$out_dir, "og_qc.rds"))

ogs_pass <- og_qc[status == "PASS", og]
ogs_fail <- og_qc[status != "PASS", og]

writeLines(ogs_pass, file.path(opt$out_dir, "ogs_pass_qc.txt"))
writeLines(ogs_fail, file.path(opt$out_dir, "ogs_fail_qc.txt"))

saveRDS(ogs_pass, file.path(opt$out_dir, "ogs_pass_qc.rds"))
saveRDS(ogs_fail, file.path(opt$out_dir, "ogs_fail_qc.rds"))

if (!isTRUE(opt$no_sample_missingness)) {
  sample_missingness <- rbindlist(
    lapply(qc_results, `[[`, "sample_missingness"),
    fill = TRUE
  )
  
  fwrite(sample_missingness, file.path(opt$out_dir, "sample_missingness.tsv"), sep = "\t")
  saveRDS(sample_missingness, file.path(opt$out_dir, "sample_missingness.rds"))
  
  if ("is_reference" %in% names(sample_missingness)) {
    reference_nonmissing_by_og <- sample_missingness[
      is_reference %in% TRUE,
      .(
        n_reference_samples = .N,
        mean_reference_nonmissing_frac = mean(nonmissing_frac, na.rm = TRUE),
        min_reference_nonmissing_frac = min(nonmissing_frac, na.rm = TRUE),
        max_reference_nonmissing_frac = max(nonmissing_frac, na.rm = TRUE)
      ),
      by = og
    ]
    
    saveRDS(
      reference_nonmissing_by_og,
      file.path(opt$out_dir, "reference_nonmissing_by_og.rds")
    )
  }
}

# =============================================================================
# Save parameters and run info
# =============================================================================

qc_params <- data.table(
  parameter = c(
    "align_dir",
    "file_glob",
    "metadata_xlsx",
    "metadata_sheet",
    "phy_dir",
    "phy_glob",
    "min_len",
    "min_nonmissing_frac",
    "min_pop_occupancy",
    "min_pop_occupancy_frac_pops",
    "max_mean_ambig",
    "exclude_regex",
    "allowed_extra_labels",
    "cores",
    "debug_n",
    "label_check"
  ),
  value = c(
    opt$align_dir,
    opt$file_glob,
    opt$metadata_xlsx,
    opt$metadata_sheet,
    opt$phy_dir,
    opt$phy_glob,
    as.character(opt$min_len),
    as.character(opt$min_nonmissing_frac),
    as.character(opt$min_pop_occupancy),
    as.character(opt$min_pop_occupancy_frac_pops),
    as.character(opt$max_mean_ambig),
    opt$exclude_regex,
    paste(allowed_extra_labels, collapse = ","),
    as.character(opt$cores),
    as.character(opt$debug_n),
    opt$label_check
  )
)

fwrite(qc_params, file.path(opt$out_dir, "qc_params.tsv"), sep = "\t")
saveRDS(qc_params, file.path(opt$out_dir, "qc_params.rds"))

run_info <- list(
  step = "step1_2_metadata_snapshot_orthologqc",
  repo_root = repo_root,
  options = opt,
  allowed_extra_labels = allowed_extra_labels,
  n_metadata_rows = nrow(metadata_clean),
  n_fasta_files = length(fasta_files),
  n_ogs_pass = length(ogs_pass),
  n_ogs_fail = length(ogs_fail),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(run_info, file.path(opt$out_dir, "step1_2_run_info.rds"))

# =============================================================================
# Final report
# =============================================================================

message("Done.")
message("Metadata written to: ", file.path(opt$out_dir, "metadata_clean.tsv"))
message("Region counts:       ", file.path(opt$out_dir, "metadata_region_counts.tsv"))
message("OG QC written to:    ", file.path(opt$out_dir, "og_qc.tsv"))
message("PASS OGs:            ", length(ogs_pass))
message("FAIL OGs:            ", length(ogs_fail))
message("R objects written to: ", opt$out_dir)

if (file.exists(file.path(opt$out_dir, "label_consistency_report_phy.txt"))) {
  message("Label report:         ", file.path(opt$out_dir, "label_consistency_report_phy.txt"))
}