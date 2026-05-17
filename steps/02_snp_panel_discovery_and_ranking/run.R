#!/usr/bin/env Rscript

# =============================================================================
# PESTFLY — Step 3: Multi-panel SNP extraction and ranking
# =============================================================================
#
# Purpose
# -------
# Select informative biallelic SNPs from FASTA alignments that passed Step 1–2 QC.
#
# Default panels:
#   P1: macroregion_3 Africa vs Asia
#   P2: subregion within Africa
#   P3: subregion within Asia
#
# Optional panel:
#   P4: population-level global panel
#
# P4 is not run by default because it is computationally heavier and more
# exploratory. To run it, use:
#
#   --run_p4
#
# Inputs
# ------
#   results/00_fasta/OG*.fasta
#   results/01_qc/metadata_clean.tsv
#   results/01_qc/ogs_pass_qc.txt
#
# Outputs
# -------
#   results/02_snp_panels/panels_index.tsv
#   results/02_snp_panels/panels_index.rds
#   results/02_snp_panels/snp_map.tsv
#   results/02_snp_panels/snp_map.rds
#   results/02_snp_panels/step3_run_info.rds
#
#   results/02_snp_panels/panels/<panel_id>/
#     reference_group_counts.tsv
#     reference_group_counts.rds
#     per_og_stats.tsv
#     per_og_stats.rds
#     snp_map.tsv
#     snp_map.rds
#     snp_matrix.rds
#     snp_matrix.tsv.gz       optional, if --write_tsv_matrix
#     og_panel_topK_<K>.txt
#     panel_run_info.rds
#
# Run default P1-P3:
#   Rscript steps/step3_snp_extraction_ranking/run.R --cores 8
#
# Test:
#   Rscript steps/step3_snp_extraction_ranking/run.R --debug_n 200 --cores 4
#
# Optional P4:
#   Rscript steps/step3_snp_extraction_ranking/run.R --cores 8 --run_p4
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
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
    help = "Directory with Step 0 FASTA alignments [default %default]"
  ),
  make_option(
    "--file_glob",
    type = "character",
    default = "OG*.fasta",
    help = "FASTA file pattern [default %default]"
  ),
  make_option(
    "--qc_dir",
    type = "character",
    default = "results/01_qc",
    help = "Step 1–2 output directory [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Step 3 output directory [default %default]"
  ),
  make_option(
    "--reference_value",
    type = "character",
    default = "auto",
    help = "Which is_reference value means reference DB: auto, TRUE, or FALSE [default %default]"
  ),
  make_option(
    "--macroregions",
    type = "character",
    default = "Africa,Asia",
    help = "Comma-separated macroregion_3 values used for P1/P2/P3 [default %default]"
  ),
  make_option(
    "--max_site_missing",
    type = "double",
    default = 0.20,
    help = "Maximum missing fraction among reference samples at a site [default %default]"
  ),
  make_option(
    "--min_mac",
    type = "integer",
    default = 2,
    help = "Minimum minor allele count among reference samples [default %default]"
  ),
  make_option(
    "--min_ref_per_group",
    type = "integer",
    default = 3,
    help = "Minimum reference samples per class/group for a panel [default %default]"
  ),
  make_option(
    "--max_snps_per_og",
    type = "integer",
    default = 1,
    help = "Maximum SNPs retained per OG per panel [default %default]"
  ),
  make_option(
    "--exclude_regex",
    type = "character",
    default = "^$",
    help = "Regex of taxa to exclude before SNP extraction, if needed [default %default]"
  ),
  make_option(
    "--allowed_extra_labels",
    type = "character",
    default = "Bdors,Blati",
    help = "Comma-separated labels allowed in FASTA but absent from metadata [default %default]"
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 0,
    help = "Parallel workers. 0 = all detected cores minus one [default %default]"
  ),
  make_option(
    "--chunk_size",
    type = "integer",
    default = 300,
    help = "Number of FASTA files per processing chunk [default %default]"
  ),
  make_option(
    "--debug_n",
    type = "integer",
    default = 0,
    help = "Process first N passing OGs only; 0 = all [default %default]"
  ),
  make_option(
    "--panel_sizes",
    type = "character",
    default = "20,50,100,200,500,1000,2000,5000",
    help = "Comma-separated Top-K values for og_panel_topK files [default %default]"
  ),
  make_option(
    "--write_tsv_matrix",
    action = "store_true",
    default = FALSE,
    help = "Also write snp_matrix.tsv.gz for each panel [default %default]"
  ),
  make_option(
    "--run_p4",
    action = "store_true",
    default = FALSE,
    help = "Also run optional global population-level panel P4 [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$align_dir <- resolve_path(opt$align_dir)
opt$qc_dir <- resolve_path(opt$qc_dir)
opt$out_dir <- resolve_path(opt$out_dir)

if (opt$cores <= 0) {
  opt$cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

panel_sizes <- as.integer(trimws(unlist(strsplit(opt$panel_sizes, ","))))
panel_sizes <- panel_sizes[!is.na(panel_sizes) & panel_sizes > 0]

allowed_extra_labels <- trimws(unlist(strsplit(opt$allowed_extra_labels, ",")))
allowed_extra_labels <- allowed_extra_labels[nzchar(allowed_extra_labels)]

macroregions <- trimws(unlist(strsplit(opt$macroregions, ",")))
macroregions <- macroregions[nzchar(macroregions)]

if (length(macroregions) != 2L) {
  stop("--macroregions must contain exactly two values, e.g. Africa,Asia")
}

message("Repo root:            ", repo_root)
message("FASTA dir:            ", opt$align_dir)
message("QC dir:               ", opt$qc_dir)
message("Output dir:           ", opt$out_dir)
message("Cores:                ", opt$cores)
message("Reference value:      ", opt$reference_value)
message("Macroregions:         ", paste(macroregions, collapse = ", "))
message("Allowed extra labels: ", paste(allowed_extra_labels, collapse = ", "))
message("Run P4:               ", isTRUE(opt$run_p4))

# =============================================================================
# Helpers
# =============================================================================

clean_text <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("[\u200B-\u200D\uFEFF]", "", x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  x
}

as_clean_logical <- function(x) {
  if (is.logical(x)) return(x)
  
  y <- clean_text(x)
  y <- tolower(y)
  
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE
  
  out
}

is_missing_base <- function(x) {
  is.na(x) | !(x %in% c("A", "C", "G", "T"))
}

safe_aln_matrix <- function(aln) {
  m <- as.matrix(aln)
  
  if (nrow(m) == length(aln)) {
    rownames(m) <- clean_text(names(aln))
    return(m)
  }
  
  if (ncol(m) == length(aln)) {
    m <- t(m)
    rownames(m) <- clean_text(names(aln))
    return(m)
  }
  
  stop(
    "Unexpected alignment matrix orientation. dim=",
    paste(dim(m), collapse = "x"),
    "; n_seq=",
    length(aln)
  )
}

hudson_fst <- function(p1, p2, n1, n2) {
  if (is.na(p1) || is.na(p2) || n1 < 2 || n2 < 2) return(NA_real_)
  
  h1 <- 2 * p1 * (1 - p1)
  h2 <- 2 * p2 * (1 - p2)
  
  between <- (p1 - p2)^2 - (h1 / (n1 - 1)) - (h2 / (n2 - 1))
  within <- between + h1 + h2
  
  if (!is.finite(within) || within <= 0) return(NA_real_)
  
  fst <- between / within
  fst <- max(0, min(1, fst))
  fst
}

score_binary_fst <- function(alt_allele, col, groups) {
  lev <- unique(groups[!is.na(groups)])
  if (length(lev) != 2L) return(NA_real_)
  
  idx1 <- which(groups == lev[1])
  idx2 <- which(groups == lev[2])
  
  c1 <- col[idx1]
  c2 <- col[idx2]
  
  nm1 <- !is_missing_base(c1)
  nm2 <- !is_missing_base(c2)
  
  n1 <- sum(nm1)
  n2 <- sum(nm2)
  
  if (n1 < 2 || n2 < 2) return(NA_real_)
  
  p1 <- sum(c1[nm1] == alt_allele) / n1
  p2 <- sum(c2[nm2] == alt_allele) / n2
  
  hudson_fst(p1, p2, n1, n2)
}

score_multiclass_maxpair <- function(alt_allele, col, groups) {
  lev <- unique(groups[!is.na(groups)])
  
  if (length(lev) < 2L) {
    return(list(score = NA_real_, best_pair = NA_character_))
  }
  
  best <- NA_real_
  best_pair <- NA_character_
  
  for (a in seq_len(length(lev) - 1L)) {
    for (b in (a + 1L):length(lev)) {
      g1 <- lev[a]
      g2 <- lev[b]
      
      idx1 <- which(groups == g1)
      idx2 <- which(groups == g2)
      
      c1 <- col[idx1]
      c2 <- col[idx2]
      
      nm1 <- !is_missing_base(c1)
      nm2 <- !is_missing_base(c2)
      
      n1 <- sum(nm1)
      n2 <- sum(nm2)
      
      if (n1 < 2 || n2 < 2) next
      
      p1 <- sum(c1[nm1] == alt_allele) / n1
      p2 <- sum(c2[nm2] == alt_allele) / n2
      
      fst <- hudson_fst(p1, p2, n1, n2)
      
      if (is.na(fst)) next
      
      if (is.na(best) || fst > best) {
        best <- fst
        best_pair <- paste(g1, g2, sep = "__vs__")
      }
    }
  }
  
  list(score = best, best_pair = best_pair)
}

detect_reference_value <- function(meta, eligible, reference_value) {
  reference_value <- toupper(reference_value)
  
  if (reference_value %in% c("TRUE", "T", "1", "YES")) {
    return(list(value = TRUE, note = "user_defined_TRUE"))
  }
  
  if (reference_value %in% c("FALSE", "F", "0", "NO")) {
    return(list(value = FALSE, note = "user_defined_FALSE"))
  }
  
  n_true <- sum(meta$is_reference %in% TRUE & eligible, na.rm = TRUE)
  n_false <- sum(meta$is_reference %in% FALSE & eligible, na.rm = TRUE)
  
  if (n_true >= n_false) {
    list(value = TRUE, note = paste0("auto_TRUE_n=", n_true))
  } else {
    list(value = FALSE, note = paste0("auto_FALSE_n=", n_false))
  }
}

# =============================================================================
# Load Step 1–2 outputs
# =============================================================================

metadata_path <- file.path(opt$qc_dir, "metadata_clean.tsv")
ogs_pass_path <- file.path(opt$qc_dir, "ogs_pass_qc.txt")

if (!file.exists(metadata_path)) {
  stop("Missing Step 1–2 metadata file: ", metadata_path)
}

if (!file.exists(ogs_pass_path)) {
  stop("Missing Step 1–2 passing OG list: ", ogs_pass_path)
}

meta <- fread(metadata_path)

required_meta_cols <- c("sample_id", "is_reference", "population", "macroregion_3", "subregion")
missing_meta_cols <- setdiff(required_meta_cols, names(meta))

if (length(missing_meta_cols) > 0L) {
  stop(
    "metadata_clean.tsv is missing required column(s): ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

for (nm in names(meta)) {
  if (!inherits(meta[[nm]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    meta[[nm]] <- clean_text(meta[[nm]])
  }
}

meta[, sample_id := clean_text(sample_id)]
meta[, is_reference := as_clean_logical(is_reference)]
meta[, population := clean_text(population)]
meta[, macroregion_3 := clean_text(macroregion_3)]
meta[, subregion := clean_text(subregion)]

if (any(is.na(meta$sample_id))) {
  stop("metadata_clean.tsv contains missing sample_id values.")
}

dup_ids <- meta[duplicated(sample_id), unique(sample_id)]

if (length(dup_ids) > 0L) {
  stop(
    "metadata_clean.tsv contains duplicated sample_id values: ",
    paste(head(dup_ids, 20), collapse = ", "),
    if (length(dup_ids) > 20) " ..." else ""
  )
}

all_samples <- meta$sample_id

ogs_pass <- clean_text(readLines(ogs_pass_path, warn = FALSE))
ogs_pass <- ogs_pass[!is.na(ogs_pass) & nzchar(ogs_pass)]

fasta_all <- sort(Sys.glob(file.path(opt$align_dir, opt$file_glob)))

if (length(fasta_all) == 0L) {
  stop("No FASTA files found in: ", opt$align_dir)
}

fasta_dt <- data.table(
  og = tools::file_path_sans_ext(basename(fasta_all)),
  file = fasta_all
)

fasta_dt <- fasta_dt[og %in% ogs_pass]

if (nrow(fasta_dt) == 0L) {
  stop("No FASTA files match passing OGs from Step 1–2.")
}

setorder(fasta_dt, og)

if (opt$debug_n > 0L) {
  fasta_dt <- head(fasta_dt, opt$debug_n)
}

message("Metadata samples: ", nrow(meta))
message("Passing OG FASTA files for Step 3: ", nrow(fasta_dt))

# =============================================================================
# Define panels
# =============================================================================

make_panels <- function(meta, macroregions, run_p4 = FALSE) {
  m1 <- macroregions[1]
  m2 <- macroregions[2]
  
  panels <- list(
    list(
      panel_id = "P1_macroregion_africa_vs_asia",
      group_col = "macroregion_3",
      panel_type = "binary",
      eligible_fun = function(meta) {
        !is.na(meta$macroregion_3) & meta$macroregion_3 %in% c(m1, m2)
      }
    ),
    list(
      panel_id = "P2_subregion_within_africa",
      group_col = "subregion",
      panel_type = "multiclass",
      eligible_fun = function(meta) {
        !is.na(meta$macroregion_3) &
          meta$macroregion_3 == m1 &
          !is.na(meta$subregion)
      }
    ),
    list(
      panel_id = "P3_subregion_within_asia",
      group_col = "subregion",
      panel_type = "multiclass",
      eligible_fun = function(meta) {
        !is.na(meta$macroregion_3) &
          meta$macroregion_3 == m2 &
          !is.na(meta$subregion)
      }
    )
  )
  
  if (isTRUE(run_p4)) {
    panels[[length(panels) + 1L]] <- list(
      panel_id = "P4_population_global",
      group_col = "population",
      panel_type = "multiclass",
      eligible_fun = function(meta) {
        !is.na(meta$population)
      }
    )
  }
  
  panels
}

panels <- make_panels(meta, macroregions, run_p4 = opt$run_p4)

# =============================================================================
# Process one OG for one panel
# =============================================================================

process_one_og_one_panel <- function(
    f,
    og_id,
    panel_id,
    panel_type,
    group_vec_named,
    exclude_regex,
    allowed_extra_labels,
    max_site_missing,
    min_mac,
    max_snps_per_og,
    all_samples
) {
  aln <- tryCatch(Biostrings::readDNAStringSet(f), error = function(e) e)
  
  if (inherits(aln, "error")) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "read_error",
        error = conditionMessage(aln),
        n_taxa_raw = NA_integer_,
        n_taxa_metadata = NA_integer_,
        n_ref = NA_integer_,
        aln_len = NA_integer_,
        n_sites_scanned = NA_integer_,
        n_sites_biallelic = NA_integer_,
        n_sites_pass = NA_integer_,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  taxa <- clean_text(names(aln))
  names(aln) <- taxa
  
  n_taxa_raw <- length(aln)
  
  if (n_taxa_raw == 0L) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "no_taxa",
        error = NA_character_,
        n_taxa_raw = 0L,
        n_taxa_metadata = 0L,
        n_ref = 0L,
        aln_len = NA_integer_,
        n_sites_scanned = 0L,
        n_sites_biallelic = 0L,
        n_sites_pass = 0L,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  if (nzchar(exclude_regex) && exclude_regex != "^$") {
    keep <- !grepl(exclude_regex, taxa, perl = TRUE)
    aln <- aln[keep]
    taxa <- clean_text(names(aln))
  }
  
  if (length(allowed_extra_labels) > 0L) {
    keep <- !(taxa %in% allowed_extra_labels)
    aln <- aln[keep]
    taxa <- clean_text(names(aln))
  }
  
  keep_meta <- taxa %in% all_samples
  aln <- aln[keep_meta]
  taxa <- clean_text(names(aln))
  
  if (length(aln) == 0L) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "no_metadata_taxa",
        error = NA_character_,
        n_taxa_raw = n_taxa_raw,
        n_taxa_metadata = 0L,
        n_ref = 0L,
        aln_len = NA_integer_,
        n_sites_scanned = 0L,
        n_sites_biallelic = 0L,
        n_sites_pass = 0L,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  mat <- tryCatch(safe_aln_matrix(aln), error = function(e) e)
  
  if (inherits(mat, "error")) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "matrix_error",
        error = conditionMessage(mat),
        n_taxa_raw = n_taxa_raw,
        n_taxa_metadata = length(aln),
        n_ref = NA_integer_,
        aln_len = NA_integer_,
        n_sites_scanned = NA_integer_,
        n_sites_biallelic = NA_integer_,
        n_sites_pass = NA_integer_,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  mat <- toupper(mat)
  aln_len <- ncol(mat)
  
  ref_taxa <- intersect(rownames(mat), names(group_vec_named))
  ref_taxa <- ref_taxa[!is.na(group_vec_named[ref_taxa])]
  
  if (length(ref_taxa) < 2L) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "too_few_ref",
        error = NA_character_,
        n_taxa_raw = n_taxa_raw,
        n_taxa_metadata = nrow(mat),
        n_ref = length(ref_taxa),
        aln_len = aln_len,
        n_sites_scanned = aln_len,
        n_sites_biallelic = 0L,
        n_sites_pass = 0L,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  ref_mat <- mat[ref_taxa, , drop = FALSE]
  ref_groups <- group_vec_named[ref_taxa]
  
  n_sites_scanned <- ncol(ref_mat)
  n_sites_biallelic <- 0L
  n_sites_pass <- 0L
  
  candidates <- list()
  
  for (pos in seq_len(ncol(ref_mat))) {
    col <- ref_mat[, pos]
    nonmiss <- !is_missing_base(col)
    n_nonmiss <- sum(nonmiss)
    
    if (n_nonmiss == 0L) next
    
    miss_frac <- 1 - (n_nonmiss / length(col))
    
    if (miss_frac > max_site_missing) next
    
    alleles <- sort(unique(col[nonmiss]))
    alleles <- alleles[alleles %in% c("A", "C", "G", "T")]
    
    if (length(alleles) != 2L) next
    
    n_sites_biallelic <- n_sites_biallelic + 1L
    
    counts <- table(factor(col[nonmiss], levels = alleles))
    mac <- min(as.integer(counts))
    
    if (mac < min_mac) next
    
    alt_allele <- names(which.min(counts))[1]
    ref_allele <- setdiff(alleles, alt_allele)[1]
    
    if (panel_type == "binary") {
      score <- score_binary_fst(alt_allele, col, ref_groups)
      best_pair <- paste(unique(ref_groups), collapse = "__vs__")
    } else {
      sc <- score_multiclass_maxpair(alt_allele, col, ref_groups)
      score <- sc$score
      best_pair <- sc$best_pair
    }
    
    if (is.na(score) || !is.finite(score)) next
    
    n_sites_pass <- n_sites_pass + 1L
    
    candidates[[length(candidates) + 1L]] <- data.table(
      panel_id = panel_id,
      og = og_id,
      file = f,
      pos = pos,
      snp_id = paste(og_id, pos, ref_allele, alt_allele, sep = "_"),
      ref_allele = ref_allele,
      alt_allele = alt_allele,
      score = score,
      best_pair = best_pair,
      missing_frac_ref = miss_frac,
      mac_ref = mac,
      n_ref_nonmissing = n_nonmiss,
      n_ref_total = length(col)
    )
  }
  
  if (length(candidates) == 0L) {
    return(list(
      stats = data.table(
        og = og_id,
        file = f,
        panel_id = panel_id,
        status = "no_snp_pass",
        error = NA_character_,
        n_taxa_raw = n_taxa_raw,
        n_taxa_metadata = nrow(mat),
        n_ref = length(ref_taxa),
        aln_len = aln_len,
        n_sites_scanned = n_sites_scanned,
        n_sites_biallelic = n_sites_biallelic,
        n_sites_pass = n_sites_pass,
        best_score = NA_real_
      ),
      snps = NULL,
      geno = NULL
    ))
  }
  
  cand_dt <- rbindlist(candidates, fill = TRUE)
  setorder(cand_dt, -score, missing_frac_ref, -mac_ref)
  
  keep_n <- min(max_snps_per_og, nrow(cand_dt))
  cand_dt <- cand_dt[seq_len(keep_n)]
  cand_dt[, rank_within_og := seq_len(.N)]
  
  geno_list <- list()
  
  for (i in seq_len(nrow(cand_dt))) {
    pos <- cand_dt$pos[i]
    refa <- cand_dt$ref_allele[i]
    alta <- cand_dt$alt_allele[i]
    sid <- cand_dt$snp_id[i]
    
    calls <- mat[, pos]
    names(calls) <- rownames(mat)
    
    g <- rep(NA_integer_, length(all_samples))
    names(g) <- all_samples
    
    ok <- !(is_missing_base(calls)) & calls %in% c(refa, alta)
    g[names(calls)[ok]] <- ifelse(calls[ok] == refa, 0L, 1L)
    
    geno_list[[sid]] <- g
  }
  
  stats <- data.table(
    og = og_id,
    file = f,
    panel_id = panel_id,
    status = "ok",
    error = NA_character_,
    n_taxa_raw = n_taxa_raw,
    n_taxa_metadata = nrow(mat),
    n_ref = length(ref_taxa),
    aln_len = aln_len,
    n_sites_scanned = n_sites_scanned,
    n_sites_biallelic = n_sites_biallelic,
    n_sites_pass = n_sites_pass,
    best_score = max(cand_dt$score, na.rm = TRUE)
  )
  
  list(
    stats = stats,
    snps = cand_dt,
    geno = geno_list
  )
}

# =============================================================================
# Parallel processing helpers
# =============================================================================

process_panel_chunk <- function(
    chunk_dt,
    panel_id,
    panel_type,
    group_vec_named,
    exclude_regex,
    allowed_extra_labels,
    max_site_missing,
    min_mac,
    max_snps_per_og,
    all_samples
) {
  out <- vector("list", nrow(chunk_dt))
  
  for (i in seq_len(nrow(chunk_dt))) {
    out[[i]] <- process_one_og_one_panel(
      f = chunk_dt$file[i],
      og_id = chunk_dt$og[i],
      panel_id = panel_id,
      panel_type = panel_type,
      group_vec_named = group_vec_named,
      exclude_regex = exclude_regex,
      allowed_extra_labels = allowed_extra_labels,
      max_site_missing = max_site_missing,
      min_mac = min_mac,
      max_snps_per_og = max_snps_per_og,
      all_samples = all_samples
    )
  }
  
  out
}

run_parallel_panel <- function(fasta_dt, panel_id, panel_type, group_vec_named) {
  idx <- split(seq_len(nrow(fasta_dt)), ceiling(seq_len(nrow(fasta_dt)) / opt$chunk_size))
  chunks <- lapply(idx, function(ii) fasta_dt[ii])
  
  message("Panel ", panel_id, ": processing ", nrow(fasta_dt), " OGs in ", length(chunks), " chunk(s).")
  
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
        "clean_text",
        "is_missing_base",
        "safe_aln_matrix",
        "hudson_fst",
        "score_binary_fst",
        "score_multiclass_maxpair",
        "process_one_og_one_panel",
        "process_panel_chunk",
        "panel_id",
        "panel_type",
        "group_vec_named",
        "opt",
        "allowed_extra_labels",
        "all_samples"
      ),
      envir = environment()
    )
    
    chunk_results <- parallel::parLapply(
      cl,
      chunks,
      function(ch) {
        process_panel_chunk(
          chunk_dt = ch,
          panel_id = panel_id,
          panel_type = panel_type,
          group_vec_named = group_vec_named,
          exclude_regex = opt$exclude_regex,
          allowed_extra_labels = allowed_extra_labels,
          max_site_missing = opt$max_site_missing,
          min_mac = opt$min_mac,
          max_snps_per_og = opt$max_snps_per_og,
          all_samples = all_samples
        )
      }
    )
    
  } else {
    chunk_results <- parallel::mclapply(
      chunks,
      function(ch) {
        process_panel_chunk(
          chunk_dt = ch,
          panel_id = panel_id,
          panel_type = panel_type,
          group_vec_named = group_vec_named,
          exclude_regex = opt$exclude_regex,
          allowed_extra_labels = allowed_extra_labels,
          max_site_missing = opt$max_site_missing,
          min_mac = opt$min_mac,
          max_snps_per_og = opt$max_snps_per_og,
          all_samples = all_samples
        )
      },
      mc.cores = opt$cores
    )
  }
  
  unlist(chunk_results, recursive = FALSE)
}

# =============================================================================
# Build genotype matrix safely
# =============================================================================

build_snp_matrix_safely <- function(geno_list_nested, snp_map, all_samples, panel_id) {
  if (length(geno_list_nested) == 0L || nrow(snp_map) == 0L) {
    return(list(
      snp_map = snp_map,
      snp_matrix = matrix(
        integer(0),
        nrow = 0,
        ncol = length(all_samples),
        dimnames = list(character(), all_samples)
      )
    ))
  }
  
  geno_flat <- list()
  
  for (lst in geno_list_nested) {
    if (is.null(lst) || length(lst) == 0L) next
    
    nms <- names(lst)
    
    if (is.null(nms)) {
      warning("Panel ", panel_id, ": found genotype sublist without names; skipping.")
      next
    }
    
    for (nm in nms) {
      if (!is.null(lst[[nm]])) {
        geno_flat[[nm]] <- lst[[nm]]
      }
    }
  }
  
  if (length(geno_flat) == 0L) {
    warning("Panel ", panel_id, ": genotype list is empty after flattening.")
    
    return(list(
      snp_map = snp_map[0],
      snp_matrix = matrix(
        integer(0),
        nrow = 0,
        ncol = length(all_samples),
        dimnames = list(character(), all_samples)
      )
    ))
  }
  
  wanted_snps <- snp_map$snp_id
  available_snps <- names(geno_flat)
  
  missing_genotypes <- setdiff(wanted_snps, available_snps)
  
  if (length(missing_genotypes) > 0L) {
    warning(
      "Panel ", panel_id, ": ",
      length(missing_genotypes),
      " SNP(s) in snp_map did not have genotype vectors. Dropping them from snp_map."
    )
    
    snp_map <- snp_map[!(snp_id %in% missing_genotypes)]
    wanted_snps <- snp_map$snp_id
  }
  
  if (length(wanted_snps) == 0L) {
    warning("Panel ", panel_id, ": no SNPs left after matching snp_map to genotype vectors.")
    
    return(list(
      snp_map = snp_map,
      snp_matrix = matrix(
        integer(0),
        nrow = 0,
        ncol = length(all_samples),
        dimnames = list(character(), all_samples)
      )
    ))
  }
  
  geno_flat <- geno_flat[wanted_snps]
  
  snp_matrix <- matrix(
    NA_integer_,
    nrow = length(wanted_snps),
    ncol = length(all_samples),
    dimnames = list(wanted_snps, all_samples)
  )
  
  for (ii in seq_along(wanted_snps)) {
    sid <- wanted_snps[ii]
    g <- geno_flat[[sid]]
    
    out <- rep(NA_integer_, length(all_samples))
    names(out) <- all_samples
    
    common <- intersect(names(g), all_samples)
    
    if (length(common) > 0L) {
      out[common] <- as.integer(g[common])
    }
    
    snp_matrix[ii, ] <- out[all_samples]
  }
  
  storage.mode(snp_matrix) <- "integer"
  
  list(
    snp_map = snp_map,
    snp_matrix = snp_matrix
  )
}

# =============================================================================
# Build panels
# =============================================================================

panels_index <- list()
global_map_list <- list()

for (panel in panels) {
  panel_id <- panel$panel_id
  group_col <- panel$group_col
  panel_type <- panel$panel_type
  
  message("\n============================================================")
  message("Building panel: ", panel_id)
  message("Group column:   ", group_col)
  message("Panel type:     ", panel_type)
  message("============================================================")
  
  panel_dir <- file.path(opt$out_dir, "panels", panel_id)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  
  eligible <- panel$eligible_fun(meta)
  ref_detect <- detect_reference_value(meta, eligible, opt$reference_value)
  ref_value <- ref_detect$value
  ref_note <- ref_detect$note
  
  ref_meta <- meta[eligible & is_reference %in% ref_value]
  ref_meta <- ref_meta[!is.na(get(group_col))]
  
  group_counts <- ref_meta[, .N, by = group_col]
  setnames(group_counts, group_col, "group")
  
  good_groups <- group_counts[N >= opt$min_ref_per_group, group]
  ref_meta <- ref_meta[get(group_col) %in% good_groups]
  
  group_counts_final <- ref_meta[, .N, by = group_col]
  setnames(group_counts_final, group_col, "group")
  
  if (nrow(group_counts_final) < 2L) {
    warning("Panel ", panel_id, " has fewer than two eligible reference groups. Skipping.")
    
    empty_index <- data.table(
      panel_id = panel_id,
      panel_dir = panel_dir,
      group_col = group_col,
      panel_type = panel_type,
      n_reference_samples = nrow(ref_meta),
      n_reference_groups = nrow(group_counts_final),
      n_snps = 0L,
      reference_value = as.character(ref_value),
      reference_note = ref_note,
      status = "SKIPPED_TOO_FEW_GROUPS"
    )
    
    panels_index[[panel_id]] <- empty_index
    fwrite(group_counts_final, file.path(panel_dir, "reference_group_counts.tsv"), sep = "\t")
    saveRDS(group_counts_final, file.path(panel_dir, "reference_group_counts.rds"))
    next
  }
  
  group_vec_named <- ref_meta[[group_col]]
  names(group_vec_named) <- ref_meta$sample_id
  
  fwrite(group_counts_final, file.path(panel_dir, "reference_group_counts.tsv"), sep = "\t")
  saveRDS(group_counts_final, file.path(panel_dir, "reference_group_counts.rds"))
  
  results <- run_parallel_panel(
    fasta_dt = fasta_dt,
    panel_id = panel_id,
    panel_type = panel_type,
    group_vec_named = group_vec_named
  )
  
  per_og_stats <- rbindlist(lapply(results, `[[`, "stats"), fill = TRUE)
  setorder(per_og_stats, og)
  
  snp_list <- lapply(results, `[[`, "snps")
  snp_list <- snp_list[!vapply(snp_list, is.null, logical(1))]
  
  geno_list_nested <- lapply(results, `[[`, "geno")
  geno_list_nested <- geno_list_nested[!vapply(geno_list_nested, is.null, logical(1))]
  
  if (length(snp_list) > 0L) {
    snp_map <- rbindlist(snp_list, fill = TRUE)
    setorder(snp_map, -score, missing_frac_ref, -mac_ref, og, pos)
    snp_map[, global_rank := seq_len(.N)]
  } else {
    snp_map <- data.table()
  }
  
  matrix_build <- build_snp_matrix_safely(
    geno_list_nested = geno_list_nested,
    snp_map = snp_map,
    all_samples = all_samples,
    panel_id = panel_id
  )
  
  snp_map <- matrix_build$snp_map
  snp_matrix <- matrix_build$snp_matrix
  
  message(
    "Panel ", panel_id,
    ": snp_map rows after matrix build = ", nrow(snp_map),
    "; snp_matrix dim = ",
    paste(dim(snp_matrix), collapse = " x ")
  )
  
  if (nrow(snp_map) > 0L) {
    setorder(snp_map, -score, missing_frac_ref, -mac_ref, og, pos)
    snp_map[, global_rank := seq_len(.N)]
  }
  
  fwrite(per_og_stats, file.path(panel_dir, "per_og_stats.tsv"), sep = "\t")
  saveRDS(per_og_stats, file.path(panel_dir, "per_og_stats.rds"))
  
  fwrite(snp_map, file.path(panel_dir, "snp_map.tsv"), sep = "\t")
  saveRDS(snp_map, file.path(panel_dir, "snp_map.rds"))
  
  saveRDS(snp_matrix, file.path(panel_dir, "snp_matrix.rds"))
  
  if (isTRUE(opt$write_tsv_matrix) && nrow(snp_matrix) > 0L) {
    matrix_dt <- data.table(snp_id = rownames(snp_matrix))
    matrix_dt <- cbind(matrix_dt, as.data.table(snp_matrix))
    fwrite(matrix_dt, gzfile(file.path(panel_dir, "snp_matrix.tsv.gz")), sep = "\t")
  }
  
  if (nrow(snp_map) > 0L) {
    rank_dt <- snp_map[order(-score, missing_frac_ref, -mac_ref)]
    
    for (K in panel_sizes) {
      K2 <- min(K, nrow(rank_dt))
      top_dt <- rank_dt[seq_len(K2)]
      writeLines(unique(top_dt$og), file.path(panel_dir, paste0("og_panel_topK_", K, ".txt")))
    }
  }
  
  panel_run_info <- list(
    panel_id = panel_id,
    group_col = group_col,
    panel_type = panel_type,
    reference_value = ref_value,
    reference_note = ref_note,
    reference_group_counts = group_counts_final,
    n_reference_samples = nrow(ref_meta),
    n_reference_groups = nrow(group_counts_final),
    n_ogs_processed = nrow(fasta_dt),
    n_snps = nrow(snp_map),
    parameters = opt,
    timestamp = Sys.time()
  )
  
  saveRDS(panel_run_info, file.path(panel_dir, "panel_run_info.rds"))
  
  panels_index[[panel_id]] <- data.table(
    panel_id = panel_id,
    panel_dir = panel_dir,
    group_col = group_col,
    panel_type = panel_type,
    n_reference_samples = nrow(ref_meta),
    n_reference_groups = nrow(group_counts_final),
    n_snps = nrow(snp_map),
    reference_value = as.character(ref_value),
    reference_note = ref_note,
    status = "OK"
  )
  
  if (nrow(snp_map) > 0L) {
    global_map_list[[panel_id]] <- snp_map
  }
  
  message("Panel ", panel_id, ": selected ", nrow(snp_map), " SNPs.")
}

# =============================================================================
# Write global outputs
# =============================================================================

panels_index_dt <- rbindlist(panels_index, fill = TRUE)
fwrite(panels_index_dt, file.path(opt$out_dir, "panels_index.tsv"), sep = "\t")
saveRDS(panels_index_dt, file.path(opt$out_dir, "panels_index.rds"))

if (length(global_map_list) > 0L) {
  global_map <- rbindlist(global_map_list, fill = TRUE)
  setorder(global_map, panel_id, global_rank)
} else {
  global_map <- data.table()
}

fwrite(global_map, file.path(opt$out_dir, "snp_map.tsv"), sep = "\t")
saveRDS(global_map, file.path(opt$out_dir, "snp_map.rds"))

run_info <- list(
  step = "step3_snp_extraction_ranking",
  repo_root = repo_root,
  align_dir = opt$align_dir,
  qc_dir = opt$qc_dir,
  out_dir = opt$out_dir,
  n_metadata_samples = nrow(meta),
  n_ogs_pass_step1_2 = length(ogs_pass),
  n_ogs_processed = nrow(fasta_dt),
  panels_index = panels_index_dt,
  parameters = opt,
  allowed_extra_labels = allowed_extra_labels,
  macroregions = macroregions,
  run_p4 = isTRUE(opt$run_p4),
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(run_info, file.path(opt$out_dir, "step3_run_info.rds"))

message("\nDone.")
message("Panels index:   ", file.path(opt$out_dir, "panels_index.tsv"))
message("Global SNP map: ", file.path(opt$out_dir, "snp_map.tsv"))
message("Step 3 outputs: ", opt$out_dir)