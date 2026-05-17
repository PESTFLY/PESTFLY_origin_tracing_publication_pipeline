#!/usr/bin/env Rscript

# =============================================================================
# PESTFLY — Step 4b: Leave-one-out validation of hierarchical assignment panels
# =============================================================================
#
# Purpose
# -------
# Validate the Step 4 assignment logic on the reference database itself using
# leave-one-out cross-validation (LOO-CV).
#
# For each reference sample:
#
#   1. remove the sample from the reference database;
#   2. estimate allele frequencies from the remaining reference samples;
#   3. assign the held-out sample as if it were unknown;
#   4. compare the predicted origin to the known metadata origin.
#
# This gives empirical estimates of:
#
#   - expected accuracy;
#   - uncertainty rate;
#   - class-specific confusion;
#   - failure modes such as too few usable SNPs;
#   - stability of calls across increasing Top-K SNP subsets.
#
# =============================================================================
# Conceptual relationship to Step 4
# =============================================================================
#
# Step 4 assigns query/intercept samples using a hierarchical multi-K framework:
#
#   P1: macroregion, Africa vs Asia
#   P2: subregion within Africa
#   P3: subregion within Asia
#
# Step 4b validates the same likelihood model and the same threshold logic on
# known reference samples.
#
# This updated Step 4b mirrors the improved Step 4 logic:
#
#   - all K values are tested and written to raw output;
#   - only usable K values enter stability calculations;
#   - usable K values are those with enough non-missing SNPs;
#   - final reliability requires posterior support, posterior gap, enough usable
#     K values, global K agreement, and large-K tail agreement.
#
# Small K values are therefore kept for transparency, but do not penalise
# stability if they fall below the minimum SNP threshold.
#
# =============================================================================
# Default panels
# =============================================================================
#
#   P1_macroregion_africa_vs_asia
#      group_col = macroregion_3
#      references = Africa + Asia
#      min SNPs = 200
#      min gap = 0.20
#
#   P2_subregion_within_africa
#      group_col = subregion
#      references = macroregion_3 == Africa
#      min SNPs = 300
#      min gap = 0.25
#
#   P3_subregion_within_asia
#      group_col = subregion
#      references = macroregion_3 == Asia
#      min SNPs = 300
#      min gap = 0.25
#
# =============================================================================
# Inputs
# =============================================================================
#
#   results/01_qc/metadata_clean.tsv
#   results/02_snp_panels/panels/<panel_id>/snp_map.tsv
#   results/02_snp_panels/panels/<panel_id>/snp_matrix.rds
#
# =============================================================================
# Outputs
# =============================================================================
#
# Written to:
#
#   results/05_loo_validation/
#
# Main outputs:
#
#   step4b_loo_panels_index.tsv / .rds
#   step4b_loo_summary_all_panels.tsv / .rds
#   step4b_loo_predictions_all_panels.tsv / .rds
#   step4b_loo_raw_all_panels.tsv.gz
#   step4b_K_grid_used.tsv / .rds
#   step4b_loo_validation_summary.xlsx
#   step4b_run_info.rds
#
# Per-panel outputs:
#
#   <prefix>_loo_predictions.tsv / .rds
#   <prefix>_loo_raw.tsv.gz
#   <prefix>_loo_summary.tsv / .rds
#   <prefix>_loo_perclass.tsv / .rds
#   <prefix>_loo_confusion.tsv / .rds
#   <prefix>_loo_reason_counts.tsv
#   FIG_<prefix>_validation_A4.pdf
#
# =============================================================================
# Run
# =============================================================================
#
#   Rscript steps/step4b_loo_validation/run.R
#
# Optional:
#
#   Rscript steps/step4b_loo_validation/run.R \
#     --tail_fraction 0.50 \
#     --min_tail_agreement 1.00 \
#     --min_agreement 0.90
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(openxlsx)
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
    "--qc_dir",
    type = "character",
    default = "results/01_qc",
    help = "Step 1-2 output directory containing metadata_clean.tsv [default %default]"
  ),
  make_option(
    "--step3_dir",
    type = "character",
    default = "results/02_snp_panels",
    help = "Step 3 output directory containing panels/ [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/05_loo_validation",
    help = "Step 4b output directory [default %default]"
  ),
  make_option(
    "--panels",
    type = "character",
    default = "P1_macroregion_africa_vs_asia,P2_subregion_within_africa,P3_subregion_within_asia",
    help = "Comma-separated panel IDs to validate [default %default]"
  ),
  
  # Minimum usable SNPs
  make_option(
    "--min_snps_query_macroregion",
    type = "integer",
    default = 200,
    help = "Minimum non-missing SNPs required for P1 macroregion LOO assignment [default %default]"
  ),
  make_option(
    "--min_snps_query_subregion",
    type = "integer",
    default = 300,
    help = "Minimum non-missing SNPs required for P2/P3 subregion LOO assignment [default %default]"
  ),
  
  # Reference-class filter
  make_option(
    "--min_ref_per_group",
    type = "integer",
    default = 3,
    help = "Minimum reference samples required per class before LOO [default %default]"
  ),
  
  # SNP filter
  make_option(
    "--keep_snps_frac",
    type = "double",
    default = 0.70,
    help = "Keep SNPs with at least this non-missing fraction across panel references [default %default]"
  ),
  
  # Likelihood model
  make_option(
    "--pseudocount",
    type = "double",
    default = 0.5,
    help = "Laplace smoothing pseudocount for reference allele frequencies [default %default]"
  ),
  make_option(
    "--epsilon",
    type = "double",
    default = 0.02,
    help = "Per-site genotype error/flip probability [default %default]"
  ),
  
  # Posterior thresholds
  make_option(
    "--high_posterior",
    type = "double",
    default = 0.95,
    help = "High posterior threshold [default %default]"
  ),
  make_option(
    "--moderate_posterior",
    type = "double",
    default = 0.85,
    help = "Moderate posterior lower bound [default %default]"
  ),
  make_option(
    "--min_gap_macroregion",
    type = "double",
    default = 0.20,
    help = "Minimum top-second posterior gap for P1 [default %default]"
  ),
  make_option(
    "--min_gap_subregion",
    type = "double",
    default = 0.25,
    help = "Minimum top-second posterior gap for P2/P3 [default %default]"
  ),
  
  # Multi-K grid
  make_option(
    "--base_K",
    type = "character",
    default = "1,2,3,4,5,10,20,50,100,200,500",
    help = "Comma-separated base K values [default %default]"
  ),
  make_option(
    "--extra_K",
    type = "character",
    default = "1000,2000,5000,10000",
    help = "Comma-separated extra K values tested when available [default %default]"
  ),
  make_option(
    "--auto_extend_K",
    type = "character",
    default = "TRUE",
    help = "Auto-extend K grid beyond base_K when enough SNPs are available: TRUE/FALSE [default %default]"
  ),
  make_option(
    "--include_all_snps_K",
    type = "character",
    default = "TRUE",
    help = "Include all available SNPs as final K value: TRUE/FALSE [default %default]"
  ),
  make_option(
    "--max_K_points",
    type = "integer",
    default = 30,
    help = "Maximum number of K points after auto-extension [default %default]"
  ),
  
  # Stability thresholds
  make_option(
    "--min_agreement",
    type = "double",
    default = 0.90,
    help = "Minimum fraction of usable K calls matching final largest-K call [default %default]"
  ),
  make_option(
    "--min_K_available",
    type = "integer",
    default = 6,
    help = "Minimum number of usable K results needed to judge stability [default %default]"
  ),
  make_option(
    "--tail_fraction",
    type = "double",
    default = 0.50,
    help = "Fraction of largest usable K values used for tail-stability check [default %default]"
  ),
  make_option(
    "--min_tail_agreement",
    type = "double",
    default = 1.00,
    help = "Required agreement with final call among largest usable-K tail values [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$qc_dir <- resolve_path(opt$qc_dir)
opt$step3_dir <- resolve_path(opt$step3_dir)
opt$out_dir <- resolve_path(opt$out_dir)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

panels_requested <- trimws(unlist(strsplit(opt$panels, ",")))
panels_requested <- panels_requested[nzchar(panels_requested)]

message("Repo root:  ", repo_root)
message("QC dir:     ", opt$qc_dir)
message("Step 3 dir: ", opt$step3_dir)
message("Output dir: ", opt$out_dir)
message("Panels:     ", paste(panels_requested, collapse = ", "))

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

parse_bool1 <- function(x, default = FALSE, argname = "") {
  if (is.null(x) || length(x) == 0L) return(as.logical(default)[1])
  
  y <- tolower(clean_text(x[1]))
  
  if (is.na(y) || !nzchar(y)) return(as.logical(default)[1])
  if (y %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (y %in% c("false", "f", "0", "no", "n")) return(FALSE)
  
  stop("Cannot parse boolean for ", argname, ": '", x, "'. Use TRUE/FALSE.")
}

parse_integer_vector <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) {
    return(integer())
  }
  
  out <- suppressWarnings(as.integer(clean_text(strsplit(x, ",")[[1]])))
  out <- out[!is.na(out) & out > 0L]
  sort(unique(out))
}

as_clean_logical <- function(x) {
  if (is.logical(x)) return(x)
  
  y <- tolower(clean_text(x))
  
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y", "reference", "ref")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n", "query", "intercept")] <- FALSE
  
  out
}

conf_label <- function(pmax, high = 0.95, mod = 0.85) {
  if (is.na(pmax)) return("NA")
  if (pmax >= high) return("High")
  if (pmax >= mod) return("Moderate")
  "Low"
}

safe_intersect <- function(x, y) {
  intersect(as.character(x), as.character(y))
}

loglik_one_group <- function(qvec, pvec, epsilon = 0.02) {
  ok <- !is.na(qvec) & !is.na(pvec)
  
  if (!any(ok)) return(NA_real_)
  
  q <- qvec[ok]
  p <- pvec[ok]
  
  p1 <- p * (1 - epsilon) + (1 - p) * epsilon
  p0 <- (1 - p) * (1 - epsilon) + p * epsilon
  
  p1 <- pmin(pmax(p1, 1e-12), 1 - 1e-12)
  p0 <- pmin(pmax(p0, 1e-12), 1 - 1e-12)
  
  sum(ifelse(q == 1, log(p1), log(p0)))
}

post_from_ll <- function(lls) {
  if (all(is.na(lls))) return(rep(NA_real_, length(lls)))
  
  m <- max(lls, na.rm = TRUE)
  w <- exp(lls - m)
  w[is.na(w)] <- 0
  
  if (sum(w) <= 0) return(rep(NA_real_, length(lls)))
  
  out <- w / sum(w)
  names(out) <- names(lls)
  out
}

make_Ks <- function(
    n_snps,
    baseK,
    extraK = c(1000L, 2000L, 5000L, 10000L),
    auto_extend = TRUE,
    include_all = TRUE,
    max_points = 30
) {
  if (is.na(n_snps) || n_snps <= 0L) {
    return(integer())
  }
  
  baseK <- sort(unique(as.integer(baseK)))
  baseK <- baseK[!is.na(baseK) & baseK > 0L]
  
  extraK <- sort(unique(as.integer(extraK)))
  extraK <- extraK[!is.na(extraK) & extraK > 0L]
  
  Ks <- baseK[baseK <= n_snps]
  
  if (isTRUE(auto_extend)) {
    Ks <- sort(unique(c(Ks, extraK[extraK <= n_snps])))
  }
  
  if (isTRUE(include_all)) {
    Ks <- sort(unique(c(Ks, as.integer(n_snps))))
  }
  
  Ks <- Ks[Ks <= n_snps]
  Ks <- sort(unique(Ks))
  
  if (length(Ks) > max_points) {
    small_K <- Ks[Ks <= 500L]
    large_K <- Ks[Ks > 500L]
    
    if (length(large_K) > 0L) {
      keep_large <- unique(round(seq(
        from = min(large_K),
        to = max(large_K),
        length.out = max(1L, max_points - length(small_K))
      )))
      
      Ks <- sort(unique(c(small_K, keep_large, max(Ks))))
    }
    
    if (length(Ks) > max_points) {
      Ks <- sort(unique(c(head(Ks, max_points - 1L), max(Ks))))
    }
  }
  
  Ks
}

stable_from_threshold <- function(dt, final_group) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]
  
  if (nrow(dt2) == 0L || is.na(final_group)) return(NA_real_)
  
  dt2 <- dt2[order(K)]
  Ks_sorted <- sort(unique(dt2$K))
  
  for (k in Ks_sorted) {
    later <- dt2[K >= k]
    if (nrow(later) > 0L && all(later$top_group == final_group)) {
      return(as.numeric(k))
    }
  }
  
  NA_real_
}

tail_agreement_with_final <- function(dt, final_group, tail_fraction = 0.50) {
  dt2 <- dt[status %in% c("ok", "low_gap") & !is.na(top_group)]
  
  if (nrow(dt2) == 0L || is.na(final_group)) {
    return(NA_real_)
  }
  
  dt2 <- dt2[order(K)]
  
  n_tail <- ceiling(nrow(dt2) * tail_fraction)
  n_tail <- max(1L, n_tail)
  
  tail_dt <- tail(dt2, n_tail)
  
  mean(tail_dt$top_group == final_group)
}

safe_sheet <- function(wb, name) {
  nm <- gsub("[\\[\\]\\*\\?/\\\\:]", "_", name)
  
  if (nchar(nm) > 31) {
    nm <- substr(nm, 1, 31)
  }
  
  existing <- names(wb)
  
  if (!(nm %in% existing)) return(nm)
  
  i <- 2
  
  repeat {
    suffix <- paste0("_", i)
    base <- substr(nm, 1, max(1, 31 - nchar(suffix)))
    nm2 <- paste0(base, suffix)
    
    if (!(nm2 %in% existing)) return(nm2)
    
    i <- i + 1
  }
}

panel_config <- function(panel_id) {
  if (panel_id == "P1_macroregion_africa_vs_asia") {
    return(list(
      panel_id = panel_id,
      group_col = "macroregion_3",
      prefix = "LOO_P1_macroregion",
      min_snps_query = opt$min_snps_query_macroregion,
      min_gap = opt$min_gap_macroregion,
      filter_fun = function(dt) dt[macroregion_3 %in% c("Africa", "Asia")]
    ))
  }
  
  if (panel_id == "P2_subregion_within_africa") {
    return(list(
      panel_id = panel_id,
      group_col = "subregion",
      prefix = "LOO_P2_africa_subregion",
      min_snps_query = opt$min_snps_query_subregion,
      min_gap = opt$min_gap_subregion,
      filter_fun = function(dt) dt[macroregion_3 == "Africa"]
    ))
  }
  
  if (panel_id == "P3_subregion_within_asia") {
    return(list(
      panel_id = panel_id,
      group_col = "subregion",
      prefix = "LOO_P3_asia_subregion",
      min_snps_query = opt$min_snps_query_subregion,
      min_gap = opt$min_gap_subregion,
      filter_fun = function(dt) dt[macroregion_3 == "Asia"]
    ))
  }
  
  stop(
    "Unknown default panel ID: ", panel_id,
    "\nKnown defaults are P1_macroregion_africa_vs_asia, ",
    "P2_subregion_within_africa, P3_subregion_within_asia."
  )
}

# =============================================================================
# Core likelihood functions
# =============================================================================

build_ref_freqs <- function(snp_mat, meta_ref, group_col, pseudocount = 0.5) {
  if (!(group_col %in% names(meta_ref))) {
    stop(
      "metadata missing group_col: ", group_col,
      "\nAvailable: ", paste(names(meta_ref), collapse = ", ")
    )
  }
  
  ref_ids <- safe_intersect(meta_ref$sample_id, colnames(snp_mat))
  
  if (length(ref_ids) < 2L) {
    return(list(freq = NULL, groups = character(0), ref_ids = ref_ids))
  }
  
  meta_ref2 <- meta_ref[sample_id %in% ref_ids]
  meta_ref2 <- meta_ref2[!is.na(get(group_col)) & clean_text(get(group_col)) != ""]
  
  if (nrow(meta_ref2) == 0L) {
    return(list(freq = NULL, groups = character(0), ref_ids = ref_ids))
  }
  
  groups <- sort(unique(meta_ref2[[group_col]]))
  groups <- groups[!is.na(groups) & clean_text(groups) != ""]
  
  if (length(groups) < 2L) {
    return(list(freq = NULL, groups = groups, ref_ids = ref_ids))
  }
  
  freq <- matrix(
    NA_real_,
    nrow = nrow(snp_mat),
    ncol = length(groups),
    dimnames = list(rownames(snp_mat), groups)
  )
  
  ids_by_group <- lapply(groups, function(g) {
    meta_ref2[get(group_col) == g, sample_id]
  })
  
  names(ids_by_group) <- groups
  
  for (g in groups) {
    ids <- safe_intersect(ids_by_group[[g]], colnames(snp_mat))
    
    if (length(ids) == 0L) next
    
    X <- snp_mat[, ids, drop = FALSE]
    
    alt_ct <- rowSums(X == 1, na.rm = TRUE)
    ref_ct <- rowSums(X == 0, na.rm = TRUE)
    tot_ct <- alt_ct + ref_ct
    
    p <- rep(NA_real_, length(tot_ct))
    ok <- tot_ct > 0
    
    p[ok] <- (alt_ct[ok] + pseudocount) / (tot_ct[ok] + 2 * pseudocount)
    
    freq[, g] <- p
  }
  
  list(freq = freq, groups = groups, ref_ids = ref_ids)
}

assign_one_K <- function(
    snp_mat,
    freq,
    groups,
    sample_id,
    panel_id,
    group_col,
    K,
    epsilon,
    min_snps_query,
    high_thr,
    mod_thr,
    min_gap
) {
  qvec <- snp_mat[, sample_id]
  
  ok_anyp <- rowSums(!is.na(freq)) > 0
  n_used <- sum(!is.na(qvec) & ok_anyp)
  
  if (is.na(n_used) || n_used < min_snps_query) {
    return(data.table(
      panel_id = panel_id,
      group_col = group_col,
      K = as.integer(K),
      sample_id = sample_id,
      top_group = NA_character_,
      top_posterior = NA_real_,
      second_group = NA_character_,
      second_posterior = NA_real_,
      posterior_gap = NA_real_,
      n_snps_used = as.integer(n_used),
      confidence = "Low",
      status = "too_few_snps",
      note = paste0("Need >= ", min_snps_query, "; have ", n_used)
    ))
  }
  
  lls <- vapply(
    groups,
    function(g) loglik_one_group(qvec, freq[, g], epsilon = epsilon),
    numeric(1)
  )
  
  post <- post_from_ll(lls)
  
  if (all(is.na(post))) {
    return(data.table(
      panel_id = panel_id,
      group_col = group_col,
      K = as.integer(K),
      sample_id = sample_id,
      top_group = NA_character_,
      top_posterior = NA_real_,
      second_group = NA_character_,
      second_posterior = NA_real_,
      posterior_gap = NA_real_,
      n_snps_used = as.integer(n_used),
      confidence = "Low",
      status = "no_likelihood",
      note = "All likelihoods NA"
    ))
  }
  
  ord <- order(post, decreasing = TRUE)
  
  top1 <- ord[1]
  top2 <- if (length(ord) >= 2L) ord[2] else NA_integer_
  
  tg <- names(post)[top1]
  tp <- unname(post[top1])
  
  sg <- if (!is.na(top2)) names(post)[top2] else NA_character_
  sp <- if (!is.na(top2)) unname(post[top2]) else NA_real_
  
  gap <- if (!is.na(sp)) tp - sp else NA_real_
  
  conf <- conf_label(tp, high = high_thr, mod = mod_thr)
  
  st <- "ok"
  nt <- ""
  
  if (!is.na(gap) && gap < min_gap) {
    st <- "low_gap"
    nt <- paste0("posterior_gap<", min_gap)
  }
  
  data.table(
    panel_id = panel_id,
    group_col = group_col,
    K = as.integer(K),
    sample_id = sample_id,
    top_group = tg,
    top_posterior = tp,
    second_group = sg,
    second_posterior = sp,
    posterior_gap = gap,
    n_snps_used = as.integer(n_used),
    confidence = conf,
    status = st,
    note = nt
  )
}

summarise_multiK_one_sample <- function(raw_dt, tail_fraction) {
  dt_all <- raw_dt
  dt_ok <- raw_dt[status %in% c("ok", "low_gap") & !is.na(top_group)]
  
  nK_total <- length(unique(dt_all$K))
  nK <- length(unique(dt_ok$K))
  maxK_used <- if (nK > 0L) max(dt_ok$K) else NA_real_
  
  final_group <- if (!is.na(maxK_used)) {
    dt_ok[K == maxK_used][1, top_group]
  } else {
    NA_character_
  }
  
  agree <- if (!is.na(final_group) && nK > 0L) {
    mean(dt_ok$top_group == final_group)
  } else {
    NA_real_
  }
  
  stable_from_K <- stable_from_threshold(dt_ok[order(K)], final_group)
  
  tail_agree <- tail_agreement_with_final(
    dt = dt_ok,
    final_group = final_group,
    tail_fraction = tail_fraction
  )
  
  maxrow <- if (!is.na(maxK_used)) {
    dt_ok[K == maxK_used][1]
  } else {
    NULL
  }
  
  data.table(
    final_group = final_group,
    stable_from_K = as.numeric(stable_from_K),
    agreement_frac = as.numeric(agree),
    tail_agreement_frac = as.numeric(tail_agree),
    tail_fraction = as.numeric(tail_fraction),
    n_K_total_tested = as.integer(nK_total),
    n_K_usable_for_stability = as.integer(nK),
    n_K_available = as.integer(nK),
    maxK_used = as.numeric(maxK_used),
    maxK_top_posterior = if (!is.null(maxrow)) as.numeric(maxrow$top_posterior) else NA_real_,
    maxK_second_posterior = if (!is.null(maxrow)) as.numeric(maxrow$second_posterior) else NA_real_,
    maxK_gap = if (!is.null(maxrow)) as.numeric(maxrow$posterior_gap) else NA_real_,
    maxK_confidence = if (!is.null(maxrow)) as.character(maxrow$confidence) else NA_character_,
    maxK_status = if (!is.null(maxrow)) as.character(maxrow$status) else NA_character_,
    maxK_n_snps_used = if (!is.null(maxrow)) as.integer(maxrow$n_snps_used) else NA_integer_,
    status = if (nK > 0L) "ok" else "no_valid_K"
  )
}

final_call <- function(
    stab_row,
    high_thr,
    mod_thr,
    min_gap,
    min_agreement,
    min_K,
    min_tail_agreement
) {
  if (is.null(stab_row) || nrow(stab_row) == 0L) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_result"))
  }
  
  s <- stab_row[1]
  
  if (!("final_group" %in% names(s)) || is.na(s$final_group)) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_final_group"))
  }
  
  if (is.na(s$maxK_n_snps_used) || s$maxK_n_snps_used <= 0) {
    return(list(call = NA_character_, label = "Uncertain", reason = "no_snps_used"))
  }
  
  p <- s$maxK_top_posterior
  gap <- s$maxK_gap
  
  if (is.na(p)) {
    return(list(call = NA_character_, label = "Uncertain", reason = "posterior_NA"))
  }
  
  conf <- conf_label(p, high = high_thr, mod = mod_thr)
  
  if (!(conf %in% c("High", "Moderate"))) {
    return(list(call = s$final_group, label = "Uncertain", reason = "low_posterior"))
  }
  
  if (!is.na(gap) && gap < min_gap) {
    return(list(call = s$final_group, label = "Uncertain", reason = paste0("gap<", min_gap)))
  }
  
  if (is.na(s$n_K_available) || s$n_K_available < min_K) {
    return(list(call = s$final_group, label = "Uncertain", reason = "too_few_usable_K"))
  }
  
  if (is.na(s$agreement_frac) || s$agreement_frac < min_agreement) {
    return(list(call = s$final_group, label = "Uncertain", reason = "low_global_K_agreement"))
  }
  
  if ("tail_agreement_frac" %in% names(s)) {
    if (is.na(s$tail_agreement_frac) || s$tail_agreement_frac < min_tail_agreement) {
      return(list(call = s$final_group, label = "Uncertain", reason = "unstable_large_K_tail"))
    }
  }
  
  list(call = s$final_group, label = conf, reason = "ok")
}

# =============================================================================
# Load metadata
# =============================================================================

meta_path <- file.path(opt$qc_dir, "metadata_clean.tsv")

if (!file.exists(meta_path)) {
  stop("Missing metadata_clean.tsv: ", meta_path, " (run Step 1-2 first)")
}

meta <- fread(meta_path)
setnames(meta, tolower(names(meta)))

required_meta <- c("sample_id", "is_reference", "macroregion_3", "subregion")
missing_meta <- setdiff(required_meta, names(meta))

if (length(missing_meta) > 0L) {
  stop(
    "metadata_clean.tsv missing required column(s): ",
    paste(missing_meta, collapse = ", "),
    "\nFound: ",
    paste(names(meta), collapse = ", ")
  )
}

for (nm in names(meta)) {
  if (!inherits(meta[[nm]], c("numeric", "integer", "logical", "Date", "POSIXct"))) {
    meta[[nm]] <- clean_text(meta[[nm]])
  }
}

meta[, sample_id := clean_text(sample_id)]
meta[, is_reference := as_clean_logical(is_reference)]
meta[, macroregion_3 := clean_text(macroregion_3)]
meta[, subregion := clean_text(subregion)]

if ("country" %in% names(meta)) {
  meta[, country := clean_text(country)]
}

if ("population" %in% names(meta)) {
  meta[, population := clean_text(population)]
}

# =============================================================================
# K-grid setup
# =============================================================================

baseK <- parse_integer_vector(opt$base_K)
extraK <- parse_integer_vector(opt$extra_K)

auto_extend_K <- parse_bool1(
  opt$auto_extend_K,
  default = TRUE,
  argname = "--auto_extend_K"
)

include_all_snps_K <- parse_bool1(
  opt$include_all_snps_K,
  default = TRUE,
  argname = "--include_all_snps_K"
)

# =============================================================================
# Run panel LOO
# =============================================================================

run_one_panel_loo <- function(cfg) {
  panel_id <- cfg$panel_id
  group_col <- cfg$group_col
  prefix <- cfg$prefix
  min_snps_query <- cfg$min_snps_query
  min_gap <- cfg$min_gap
  
  message("\n============================================================")
  message("Running Step 4b LOO panel: ", panel_id)
  message("Group column: ", group_col)
  message("Output prefix: ", prefix)
  message("============================================================")
  
  panel_dir <- file.path(opt$step3_dir, "panels", panel_id)
  snp_matrix_path <- file.path(panel_dir, "snp_matrix.rds")
  snp_map_path <- file.path(panel_dir, "snp_map.tsv")
  
  if (!file.exists(snp_matrix_path)) {
    warning("Skipping panel because snp_matrix.rds is missing: ", snp_matrix_path)
    return(NULL)
  }
  
  if (!file.exists(snp_map_path)) {
    warning("Skipping panel because snp_map.tsv is missing: ", snp_map_path)
    return(NULL)
  }
  
  X0 <- readRDS(snp_matrix_path)
  
  if (!is.matrix(X0)) {
    X0 <- as.matrix(X0)
  }
  
  storage.mode(X0) <- "integer"
  
  if (is.null(rownames(X0)) || is.null(colnames(X0))) {
    stop("snp_matrix.rds must have rownames=snp_id and colnames=sample_id for ", panel_id)
  }
  
  snp_map <- fread(snp_map_path)
  
  if (!("snp_id" %in% names(snp_map))) {
    stop("snp_map.tsv missing snp_id column for ", panel_id)
  }
  
  if ("global_rank" %in% names(snp_map)) {
    snp_map <- snp_map[order(global_rank)]
  } else if ("score" %in% names(snp_map)) {
    snp_map <- snp_map[order(-score)]
  }
  
  snp_map <- snp_map[snp_id %in% rownames(X0)]
  
  if (nrow(snp_map) == 0L) {
    warning("Skipping panel because no SNPs in snp_map match snp_matrix rownames: ", panel_id)
    return(NULL)
  }
  
  ref_meta <- meta[is_reference %in% TRUE]
  ref_meta <- cfg$filter_fun(ref_meta)
  
  if (!(group_col %in% names(ref_meta))) {
    stop("Metadata missing group column ", group_col, " for panel ", panel_id)
  }
  
  ref_meta <- ref_meta[!is.na(get(group_col)) & clean_text(get(group_col)) != ""]
  ref_meta <- ref_meta[sample_id %in% colnames(X0)]
  
  group_counts <- ref_meta[, .N, by = group_col]
  setnames(group_counts, group_col, "group")
  
  keep_groups <- group_counts[N >= opt$min_ref_per_group, group]
  ref_meta <- ref_meta[get(group_col) %in% keep_groups]
  
  group_counts_final <- ref_meta[, .N, by = group_col]
  setnames(group_counts_final, group_col, "group")
  
  if (nrow(group_counts_final) < 2L) {
    warning("Skipping panel because fewer than two groups remain after filtering: ", panel_id)
    return(NULL)
  }
  
  ref_ids <- ref_meta$sample_id
  
  keep_snps <- rowMeans(!is.na(X0[, ref_ids, drop = FALSE])) >= opt$keep_snps_frac
  X <- X0[keep_snps, , drop = FALSE]
  
  snp_map <- snp_map[snp_id %in% rownames(X)]
  
  if (nrow(snp_map) == 0L) {
    warning("Skipping panel because no SNPs remain after completeness filter: ", panel_id)
    return(NULL)
  }
  
  snp_order <- snp_map$snp_id
  n_snps <- length(snp_order)
  
  Ks <- make_Ks(
    n_snps = n_snps,
    baseK = baseK,
    extraK = extraK,
    auto_extend = auto_extend_K,
    include_all = include_all_snps_K,
    max_points = opt$max_K_points
  )
  
  message("Reference samples: ", nrow(ref_meta))
  message("Groups: ", paste(group_counts_final$group, group_counts_final$N, sep = "=", collapse = ", "))
  message("SNPs retained: ", n_snps)
  message("K grid: ", paste(Ks, collapse = ", "))
  
  if (length(Ks) == 0L) {
    warning("Skipping panel because K grid is empty: ", panel_id)
    return(NULL)
  }
  
  loo_raw_list <- vector("list", nrow(ref_meta))
  loo_stab_list <- vector("list", nrow(ref_meta))
  
  for (i in seq_len(nrow(ref_meta))) {
    sid <- ref_meta$sample_id[i]
    true_g <- ref_meta[[group_col]][i]
    
    train_meta <- ref_meta[sample_id != sid]
    
    rf <- build_ref_freqs(
      snp_mat = X[snp_order, , drop = FALSE],
      meta_ref = train_meta,
      group_col = group_col,
      pseudocount = opt$pseudocount
    )
    
    if (is.null(rf$freq) || length(rf$groups) < 2L) {
      raw_i <- data.table(
        panel_id = panel_id,
        group_col = group_col,
        K = Ks,
        sample_id = sid,
        top_group = NA_character_,
        top_posterior = NA_real_,
        second_group = NA_character_,
        second_posterior = NA_real_,
        posterior_gap = NA_real_,
        n_snps_used = NA_integer_,
        confidence = "Low",
        status = "no_freqs",
        note = "Fewer than two train groups after leave-one-out"
      )
      
      stab_i <- summarise_multiK_one_sample(raw_i, tail_fraction = opt$tail_fraction)
      
    } else {
      raw_by_K <- vector("list", length(Ks))
      
      for (k_i in seq_along(Ks)) {
        K <- Ks[k_i]
        snp_ids_K <- snp_order[seq_len(min(K, length(snp_order)))]
        
        raw_by_K[[k_i]] <- assign_one_K(
          snp_mat = X[snp_ids_K, , drop = FALSE],
          freq = rf$freq[snp_ids_K, , drop = FALSE],
          groups = rf$groups,
          sample_id = sid,
          panel_id = panel_id,
          group_col = group_col,
          K = K,
          epsilon = opt$epsilon,
          min_snps_query = min_snps_query,
          high_thr = opt$high_posterior,
          mod_thr = opt$moderate_posterior,
          min_gap = min_gap
        )
      }
      
      raw_i <- rbindlist(raw_by_K, fill = TRUE)
      stab_i <- summarise_multiK_one_sample(raw_i, tail_fraction = opt$tail_fraction)
    }
    
    fc <- final_call(
      stab_i,
      high_thr = opt$high_posterior,
      mod_thr = opt$moderate_posterior,
      min_gap = min_gap,
      min_agreement = opt$min_agreement,
      min_K = opt$min_K_available,
      min_tail_agreement = opt$min_tail_agreement
    )
    
    # Critical: explicitly attach these columns to avoid data.table resolving
    # panel_id/group_col from the parent function environment.
    stab_i[, `:=`(
      panel_id = panel_id,
      group_col = group_col,
      sample_id = sid,
      true_group = true_g,
      predicted_group = fc$call,
      final_confidence = fc$label,
      final_reason = fc$reason,
      correct = !is.na(fc$call) && fc$call == true_g
    )]
    
    raw_i[, `:=`(
      panel_id = panel_id,
      group_col = group_col,
      sample_id = sid,
      true_group = true_g
    )]
    
    loo_raw_list[[i]] <- raw_i
    loo_stab_list[[i]] <- stab_i
  }
  
  loo_raw <- rbindlist(loo_raw_list, fill = TRUE)
  loo_pred <- rbindlist(loo_stab_list, fill = TRUE)
  
  setcolorder(
    loo_pred,
    intersect(
      c(
        "panel_id",
        "group_col",
        "sample_id",
        "true_group",
        "predicted_group",
        "final_confidence",
        "final_reason",
        "correct",
        "final_group",
        "maxK_used",
        "maxK_top_posterior",
        "maxK_second_posterior",
        "maxK_gap",
        "maxK_confidence",
        "maxK_status",
        "maxK_n_snps_used",
        "agreement_frac",
        "tail_agreement_frac",
        "tail_fraction",
        "stable_from_K",
        "n_K_total_tested",
        "n_K_usable_for_stability",
        "n_K_available",
        "status"
      ),
      names(loo_pred)
    )
  )
  
  loo_summary <- loo_pred[
    ,
    .(
      n = .N,
      attempted = sum(!is.na(predicted_group)),
      correct = sum(correct, na.rm = TRUE),
      accuracy_all = mean(correct, na.rm = TRUE),
      accuracy_attempted = ifelse(
        sum(!is.na(predicted_group)) > 0,
        mean(correct[!is.na(predicted_group)], na.rm = TRUE),
        NA_real_
      ),
      high_n = sum(final_confidence == "High", na.rm = TRUE),
      moderate_n = sum(final_confidence == "Moderate", na.rm = TRUE),
      uncertain_n = sum(final_confidence == "Uncertain", na.rm = TRUE),
      high_frac = mean(final_confidence == "High", na.rm = TRUE),
      high_or_moderate_frac = mean(final_confidence %in% c("High", "Moderate"), na.rm = TRUE),
      uncertainty_rate = mean(final_confidence == "Uncertain", na.rm = TRUE),
      mean_maxK_posterior = mean(maxK_top_posterior, na.rm = TRUE),
      mean_maxK_gap = mean(maxK_gap, na.rm = TRUE),
      mean_agreement = mean(agreement_frac, na.rm = TRUE),
      mean_tail_agreement = mean(tail_agreement_frac, na.rm = TRUE),
      too_few_usable_K_rate = mean(final_reason == "too_few_usable_K", na.rm = TRUE),
      too_few_snps_raw_rate = mean(loo_raw$status == "too_few_snps", na.rm = TRUE)
    ),
    by = .(panel_id, group_col)
  ]
  
  loo_perclass <- loo_pred[
    ,
    .(
      n = .N,
      attempted = sum(!is.na(predicted_group)),
      correct = sum(correct, na.rm = TRUE),
      accuracy_all = mean(correct, na.rm = TRUE),
      accuracy_attempted = ifelse(
        sum(!is.na(predicted_group)) > 0,
        mean(correct[!is.na(predicted_group)], na.rm = TRUE),
        NA_real_
      ),
      high_or_moderate_frac = mean(final_confidence %in% c("High", "Moderate"), na.rm = TRUE),
      uncertainty_rate = mean(final_confidence == "Uncertain", na.rm = TRUE)
    ),
    by = .(panel_id, true_group)
  ][order(true_group)]
  
  confusion <- loo_pred[
    !is.na(predicted_group),
    .N,
    by = .(panel_id, true_group, predicted_group)
  ]
  
  if (nrow(confusion) > 0L) {
    confusion[, prop_within_true := N / sum(N), by = .(panel_id, true_group)]
  }
  
  reason_counts <- loo_pred[
    ,
    .N,
    by = .(panel_id, final_reason)
  ][order(panel_id, -N)]
  
  K_grid_panel <- data.table(
    panel_id = panel_id,
    n_snps_available = n_snps,
    K = Ks
  )
  
  pred_file <- file.path(opt$out_dir, paste0(prefix, "_loo_predictions.tsv"))
  raw_file <- file.path(opt$out_dir, paste0(prefix, "_loo_raw.tsv.gz"))
  summary_file <- file.path(opt$out_dir, paste0(prefix, "_loo_summary.tsv"))
  perclass_file <- file.path(opt$out_dir, paste0(prefix, "_loo_perclass.tsv"))
  conf_file <- file.path(opt$out_dir, paste0(prefix, "_loo_confusion.tsv"))
  reason_file <- file.path(opt$out_dir, paste0(prefix, "_loo_reason_counts.tsv"))
  
  fwrite(loo_pred, pred_file, sep = "\t")
  saveRDS(loo_pred, sub("\\.tsv$", ".rds", pred_file))
  
  fwrite(loo_raw, raw_file, sep = "\t")
  
  fwrite(loo_summary, summary_file, sep = "\t")
  saveRDS(loo_summary, sub("\\.tsv$", ".rds", summary_file))
  
  fwrite(loo_perclass, perclass_file, sep = "\t")
  saveRDS(loo_perclass, sub("\\.tsv$", ".rds", perclass_file))
  
  fwrite(confusion, conf_file, sep = "\t")
  saveRDS(confusion, sub("\\.tsv$", ".rds", conf_file))
  
  fwrite(reason_counts, reason_file, sep = "\t")
  
  # ---------------------------------------------------------------------------
  # Base R A4 validation figure
  # ---------------------------------------------------------------------------
  
  fig_file <- file.path(opt$out_dir, paste0("FIG_", prefix, "_validation_A4.pdf"))
  
  grDevices::pdf(fig_file, width = 8.27, height = 11.69)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  
  par(mfrow = c(3, 1), mar = c(5, 5, 4, 2))
  
  if (nrow(confusion) > 0L) {
    conf_mat <- dcast(confusion, true_group ~ predicted_group, value.var = "N", fill = 0)
    true_names <- conf_mat$true_group
    mat <- as.matrix(conf_mat[, -"true_group"])
    
    rownames(mat) <- true_names
    
    image(
      x = seq_len(ncol(mat)),
      y = seq_len(nrow(mat)),
      z = t(mat[nrow(mat):1, , drop = FALSE]),
      axes = FALSE,
      main = paste0(panel_id, "\nLOO confusion matrix"),
      xlab = "Predicted group",
      ylab = "True group"
    )
    
    axis(1, at = seq_len(ncol(mat)), labels = colnames(mat), las = 2)
    axis(2, at = seq_len(nrow(mat)), labels = rev(rownames(mat)), las = 2)
    
    for (rr in seq_len(nrow(mat))) {
      for (cc in seq_len(ncol(mat))) {
        text(
          x = cc,
          y = nrow(mat) - rr + 1,
          labels = mat[rr, cc],
          cex = 0.9
        )
      }
    }
  } else {
    plot.new()
    title("No attempted predictions")
  }
  
  metrics <- c(
    accuracy_all = loo_summary$accuracy_all[1],
    accuracy_attempted = loo_summary$accuracy_attempted[1],
    uncertainty_rate = loo_summary$uncertainty_rate[1]
  )
  
  barplot(
    metrics,
    ylim = c(0, 1),
    las = 2,
    ylab = "Proportion",
    main = "Accuracy and uncertainty"
  )
  abline(h = seq(0, 1, 0.25), lty = 3, col = "grey70")
  
  if (nrow(loo_perclass) > 0L) {
    barplot(
      loo_perclass$accuracy_all,
      names.arg = loo_perclass$true_group,
      ylim = c(0, 1),
      las = 2,
      ylab = "Accuracy",
      main = "Accuracy by true group"
    )
    abline(h = loo_summary$accuracy_all[1], lty = 2)
  } else {
    plot.new()
    title("No per-class summary")
  }
  
  dev.off()
  
  panel_index <- data.table(
    panel_id = panel_id,
    group_col = group_col,
    prefix = prefix,
    panel_dir = panel_dir,
    n_reference_samples = nrow(ref_meta),
    n_reference_groups = nrow(group_counts_final),
    n_snps_available = n_snps,
    n_K_values = length(Ks),
    min_snps_query = min_snps_query,
    min_gap = min_gap,
    accuracy_all = loo_summary$accuracy_all[1],
    accuracy_attempted = loo_summary$accuracy_attempted[1],
    uncertainty_rate = loo_summary$uncertainty_rate[1],
    high_or_moderate_frac = loo_summary$high_or_moderate_frac[1],
    predictions_file = pred_file,
    raw_file = raw_file,
    summary_file = summary_file,
    perclass_file = perclass_file,
    confusion_file = conf_file,
    figure_file = fig_file
  )
  
  list(
    index = panel_index,
    predictions = loo_pred,
    raw = loo_raw,
    summary = loo_summary,
    perclass = loo_perclass,
    confusion = confusion,
    reasons = reason_counts,
    K_grid = K_grid_panel,
    group_counts = group_counts_final
  )
}

# =============================================================================
# Execute selected panels
# =============================================================================

results <- list()

for (pid in panels_requested) {
  cfg <- panel_config(pid)
  res <- run_one_panel_loo(cfg)
  
  if (!is.null(res)) {
    results[[pid]] <- res
  }
}

if (length(results) == 0L) {
  stop("No Step 4b panels completed successfully.")
}

# =============================================================================
# Global outputs
# =============================================================================

panels_index <- rbindlist(lapply(results, `[[`, "index"), fill = TRUE)
all_predictions <- rbindlist(lapply(results, `[[`, "predictions"), fill = TRUE)
all_summary <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
all_perclass <- rbindlist(lapply(results, `[[`, "perclass"), fill = TRUE)
all_confusion <- rbindlist(lapply(results, `[[`, "confusion"), fill = TRUE)
all_reasons <- rbindlist(lapply(results, `[[`, "reasons"), fill = TRUE)
all_K_grid <- rbindlist(lapply(results, `[[`, "K_grid"), fill = TRUE)
all_raw <- rbindlist(lapply(results, `[[`, "raw"), fill = TRUE)

fwrite(panels_index, file.path(opt$out_dir, "step4b_loo_panels_index.tsv"), sep = "\t")
saveRDS(panels_index, file.path(opt$out_dir, "step4b_loo_panels_index.rds"))

fwrite(all_predictions, file.path(opt$out_dir, "step4b_loo_predictions_all_panels.tsv"), sep = "\t")
saveRDS(all_predictions, file.path(opt$out_dir, "step4b_loo_predictions_all_panels.rds"))

fwrite(all_summary, file.path(opt$out_dir, "step4b_loo_summary_all_panels.tsv"), sep = "\t")
saveRDS(all_summary, file.path(opt$out_dir, "step4b_loo_summary_all_panels.rds"))

fwrite(all_perclass, file.path(opt$out_dir, "step4b_loo_perclass_all_panels.tsv"), sep = "\t")
saveRDS(all_perclass, file.path(opt$out_dir, "step4b_loo_perclass_all_panels.rds"))

fwrite(all_confusion, file.path(opt$out_dir, "step4b_loo_confusion_all_panels.tsv"), sep = "\t")
saveRDS(all_confusion, file.path(opt$out_dir, "step4b_loo_confusion_all_panels.rds"))

fwrite(all_reasons, file.path(opt$out_dir, "step4b_loo_reason_counts_all_panels.tsv"), sep = "\t")

fwrite(all_K_grid, file.path(opt$out_dir, "step4b_K_grid_used.tsv"), sep = "\t")
saveRDS(all_K_grid, file.path(opt$out_dir, "step4b_K_grid_used.rds"))

fwrite(
  all_raw,
  file.path(opt$out_dir, "step4b_loo_raw_all_panels.tsv.gz"),
  sep = "\t"
)

params <- data.table(
  parameter = names(opt),
  value = vapply(opt, function(x) paste(x, collapse = ","), character(1))
)

params <- rbind(
  params,
  data.table(parameter = "baseK_parsed", value = paste(baseK, collapse = ",")),
  data.table(parameter = "extraK_parsed", value = paste(extraK, collapse = ",")),
  data.table(parameter = "auto_extend_K_parsed", value = as.character(auto_extend_K)),
  data.table(parameter = "include_all_snps_K_parsed", value = as.character(include_all_snps_K))
)

fwrite(params, file.path(opt$out_dir, "step4b_params.tsv"), sep = "\t")
saveRDS(params, file.path(opt$out_dir, "step4b_params.rds"))

# =============================================================================
# Excel workbook
# =============================================================================

wb <- createWorkbook()

sn <- safe_sheet(wb, "PANELS_INDEX")
addWorksheet(wb, sn)
writeDataTable(wb, sn, panels_index)

sn <- safe_sheet(wb, "SUMMARY_ALL")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_summary)

sn <- safe_sheet(wb, "PREDICTIONS_ALL")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_predictions)

sn <- safe_sheet(wb, "PERCLASS_ALL")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_perclass)

sn <- safe_sheet(wb, "CONFUSION_ALL")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_confusion)

sn <- safe_sheet(wb, "REASONS_ALL")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_reasons)

sn <- safe_sheet(wb, "K_GRID_USED")
addWorksheet(wb, sn)
writeDataTable(wb, sn, all_K_grid)

sn <- safe_sheet(wb, "PARAMS")
addWorksheet(wb, sn)
writeDataTable(wb, sn, params)

xlsx_out <- file.path(opt$out_dir, "step4b_loo_validation_summary.xlsx")
saveWorkbook(wb, xlsx_out, overwrite = TRUE)

# =============================================================================
# Run info
# =============================================================================

run_info <- list(
  step = "step4b_loo_validation",
  repo_root = repo_root,
  qc_dir = opt$qc_dir,
  step3_dir = opt$step3_dir,
  out_dir = opt$out_dir,
  panels_requested = panels_requested,
  panels_index = panels_index,
  baseK = baseK,
  extraK = extraK,
  auto_extend_K = auto_extend_K,
  include_all_snps_K = include_all_snps_K,
  K_grid = all_K_grid,
  parameters = opt,
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(run_info, file.path(opt$out_dir, "step4b_run_info.rds"))

message("\nDone.")
message("Step 4b panel index: ", file.path(opt$out_dir, "step4b_loo_panels_index.tsv"))
message("Step 4b summary:     ", file.path(opt$out_dir, "step4b_loo_summary_all_panels.tsv"))
message("Step 4b predictions: ", file.path(opt$out_dir, "step4b_loo_predictions_all_panels.tsv"))
message("Step 4b workbook:    ", xlsx_out)