#!/usr/bin/env Rscript

# =============================================================================
# PESTFLY — Step 0: Convert OMA/read2tree interleaved PHYLIP MSAs to FASTA
# =============================================================================
#
# Purpose
# -------
# Convert per-ortholog OMA/read2tree PHYLIP alignments into FASTA alignments.
#
# Input
# -----
#   data/000_input_data/phy/OG*.fa
#
# Output
# ------
#   results/00_fasta/OG*.fasta
#   results/00_fasta/convert_summary.tsv
#   results/00_fasta/convert_summary.rds
#   results/00_fasta/convert_errors.tsv      if errors occur
#   results/00_fasta/convert_errors.rds      if errors occur
#   results/00_fasta/step0_run_info.rds
#
# Notes
# -----
# This parser is designed for OMA/read2tree-style interleaved PHYLIP:
#
#   ntax nsites
#   label_1   sequence_chunk_1
#   label_2   sequence_chunk_1
#   ...
#   label_N   sequence_chunk_1
#
#             continuation_chunk_for_label_1
#             continuation_chunk_for_label_2
#             ...
#
# It reconstructs full sequences, checks them against the header n_sites,
# and writes one FASTA file per OG.
#
# Run
# ---
#   Rscript steps/step0_convert_phylip_to_fasta/run.R
#
# Test
# ----
#   Rscript steps/step0_convert_phylip_to_fasta/run.R --max_files 20 --cores 4
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
})

# =============================================================================
# Repo-root detection
# =============================================================================

find_repo_root <- function(max_up = 8L) {
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
    "--in_dir",
    type = "character",
    default = "data/000_input_data/phy",
    help = "Input directory with OMA/read2tree PHYLIP MSA files [default %default]"
  ),
  make_option(
    "--out_dir",
    type = "character",
    default = "results/00_fasta",
    help = "Output directory for converted FASTA files [default %default]"
  ),
  make_option(
    "--file_glob",
    type = "character",
    default = "OG*.fa",
    help = "Input PHYLIP file pattern [default %default]"
  ),
  make_option(
    "--max_files",
    type = "integer",
    default = 0,
    help = "For testing: process first N files only; 0 = all [default %default]"
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 0,
    help = "Parallel workers. 0 = all detected cores minus one [default %default]"
  ),
  make_option(
    "--wrap_width",
    type = "integer",
    default = 80,
    help = "FASTA line width [default %default]"
  ),
  make_option(
    "--strict_length",
    action = "store_true",
    default = FALSE,
    help = "If TRUE, stop when reconstructed sequence length != header n_sites [default %default]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$in_dir <- resolve_path(opt$in_dir)
opt$out_dir <- resolve_path(opt$out_dir)

if (opt$cores <= 0) {
  opt$cores <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

message("Repo root:  ", repo_root)
message("Input dir:  ", opt$in_dir)
message("Output dir: ", opt$out_dir)
message("File glob:  ", opt$file_glob)
message("Cores:      ", opt$cores)

# =============================================================================
# PHYLIP parsing helpers
# =============================================================================

clean_seq <- function(x) {
  toupper(gsub("\\s+", "", x))
}

is_header_line <- function(x) {
  grepl("^\\s*\\d+\\s+\\d+\\s*$", x)
}

parse_header <- function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  
  if (length(parts) < 2) {
    stop("Malformed PHYLIP header: ", line)
  }
  
  ntax <- suppressWarnings(as.integer(parts[1]))
  nsites <- suppressWarnings(as.integer(parts[2]))
  
  if (is.na(ntax) || is.na(nsites) || ntax <= 0 || nsites <= 0) {
    stop("Invalid PHYLIP header: ", line)
  }
  
  list(ntax = ntax, nsites = nsites)
}

split_labeled_line <- function(line) {
  m <- regexec("^\\s*(\\S+)\\s+(.+?)\\s*$", line)
  hit <- regmatches(line, m)[[1]]
  
  if (length(hit) < 3) {
    stop("Could not parse labelled PHYLIP line: ", line)
  }
  
  list(
    label = hit[2],
    seq = clean_seq(hit[3])
  )
}

read_interleaved_phylip <- function(path, strict_length = FALSE) {
  x <- readLines(path, warn = FALSE)
  x <- sub("\r$", "", x)
  
  nonempty <- which(nzchar(trimws(x)))
  if (length(nonempty) == 0L) {
    stop("Empty file")
  }
  
  header_index <- nonempty[1]
  header_line <- trimws(sub("^\ufeff", "", x[header_index]))
  
  if (!is_header_line(header_line)) {
    stop(
      "First non-empty line is not a PHYLIP header in ",
      basename(path),
      ": ",
      header_line
    )
  }
  
  h <- parse_header(header_line)
  ntax <- h$ntax
  nsites <- h$nsites
  
  body <- x[(header_index + 1L):length(x)]
  body <- body[nzchar(trimws(body))]
  
  if (length(body) < ntax) {
    stop(
      basename(path),
      ": fewer body lines than ntax. ntax = ",
      ntax,
      "; body lines = ",
      length(body)
    )
  }
  
  first_block <- body[seq_len(ntax)]
  rest <- body[-seq_len(ntax)]
  
  labels <- character(ntax)
  seqs <- character(ntax)
  
  for (i in seq_len(ntax)) {
    parsed <- split_labeled_line(first_block[i])
    labels[i] <- parsed$label
    seqs[i] <- parsed$seq
  }
  
  # Continuation blocks are usually sequence-only, in the same taxon order.
  # Some PHYLIP writers repeat labels in continuation blocks; this handles both.
  if (length(rest) > 0L) {
    pos <- 1L
    
    while (pos <= length(rest) && any(nchar(seqs) < nsites)) {
      for (k in seq_len(ntax)) {
        if (pos > length(rest)) break
        
        line <- rest[pos]
        pos <- pos + 1L
        
        tokens <- strsplit(trimws(line), "\\s+")[[1]]
        
        if (length(tokens) >= 2L && tokens[1] %in% labels) {
          lab <- tokens[1]
          seq_part <- clean_seq(paste(tokens[-1], collapse = ""))
          idx <- match(lab, labels)
          seqs[idx] <- paste0(seqs[idx], seq_part)
        } else {
          seq_part <- clean_seq(line)
          seqs[k] <- paste0(seqs[k], seq_part)
        }
      }
    }
  }
  
  lens_raw <- nchar(seqs)
  had_length_problem <- any(lens_raw != nsites)
  
  if (had_length_problem) {
    msg <- paste0(
      basename(path),
      ": reconstructed sequence length range ",
      min(lens_raw),
      "-",
      max(lens_raw),
      "; header n_sites = ",
      nsites
    )
    
    if (isTRUE(strict_length)) {
      stop(msg)
    }
    
    warning(msg, "; padding/trimming to n_sites")
    
    too_short <- lens_raw < nsites
    if (any(too_short)) {
      seqs[too_short] <- paste0(
        seqs[too_short],
        strrep("-", nsites - lens_raw[too_short])
      )
    }
    
    seqs <- substr(seqs, 1L, nsites)
  }
  
  list(
    taxa = labels,
    seqs = seqs,
    n_taxa = ntax,
    n_sites = nsites,
    min_len_raw = min(lens_raw),
    max_len_raw = max(lens_raw),
    had_length_problem = had_length_problem
  )
}

# =============================================================================
# FASTA writing
# =============================================================================

write_fasta <- function(taxa, seqs, out_file, width = 80L) {
  con <- file(out_file, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_along(taxa)) {
    writeLines(paste0(">", taxa[i]), con)
    
    s <- seqs[i]
    starts <- seq.int(1L, nchar(s), by = width)
    
    for (st in starts) {
      writeLines(substr(s, st, min(st + width - 1L, nchar(s))), con)
    }
  }
}

# =============================================================================
# Convert one file
# =============================================================================

convert_one <- function(f, out_dir, wrap_width, strict_length) {
  og_id <- tools::file_path_sans_ext(basename(f))
  out_f <- file.path(out_dir, paste0(og_id, ".fasta"))
  
  res <- tryCatch(
    read_interleaved_phylip(f, strict_length = strict_length),
    error = function(e) e
  )
  
  if (inherits(res, "error")) {
    return(list(
      summary = data.table(
        og = og_id,
        input_file = f,
        output_file = out_f,
        n_taxa = NA_integer_,
        n_sites = NA_integer_,
        min_len_raw = NA_integer_,
        max_len_raw = NA_integer_,
        had_length_problem = NA,
        status = "ERROR"
      ),
      error = data.table(
        og = og_id,
        file = f,
        error = conditionMessage(res)
      )
    ))
  }
  
  write_fasta(
    taxa = res$taxa,
    seqs = res$seqs,
    out_file = out_f,
    width = wrap_width
  )
  
  list(
    summary = data.table(
      og = og_id,
      input_file = f,
      output_file = out_f,
      n_taxa = res$n_taxa,
      n_sites = res$n_sites,
      min_len_raw = res$min_len_raw,
      max_len_raw = res$max_len_raw,
      had_length_problem = res$had_length_problem,
      status = "OK"
    ),
    error = NULL
  )
}

# =============================================================================
# Locate files
# =============================================================================

files <- sort(Sys.glob(file.path(opt$in_dir, opt$file_glob)))

if (length(files) == 0L) {
  stop("No input PHYLIP files found in: ", opt$in_dir)
}

if (opt$max_files > 0L) {
  files <- head(files, opt$max_files)
}

message("Found ", length(files), " PHYLIP file(s).")

# =============================================================================
# Parallel conversion
# =============================================================================

process_one_file <- function(f) {
  convert_one(
    f = f,
    out_dir = opt$out_dir,
    wrap_width = opt$wrap_width,
    strict_length = opt$strict_length
  )
}

message("Converting in parallel...")

if (.Platform$OS.type == "windows") {
  cl <- parallel::makeCluster(opt$cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  parallel::clusterEvalQ(cl, {
    library(data.table)
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "opt",
      "clean_seq",
      "is_header_line",
      "parse_header",
      "split_labeled_line",
      "read_interleaved_phylip",
      "write_fasta",
      "convert_one",
      "process_one_file"
    ),
    envir = environment()
  )
  
  results <- parallel::parLapply(cl, files, process_one_file)
  
} else {
  results <- parallel::mclapply(
    files,
    process_one_file,
    mc.cores = opt$cores
  )
}

# =============================================================================
# Write logs and R objects
# =============================================================================

summary_dt <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)

error_list <- lapply(results, `[[`, "error")
error_list <- error_list[!vapply(error_list, is.null, logical(1))]

if (length(error_list) > 0L) {
  errors_dt <- rbindlist(error_list, fill = TRUE)
} else {
  errors_dt <- data.table(
    og = character(),
    file = character(),
    error = character()
  )
}

fwrite(
  summary_dt,
  file.path(opt$out_dir, "convert_summary.tsv"),
  sep = "\t"
)

saveRDS(
  summary_dt,
  file.path(opt$out_dir, "convert_summary.rds")
)

if (nrow(errors_dt) > 0L) {
  fwrite(
    errors_dt,
    file.path(opt$out_dir, "convert_errors.tsv"),
    sep = "\t"
  )
  
  saveRDS(
    errors_dt,
    file.path(opt$out_dir, "convert_errors.rds")
  )
}

run_info <- list(
  step = "step0_convert_phylip_to_fasta",
  repo_root = repo_root,
  in_dir = opt$in_dir,
  out_dir = opt$out_dir,
  file_glob = opt$file_glob,
  n_files_requested = length(files),
  cores = opt$cores,
  wrap_width = opt$wrap_width,
  strict_length = opt$strict_length,
  timestamp = Sys.time(),
  session_info = sessionInfo()
)

saveRDS(
  run_info,
  file.path(opt$out_dir, "step0_run_info.rds")
)

# =============================================================================
# Final messages
# =============================================================================

n_ok <- sum(summary_dt$status == "OK", na.rm = TRUE)
n_error <- sum(summary_dt$status == "ERROR", na.rm = TRUE)
n_length_problem <- sum(summary_dt$had_length_problem %in% TRUE, na.rm = TRUE)

message("Done.")
message("Successful conversions: ", n_ok)
message("Errors:                 ", n_error)
message("Length warnings:        ", n_length_problem)
message("Summary:                ", file.path(opt$out_dir, "convert_summary.tsv"))

if (n_error > 0L) {
  warning("Some files failed. Check: ", file.path(opt$out_dir, "convert_errors.tsv"))
}