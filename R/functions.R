#===============================================================================
# Helper Functions for GWAS Results Processing Pipeline
#===============================================================================

library(data.table)

#-------------------------------------------------------------------------------
# Read and bind files matching a pattern
#-------------------------------------------------------------------------------
read_and_bind <- function(path, pattern, transform_fun = identity) {
  files <- list.files(path = path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop("No files found matching ", pattern, " in ", path)
  }
  dt_list <- lapply(files, function(f) {
    dt <- fread(f)
    transform_fun(dt)
  })
  rbindlist(dt_list, use.names = TRUE)
}

#-------------------------------------------------------------------------------
# Read GWAS results for a specific model
#-------------------------------------------------------------------------------
read_gwas_results <- function(gwas_folder, pattern) {
  read_and_bind(
    path = gwas_folder,
    pattern = pattern,
    transform_fun = function(dt) {
      dt[, .(
        SNPID = SNP,
        n_total = NMISS,
        Chr = CHR,
        Position = BP,
        Beta = BETA,
        SE = SE,
        Pval = P
      )]
    }
  )
}

#-------------------------------------------------------------------------------
# Read all QC metrics (frequency, HWE, callrate)
#-------------------------------------------------------------------------------
read_qc_metrics <- function(gwas_folder) {
  freq_info <- read_and_bind(
    path = gwas_folder,
    pattern = "_allele_freq\\.frq$",
    transform_fun = function(dt) {
      dt[, .(SNPID = SNP, AF_coded_all = MAF)]
    }
  )

  hwe_info <- read_and_bind(
    path = gwas_folder,
    pattern = "_hwe\\.hwe$",
    transform_fun = function(dt) {
      dt[, .(SNPID = SNP, Coded_all = A1, Noncoded_all = A2, HWE_pval = P)]
    }
  )

  callrate_info <- read_and_bind(
    path = gwas_folder,
    pattern = "_callrate\\.lmiss$",
    transform_fun = function(dt) {
      dt[, .(SNPID = SNP, Callrate = 1 - F_MISS)]
    }
  )

  list(
    freq = freq_info,
    hwe = hwe_info,
    callrate = callrate_info
  )
}

#-------------------------------------------------------------------------------
# Merge QC metrics into single data.table
#-------------------------------------------------------------------------------
merge_qc_metrics <- function(qc_list) {
  merged_data <- merge(qc_list$hwe, qc_list$freq, by = "SNPID")
  merged_data <- merge(merged_data, qc_list$callrate, by = "SNPID")
  merged_data
}

#-------------------------------------------------------------------------------
# Load or build imputation info (with caching)
#-------------------------------------------------------------------------------
load_imputation_info <- function(imp_folder, cache_path, pos_info, cores = 11) {
  if (file.exists(cache_path)) {
    message("Loading cached imputation data from: ", cache_path)
    imputation_data <- readRDS(cache_path)
  } else {
    message("Building imputation data from .info files...")

    # Read all .info files in imp_folder
    info_files <- list.files(
      path = imp_folder,
      pattern = "\\.info$",
      full.names = TRUE
    )

    if (length(info_files) == 0L) {
      stop("No .info files found in ", imp_folder)
    }

    library(parallel)
    imputation_data <- rbindlist(
      mclapply(
        info_files,
        function(f) {
          dt <- fread(f)
          # Split SNP column into Chr and Position
          tmp <- tstrsplit(dt$SNP, ":", keep = 1:2)
          dt[, Chr := as.integer(tmp[[1]])]
          dt[, Position := as.integer(tmp[[2]])]
          # Recode imputation status
          dt[, Imputed := fifelse(Genotyped == "Imputed", 1L, 0L)]
          dt[, Used_for_imp := fifelse(Genotyped == "Genotyped", 1L, 0L)]
          # Rename allele columns
          setnames(
            dt,
            old = c("ALT(1)", "REF(0)"),
            new = c("Coded_all", "Noncoded_all")
          )
          dt[, .(Chr, Position, Imputed, Used_for_imp, Coded_all, Noncoded_all)]
        },
        mc.cores = cores
      ),
      use.names = TRUE,
      fill = TRUE
    )

    # Join with position info to get SNPID
    imputation_data <- imputation_data[
      pos_info,
      on = .(Chr, Position),
      nomatch = 0
    ]

    # Save cache
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(imputation_data, cache_path)
    message("Cached imputation data to: ", cache_path)
  }

  imputation_data
}

#-------------------------------------------------------------------------------
# Create final metadata by merging imputation and QC data
#-------------------------------------------------------------------------------
create_final_metadata <- function(imputation_data, merged_qc) {
  # Filter to biallelic SNPs only (single nucleotide)
  imputation_data <- imputation_data[
    nchar(Coded_all) == 1 & nchar(Noncoded_all) == 1
  ]

  # Set keys for efficient merging
  setkey(imputation_data, Chr, Position)
  setkey(merged_qc, SNPID)

  # Join by SNPID
  final_meta <- merge(imputation_data, merged_qc, by = "SNPID", all = FALSE)[,
    -c("Chr", "Position")
  ]

  # Add constant fields
  final_meta[, Strand_genome := "+"]
  final_meta[, Oevar_imp := NA]

  # Set HWE p-value to NA for imputed SNPs (HWE doesn't apply to imputed variants)
  final_meta[Imputed == 1, HWE_pval := NA]

  final_meta
}

#-------------------------------------------------------------------------------
# Merge model results with metadata and format output
#-------------------------------------------------------------------------------
merge_model_metadata <- function(model_results, final_meta, column_order) {
  # Set key for efficient merging
  setkey(model_results, SNPID)

  # Merge model results with metadata
  merged <- merge(model_results, final_meta, by = "SNPID")

  # Blank out non-rsID SNPIDs (we don't have rsIDs in this dataset)
  merged[!grepl("^rs", SNPID), SNPID := ""]

  # Reorder columns
  setcolorder(merged, column_order)

  # Sort by chromosome and position
  setorder(merged, Chr, Position)

  merged
}

#-------------------------------------------------------------------------------
# Write GWAS output file
#-------------------------------------------------------------------------------
write_gwas_output <- function(data, output_path) {
  # Ensure directory exists
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  fwrite(
    data,
    file = output_path,
    sep = "\t",
    na = "NA",
    quote = FALSE
  )

  message("Wrote output to: ", output_path)
  output_path
}

#-------------------------------------------------------------------------------
# Create Manhattan plot (optional)
#-------------------------------------------------------------------------------
create_manhattan_plot <- function(gwas_data, output_path) {
  if (!requireNamespace("qqman", quietly = TRUE)) {
    stop(
      "Package 'qqman' is required for manhattan plots. Install with: install.packages('qqman')"
    )
  }

  # Ensure directory exists
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

  png(output_path, width = 1200, height = 600, res = 100)
  qqman::manhattan(
    gwas_data,
    chr = "Chr",
    bp = "Position",
    p = "Pval",
    snp = "SNPID"
  )
  dev.off()

  message("Created Manhattan plot: ", output_path)
  output_path
}
