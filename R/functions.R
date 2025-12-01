#-------------------------------------------------------------------------------
# Fetch REDCap data and create covariate/phenotype files
#-------------------------------------------------------------------------------
fetch_redcap_data <- function(redcap, fam_file, output_dir = ".") {
  result <- exportRecordsTyped(
    rcon = redcap$bach,
    fields = c("idno", "sex", "age", "education", "vrii_total_raw"),
    events = c("baseline_arm_1"),
  )

  records <- as.data.table(result)

  # Remove invalid records
  records <- records[grepl("--1$", idno)]
  records[, idno := sub("--1$", "", idno)]
  records <- records[grepl("^0.{3}$", idno)]

  records[, sex := as.numeric(sex)]

  fam_data <- fread(
    fam_file,
    header = FALSE,
    col.names = c("FID", "IID", "V3", "V4", "V5", "V6")
  )
  fam_data[, idno := tstrsplit(IID, "BACH", fixed = TRUE, keep = 2)]
  matched_records <- records[fam_data, on = .(idno), nomatch = 0]
  setorder(matched_records, FID)

  pheno_file <- file.path(output_dir, "pheno.txt")
  fwrite(
    matched_records[, .(FID, IID, vrii_total_raw)],
    file = pheno_file,
    sep = " ",
    col.names = FALSE,
    quote = FALSE
  )

  # model 1: sex, age
  covars_file <- file.path(output_dir, "covars.txt")
  fwrite(
    matched_records[, .(FID, IID, sex, age)],
    file = covars_file,
    sep = " ",
    col.names = FALSE,
    quote = FALSE
  )

  # model 2: sex, age, education
  covars2_file <- file.path(output_dir, "covars2.txt")
  fwrite(
    matched_records[, .(FID, IID, sex, age, education)],
    file = covars2_file,
    sep = " ",
    col.names = FALSE,
    quote = FALSE
  )

  c(pheno_file, covars_file, covars2_file)
}

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
# Read all QC metrics
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
# Merge QC metrics
#-------------------------------------------------------------------------------
merge_qc_metrics <- function(qc_list) {
  merged_data <- merge(qc_list$hwe, qc_list$freq, by = "SNPID")
  merged_data <- merge(merged_data, qc_list$callrate, by = "SNPID")
  merged_data
}

#-------------------------------------------------------------------------------
# Load or build imputation info
#-------------------------------------------------------------------------------
load_imputation_info <- function(imp_folder, cache_path, pos_info, cores = 11) {
  if (file.exists(cache_path)) {
    imputation_data <- readRDS(cache_path)
  } else {
    # Read all .info files in imp_folder
    info_files <- list.files(
      path = imp_folder,
      pattern = "\\.info$",
      full.names = TRUE
    )

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
          dt[, .(Chr, Position, Imputed, Used_for_imp)]
        },
        mc.cores = cores
      ),
      use.names = TRUE,
      fill = TRUE
    )

    imputation_data <- imputation_data[
      pos_info,
      on = .(Chr, Position),
      nomatch = 0
    ]

    # Save cache
    dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
    saveRDS(imputation_data, cache_path)
  }

  imputation_data[, .(Chr, Position, Imputed, Used_for_imp, SNPID)]
}

#-------------------------------------------------------------------------------
# Create final metadata
#-------------------------------------------------------------------------------
create_final_metadata <- function(imputation_data, merged_qc) {
  # Filter to biallelic SNPs only
  merged_qc <- merged_qc[
    nchar(Coded_all) == 1 & nchar(Noncoded_all) == 1
  ]

  setkey(imputation_data, Chr, Position)
  setkey(merged_qc, SNPID)

  final_meta <- merge(imputation_data, merged_qc, by = "SNPID", all = FALSE)[,
    -c("Chr", "Position")
  ]

  # Add constant fields
  final_meta[, Strand_genome := "+"]
  final_meta[, Oevar_imp := NA]
  # HWE doesn't apply to imputed variants
  final_meta[Imputed == 1, HWE_pval := NA]

  final_meta
}

#-------------------------------------------------------------------------------
# Merge model results with metadata and format output
#-------------------------------------------------------------------------------
merge_model_metadata <- function(model_results, final_meta, column_order) {
  setkey(model_results, SNPID)
  merged <- merge(model_results, final_meta, by = "SNPID")
  merged[!grepl("^rs", SNPID), SNPID := ""]
  setcolorder(merged, column_order)
  setorder(merged, Chr, Position)
  merged
}

#-------------------------------------------------------------------------------
# Write GWAS output file
#-------------------------------------------------------------------------------
write_gwas_output <- function(data, output_path) {
  dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
  fwrite(
    data,
    file = output_path,
    sep = "\t",
    na = "NA",
    quote = FALSE
  )
  output_path
}

#-------------------------------------------------------------------------------
# Create Manhattan plot (optional)
#-------------------------------------------------------------------------------
create_manhattan_plot <- function(gwas_data, output_path) {
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
  output_path
}
