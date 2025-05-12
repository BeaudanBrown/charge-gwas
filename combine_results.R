# load libraries
library(data.table)
library(dotenv)

# load env vars
dotenv::load_dot_env()
imp_folder  <- Sys.getenv("IMPUTATION_FOLDER")
gwas_folder <- Sys.getenv("GWAS_FOLDER")

#-------------------------------------------------------------------------------
# 1) Helper to read & bind a set of files matching a pattern
#-------------------------------------------------------------------------------
read_and_bind <- function(path, pattern, transform_fun = identity) {
  files <- list.files(path = path,
                      pattern = pattern,
                      full.names = TRUE)
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
# 2) Read & reshape model1 results across all chrs
#-------------------------------------------------------------------------------
model1_results <- read_and_bind(

  path    = gwas_folder,
  pattern = "_model1_results\\.assoc\\.linear$",
  transform_fun = function(dt) {
    dt[, .(SNPID    = SNP,
           n_total  = NMISS,
           Chr      = CHR,
           Position = BP,
           Beta     = BETA,
           SE       = SE,
           Pval     = P)]
  }
)
# library(qqman) # nolint: commented_code_linter.
# manhattan(
#   model1_results,
#   chr = "Chr",
#   bp = "Position",
#   p = "Pval",
#   snp = "SNPID",
# )

#-------------------------------------------------------------------------------
# 3) Read & reshape model2 results across all chrs
#-------------------------------------------------------------------------------
model2_results <- read_and_bind(
  path    = gwas_folder,
  pattern = "_model2_results\\.assoc\\.linear$",
  transform_fun = function(dt) {
    dt[, .( SNPID    = SNP,
            n_total  = NMISS,
            Chr      = CHR,
            Position = BP,
            Beta     = BETA,
            SE       = SE,
            Pval     = P )]
  }
)

#-------------------------------------------------------------------------------
# 4) Read & reshape allele frequency, HWE and call‐rate across all chrs
#-------------------------------------------------------------------------------
freq_info <- read_and_bind(
  path    = gwas_folder,
  pattern = "_allele_freq\\.frq$",
  transform_fun = function(dt) {
    dt[, .( SNPID         = SNP,
            AF_coded_all  = MAF )]
  }
)

hwe_info <- read_and_bind(
  path    = gwas_folder,
  pattern = "_hwe\\.hwe$",
  transform_fun = function(dt) {
    dt[, .( SNPID    = SNP,
            HWE_pval = P )]
  }
)

callrate_info <- read_and_bind(
  path    = gwas_folder,
  pattern = "_callrate\\.lmiss$",
  transform_fun = function(dt) {
    dt[, .( SNPID   = SNP,
            Callrate = 1 - F_MISS )]
  }
)

# merge frequency/HWE/callrate
merged_data <- merge(hwe_info,    freq_info,    by = "SNPID")
merged_data <- merge(merged_data, callrate_info, by = "SNPID")

#-------------------------------------------------------------------------------
# 5) Load or build imputation info
#-------------------------------------------------------------------------------
imputation_rds <- file.path("./tmp/imputation_data.rds")
if (file.exists(imputation_rds)) {
  imputation_data <- readRDS(imputation_rds)
} else {
  # read all .info files in imp_folder
  info_files <- list.files(path       = imp_folder,
                           pattern    = "\\.info$",
                           full.names = TRUE)
  imputation_data <- rbindlist(lapply(info_files, function(f) {
    dt <- fread(f)
    # split SNP column into Chr and Position
    tmp <- tstrsplit(dt$SNP, ":", keep = 1:2)
    dt[, Chr      := as.integer(tmp[[1]])]
    dt[, Position := as.integer(tmp[[2]])]
    # recode
    dt[, Imputed      := fifelse(Genotyped == "Imputed", 1L, 0L)]
    dt[, Used_for_imp := fifelse(Genotyped == "Genotyped", 1L, 0L)]
    setnames(dt,
             old = c("ALT(1)", "REF(0)"),
             new = c("Coded_all","Noncoded_all"))
    dt[, .(Chr,
           Position,
           Imputed,
           Used_for_imp,
           Coded_all,
           Noncoded_all)]
  }), use.names = TRUE, fill = TRUE)

  pos_info <- model1_results[, .(SNPID = SNPID,
                                 Chr = Chr,
                                 Position = Position)]
  imputation_data <-
    imputation_data[pos_info, on = .(Chr, Position), nomatch = 0]

  dir.create(dirname(imputation_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(imputation_data, imputation_rds)
}

#-------------------------------------------------------------------------------
# 6) Final merges & column ordering
#-------------------------------------------------------------------------------

# combine imputation + QC metrics
setkey(imputation_data,   Chr, Position)
setkey(merged_data,       SNPID)
setkey(model1_results,    SNPID)
setkey(model2_results,    SNPID)

# join by SNPID
final_meta <- merge(imputation_data,
                    merged_data,
                    by = "SNPID",
                    all = FALSE)[, -c("Chr", "Position")]
final_meta[, Strand_genome := "+"]
final_meta[, Oevar_imp := NA]
final_meta[Imputed == 1, HWE_pval := NA]

# rows for each model
model1_merged <- merge(model1_results, final_meta, by = "SNPID")
model2_merged <- merge(model2_results, final_meta, by = "SNPID")

# Only have SNPID for rsIDs
model1_merged[!grepl("^rs", SNPID), SNPID := ""]
model2_merged[!grepl("^rs", SNPID), SNPID := ""]

# column order
column_order <- c(
  "SNPID", "Chr", "Position",
  "Coded_all","Noncoded_all","Strand_genome",
  "Beta","SE","Pval","AF_coded_all","HWE_pval","Callrate","n_total",
  "Imputed","Used_for_imp","Oevar_imp"
)

# Reorder the columns
setcolorder(model1_merged, column_order)
setorder(model1_merged, Chr, Position)

setcolorder(model2_merged, column_order)
setorder(model2_merged, Chr, Position)
test <- model1_merged[model1_merged$Pval < 0.000005, ]

#-------------------------------------------------------------------------------
# 7) Write out final tables
#-------------------------------------------------------------------------------
fwrite(model1_merged,
       file   = file.path(gwas_folder, "model1.txt"),
       sep    = "\t", na = "NA", quote = FALSE)

fwrite(model2_merged,
       file   = file.path(gwas_folder, "model2.txt"),
       sep    = "\t", na = "NA", quote = FALSE)
