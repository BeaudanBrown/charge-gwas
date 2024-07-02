library(data.table)
library(dotenv)

dotenv::load_dot_env()
imp_folder <- Sys.getenv("IMPUTATION_FOLDER")
gwas_folder <- Sys.getenv("GWAS_FOLDER")

model1_results <-
  fread(file.path(gwas_folder, "model1_results.assoc.linear"), header = TRUE)
model1_results <-
  model1_results[, .(SNPID = SNP,
                     n_total = NMISS,
                     Chr = CHR,
                     Position = BP,
                     Beta = BETA,
                     SE = SE,
                     Pval = P)]

model2_results <-
  fread(file.path(gwas_folder, "model2_results.assoc.linear"),
        header = TRUE)
model2_results <-
  model2_results[, .(SNPID = SNP,
                     n_total = NMISS,
                     Chr = CHR,
                     Position = BP,
                     Beta = BETA,
                     SE = SE,
                     Pval = P)]

# library(qqman) # nolint: commented_code_linter.
# manhattan(model1_results) # nolint: commented_code_linter.
# manhattan(model2_results) # nolint: commented_code_linter.

freq_data <- fread(file.path(gwas_folder, "allele_freq.frq"), header = TRUE)
hwe_data <- fread(file.path(gwas_folder, "hwe.hwe"), header = TRUE)
callrate_data <- fread(file.path(gwas_folder, "callrate.lmiss"), header = TRUE)

freq_info <- freq_data[, .(SNPID = SNP,
                           AF_coded_all = MAF)]

hwe_info <- hwe_data[, .(SNPID = SNP,
                         HWE_pval = P)]

callrate_info <- callrate_data[, .(SNPID = SNP,
                                   Callrate = 1 - F_MISS)]

# Merge all plink data
merged_data <- merge(hwe_info, freq_info, by = "SNPID")
merged_data <- merge(merged_data, callrate_info, by = "SNPID")

imputation_path <- file.path(gwas_folder, "imputation_data.rds")
if (file.exists(imputation_path)) {
  imputation_data <- readRDS(imputation_path)
} else {
  imputation_data <- data.table()
  # Parse all chromosome info files to get imputation statistics
  imputation_info_files <- list.files(
    path = imp_folder,
    pattern = "\\.info$",
    full.names = TRUE
  )
  for (file in imputation_info_files) {
    imputation_info <- fread(file)
    # Extract the CHR and POS columns
    split_snp <- strsplit(as.character(imputation_info$SNP), ":")
    imputation_info$Chr <- sapply(split_snp, `[`, 1)
    imputation_info$Position <- sapply(split_snp, `[`, 2)

    imputation_info <- imputation_info[, .(
      Chr = as.integer(Chr),
      Position = as.integer(Position),
      Imputed = ifelse(Genotyped == "Imputed", 1, 0),
      Used_for_imp = ifelse(Genotyped == "Genotyped", 1, 0),
      Coded_all = `ALT(1)`,
      Noncoded_all = `REF(0)`
    )]

    # Filter to only include SNPs which passed QC and include SNPID
    pos_info <- model1_results[, .(SNPID = SNPID,
                                   Chr = Chr,
                                   Position = Position)]
    filtered_imputation_info <-
      imputation_info[pos_info, on = .(Chr, Position), nomatch = 0]

    imputation_data <-
      rbind(imputation_data, filtered_imputation_info, fill = TRUE)
  }
  # Save RDS to prevent recalculation
  saveRDS(imputation_data, file.path("./tmp/imputation_data.rds"))
}

# Remove SNPs that don't have an rsID
imputation_data <- imputation_data[, .(SNPID,
                                       Imputed,
                                       Used_for_imp,
                                       Coded_all,
                                       Noncoded_all)]

merged_data <- merge(imputation_data, merged_data, by = c("SNPID"))
merged_data[, Strand_genome := "+"]
merged_data[, Oevar_imp := NA]
merged_data[Imputed == 1, HWE_pval := NA]

# Compine with gwas results
model1_merged <- merge(model1_results, merged_data, by = c("SNPID"))
model2_merged <- merge(model2_results, merged_data, by = c("SNPID"))

# Only have SNPID for rsIDs
model1_merged[!grepl("^rs", SNPID), SNPID := ""]
model2_merged[!grepl("^rs", SNPID), SNPID := ""]

column_order <- c("SNPID",
                  "Chr",
                  "Position",
                  "Coded_all",
                  "Noncoded_all",
                  "Strand_genome",
                  "Beta",
                  "SE",
                  "Pval",
                  "AF_coded_all",
                  "HWE_pval",
                  "Callrate",
                  "n_total",
                  "Imputed",
                  "Used_for_imp",
                  "Oevar_imp")

# Reorder the columns
setcolorder(model1_merged, column_order)
setorder(model1_merged, Chr, Position)

setcolorder(model2_merged, column_order)
setorder(model2_merged, Chr, Position)

fwrite(model1_merged,
       file = file.path(gwas_folder, "model1.txt"),
       sep = "\t",
       na = "NA",
       quote = FALSE)
fwrite(model2_merged,
       file = file.path(gwas_folder, "model2.txt"),
       sep = "\t",
       na = "NA",
       quote = FALSE)
