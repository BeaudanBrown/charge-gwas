library(targets)
library(data.table)
library(dotenv)
library(redcapAPI)

source("R/functions.R")

tar_option_set(
  packages = c("data.table", "parallel", "redcapAPI"),
  format = "rds"
)

load_dot_env()
api_url <- Sys.getenv("REDCAP_API_URL")
redcap <- unlockREDCap(c(bach = "bach"), keyring = "bach", url = api_url)

list(
  tar_target(
    config,
    list(
      imp_folder = Sys.getenv("IMPUTATION_FOLDER"),
      cache_path = Sys.getenv("CACHE_PATH"),
      output_dir = Sys.getenv("OUTPUT_DIR"),
      meta_dir = file.path(Sys.getenv("MERGED_FOLDER"), "meta"),
      bed_dir = file.path(Sys.getenv("MERGED_FOLDER"), "bed"),
      final_dir = file.path(Sys.getenv("MERGED_FOLDER"), "all_chrs"),
      patterns = list(
        model1 = "_model1_results\\.P1\\.assoc\\.linear$",
        model2 = "_model2_results\\.P1\\.assoc\\.linear$",
        freq = "_allele_freq\\.frq$",
        hwe = "_hwe\\.hwe$",
        callrate = "_callrate\\.lmiss$",
        info = "\\.info$"
      ),
      column_order = c(
        "SNPID",
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
        "Oevar_imp"
      )
    )
  ),

  #-----------------------------------------------------------------------------
  # Bash Stage Scripts
  #-----------------------------------------------------------------------------
  tar_target(
    stage1_script,
    "bash/stage1_normalize_merge.sh",
    format = "file"
  ),

  tar_target(
    stage2_script,
    "bash/stage2_prepare_covariates.sh",
    format = "file"
  ),

  tar_target(
    stage3_script,
    "bash/stage3_merge_qc.sh",
    format = "file"
  ),

  tar_target(
    stage4_script,
    "bash/stage4_gwas.sh",
    format = "file"
  ),

  #-----------------------------------------------------------------------------
  # Run Bash Stages
  #-----------------------------------------------------------------------------
  tar_target(
    stage1_outputs,
    {
      result <- system2(
        "bash",
        args = c(stage1_script),
        stdout = TRUE,
        stderr = TRUE
      )
      if (!is.null(attr(result, "status"))) {
        stop("Stage 1 failed with exit code: ", attr(result, "status"))
      }
      list.files(
        config$bed_dir,
        pattern = "chr.*\\.(bed|bim|fam)$",
        full.names = TRUE
      )
    },
    format = "file"
  ),

  #-----------------------------------------------------------------------------
  # REDCap Data Fetch
  #-----------------------------------------------------------------------------
  tar_target(
    redcap_outputs,
    {
      stage1_outputs
      fetch_redcap_data(
        redcap = redcap,
        fam_file = file.path(config$bed_dir, "chr1.fam"),
        output_dir = config$meta_dir
      )
    },
    format = "file"
  ),

  tar_target(
    stage2_outputs,
    {
      stage1_outputs
      redcap_outputs
      result <- system2(
        "bash",
        args = c(stage2_script),
        stdout = TRUE,
        stderr = TRUE
      )
      c(
        file.path(config$meta_dir, "covars_with_batch.txt"),
        file.path(config$meta_dir, "covars2_with_batch.txt"),
        file.path(config$meta_dir, "pheno.txt")
      )
    },
    format = "file"
  ),

  tar_target(
    stage3_outputs,
    {
      stage2_outputs
      result <- system2(
        "bash",
        args = c(stage3_script),
        stdout = TRUE,
        stderr = TRUE
      )
      c(
        file.path(config$final_dir, "all_chrs_clean.bed"),
        file.path(config$final_dir, "all_chrs_clean.bim"),
        file.path(config$final_dir, "all_chrs_clean.fam")
      )
    },
    format = "file"
  ),

  tar_target(
    stage4_outputs,
    {
      stage3_outputs
      result <- system2(
        "bash",
        args = c(stage4_script),
        stdout = TRUE,
        stderr = TRUE
      )
      prefix <- file.path(config$final_dir, "all_chrs")
      c(
        paste0(prefix, "_model1_results.P1.assoc.linear"),
        paste0(prefix, "_model2_results.P1.assoc.linear"),
        paste0(prefix, "_allele_freq.frq"),
        paste0(prefix, "_hwe.hwe"),
        paste0(prefix, "_callrate.lmiss")
      )
    },
    format = "file"
  ),

  #-----------------------------------------------------------------------------
  # Input File Tracking
  #-----------------------------------------------------------------------------
  tar_target(
    gwas_files_model1,
    stage4_outputs[grepl(config$patterns$model1, stage4_outputs)],
    format = "file"
  ),

  tar_target(
    gwas_files_model2,
    stage4_outputs[grepl(config$patterns$model2, stage4_outputs)],
    format = "file"
  ),

  tar_target(
    gwas_files_qc,
    c(
      stage4_outputs[grepl(config$patterns$freq, stage4_outputs)],
      stage4_outputs[grepl(config$patterns$hwe, stage4_outputs)],
      stage4_outputs[grepl(config$patterns$callrate, stage4_outputs)]
    ),
    format = "file"
  ),

  tar_target(
    imputation_files,
    list.files(
      config$imp_folder,
      pattern = config$patterns$info,
      full.names = TRUE
    ),
    format = "file"
  ),

  #-----------------------------------------------------------------------------
  # Read GWAS Results
  #-----------------------------------------------------------------------------
  tar_target(
    model1_results,
    {
      gwas_files_model1
      read_gwas_results(config$bed_dir, config$patterns$model1)
    }
  ),

  tar_target(
    model2_results,
    {
      gwas_files_model2
      read_gwas_results(config$bed_dir, config$patterns$model2)
    }
  ),

  #-----------------------------------------------------------------------------
  # Read and Merge QC Metrics
  #-----------------------------------------------------------------------------
  tar_target(
    qc_metrics,
    {
      gwas_files_qc
      read_qc_metrics(config$bed_dir)
    }
  ),

  tar_target(
    merged_qc,
    merge_qc_metrics(qc_metrics)
  ),

  #-----------------------------------------------------------------------------
  # Load Imputation Info
  #-----------------------------------------------------------------------------
  tar_target(
    imputation_data,
    {
      imputation_files
      pos_info <- model1_results[, .(SNPID, Chr, Position)]
      load_imputation_info(
        config$imp_folder,
        config$cache_path,
        pos_info
      )
    }
  ),

  #-----------------------------------------------------------------------------
  # Create Final Metadata
  #-----------------------------------------------------------------------------
  tar_target(
    final_meta,
    create_final_metadata(imputation_data, merged_qc)
  ),

  #-----------------------------------------------------------------------------
  # Merge Models with Metadata
  #-----------------------------------------------------------------------------
  tar_target(
    model1_merged,
    merge_model_metadata(model1_results, final_meta, config$column_order)
  ),

  tar_target(
    model2_merged,
    merge_model_metadata(model2_results, final_meta, config$column_order)
  ),

  #-----------------------------------------------------------------------------
  # Write Output Files
  #-----------------------------------------------------------------------------
  tar_target(
    model1_output,
    write_gwas_output(
      model1_merged,
      file.path(config$output_dir, "model1.txt")
    ),
    format = "file"
  ),

  tar_target(
    model2_output,
    write_gwas_output(
      model2_merged,
      file.path(config$output_dir, "model2.txt")
    ),
    format = "file"
  )
)
