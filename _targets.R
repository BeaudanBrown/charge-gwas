#===============================================================================
# Targets Pipeline for GWAS Results Processing
#===============================================================================
#
# This pipeline processes GWAS results from PLINK, merges them with QC metrics
# and imputation info, and produces final formatted output tables.
#
# Usage:
#   targets::tar_make()                      # Run main pipeline
#   targets::tar_make(names = "manhattan_model1")  # Include manhattan plot
#   targets::tar_visnetwork()                # Visualize pipeline
#   targets::tar_read(model1_merged)         # Read results
#
#===============================================================================

library(targets)
library(data.table)
library(dotenv)

# Load helper functions
source("R/functions.R")

# Set target options
tar_option_set(
  packages = c("data.table", "parallel"),
  format = "rds"
)

# Load configuration from .env
load_dot_env()

#===============================================================================
# Pipeline Definition
#===============================================================================

list(
  #-----------------------------------------------------------------------------
  # Configuration
  #-----------------------------------------------------------------------------
  tar_target(
    config,
    list(
      imp_folder = Sys.getenv("IMPUTATION_FOLDER"),
      gwas_folder = Sys.getenv("GWAS_FOLDER"),
      cache_path = "./tmp/imputation_data.rds",
      output_dir = "./results",
      cores = 11,
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
  # Input File Tracking (automatically detects file changes)
  #-----------------------------------------------------------------------------
  tar_target(
    gwas_files_model1,
    list.files(
      config$gwas_folder,
      pattern = config$patterns$model1,
      full.names = TRUE
    ),
    format = "file"
  ),

  tar_target(
    gwas_files_model2,
    list.files(
      config$gwas_folder,
      pattern = config$patterns$model2,
      full.names = TRUE
    ),
    format = "file"
  ),

  tar_target(
    gwas_files_qc,
    c(
      list.files(
        config$gwas_folder,
        pattern = config$patterns$freq,
        full.names = TRUE
      ),
      list.files(
        config$gwas_folder,
        pattern = config$patterns$hwe,
        full.names = TRUE
      ),
      list.files(
        config$gwas_folder,
        pattern = config$patterns$callrate,
        full.names = TRUE
      )
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
  # Read GWAS Results (can run in parallel)
  #-----------------------------------------------------------------------------
  tar_target(
    model1_results,
    {
      gwas_files_model1 # Dependency to track file changes
      read_gwas_results(config$gwas_folder, config$patterns$model1)
    }
  ),

  tar_target(
    model2_results,
    {
      gwas_files_model2 # Dependency to track file changes
      read_gwas_results(config$gwas_folder, config$patterns$model2)
    }
  ),

  #-----------------------------------------------------------------------------
  # Read and Merge QC Metrics
  #-----------------------------------------------------------------------------
  tar_target(
    qc_metrics,
    {
      gwas_files_qc # Dependency to track file changes
      read_qc_metrics(config$gwas_folder)
    }
  ),

  tar_target(
    merged_qc,
    merge_qc_metrics(qc_metrics)
  ),

  #-----------------------------------------------------------------------------
  # Load Imputation Info (with caching)
  #-----------------------------------------------------------------------------
  tar_target(
    imputation_data,
    {
      imputation_files # Dependency to track file changes
      pos_info <- model1_results[, .(SNPID, Chr, Position)]
      load_imputation_info(
        config$imp_folder,
        config$cache_path,
        pos_info,
        config$cores
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
  ),

  #-----------------------------------------------------------------------------
  # Manhattan Plots (Optional - only run if explicitly requested)
  #-----------------------------------------------------------------------------
  tar_target(
    manhattan_model1,
    create_manhattan_plot(
      model1_merged,
      file.path(config$output_dir, "manhattan_model1.png")
    ),
    format = "file",
    cue = tar_cue("never") # Never runs automatically
  ),

  tar_target(
    manhattan_model2,
    create_manhattan_plot(
      model2_merged,
      file.path(config$output_dir, "manhattan_model2.png")
    ),
    format = "file",
    cue = tar_cue("never") # Never runs automatically
  )
)
