#!/usr/bin/env bash
set -euo pipefail

# adjust these if needed
WORKDIR=bed
COV1=covars.txt
COV2=covars2.txt
PHENO=pheno.txt

process_chromosome() {
  chr=$1
  echo "---- Processing chr${chr} ----"
  PREFIX=${WORKDIR}/chr${chr}

  # 1) sample‐missing filter
  plink \
    --bfile ${PREFIX} \
    --mind 0.05 \
    --make-bed \
    --out ${PREFIX}_qc1

  # 2) marker‐missingness filter
  plink \
    --bfile ${PREFIX}_qc1 \
    --geno 0.05 \
    --make-bed \
    --out ${PREFIX}_qc2

  # 3) MAF filter
  plink \
    --bfile ${PREFIX}_qc2 \
    --maf 0.01 \
    --make-bed \
    --out ${PREFIX}_qc3

  # 4) HWE filter
  plink \
    --bfile ${PREFIX}_qc3 \
    --hwe 1e-6 \
    --make-bed \
    --out ${PREFIX}_qc4

  # 5) keep autosomes only
  plink \
    --bfile ${PREFIX}_qc4 \
    --autosome \
    --make-bed \
    --out ${PREFIX}_qc5

  # ===== final analyses on qc5 =====
  plink \
    --bfile ${PREFIX}_qc5 \
    --covar ${COV1} \
    --pheno ${PHENO} \
    --linear no-x-sex hide-covar \
    --ci 0.95 \
    --xchr-model 0 \
    --allow-no-sex \
    --out ${PREFIX}_model1_results

  plink \
    --bfile ${PREFIX}_qc5 \
    --covar ${COV2} \
    --pheno ${PHENO} \
    --linear no-x-sex hide-covar \
    --ci 0.95 \
    --xchr-model 0 \
    --allow-no-sex \
    --out ${PREFIX}_model2_results

  plink \
    --bfile ${PREFIX}_qc5 \
    --freq \
    --out ${PREFIX}_allele_freq

  plink \
    --bfile ${PREFIX}_qc5 \
    --hardy \
    --out ${PREFIX}_hwe

  plink \
    --bfile ${PREFIX}_qc5 \
    --missing \
    --out ${PREFIX}_callrate

  # ===== clean up all QC intermediates =====
  rm -f ${PREFIX}_qc1.* ${PREFIX}_qc2.* ${PREFIX}_qc3.* ${PREFIX}_qc4.* ${PREFIX}_qc5.*
}

# Export the function and variables so parallel can see them
export -f process_chromosome
export WORKDIR COV1 COV2 PHENO

# Run in parallel
parallel process_chromosome ::: {1..22}
