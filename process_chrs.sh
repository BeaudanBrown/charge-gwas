#!/usr/bin/env bash
set -euo pipefail

# adjust these if needed
WORKDIR=qc
COV1=covars.txt
COV2=covars2.txt
PHENO=pheno.txt

for chr in {1..22}; do
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

  # 5) keep autosomes only (optional if your files are already per‐chr)
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

done

echo "Done.  For each chromosome 1–22 you now have exactly these five files:"
echo "  chrN_model1_results.*"
echo "  chrN_model2_results.*"
echo "  chrN_allele_freq.*"
echo "  chrN_hwe.*"
echo "  chrN_callrate.*"
