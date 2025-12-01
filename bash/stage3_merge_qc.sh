#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# STAGE 3: CHROMOSOME MERGING AND QUALITY CONTROL
# ==============================================================================
# This stage:
# 1. Merges all 22 chromosomes into a single PLINK fileset
# 2. Applies quality control filters:
#    - mind 0.05 (individual missingness < 5%)
#    - geno 0.05 (SNP missingness < 5%)
#    - maf 0.01 (minor allele frequency > 1%)
#    - hwe 1e-6 (Hardy-Weinberg equilibrium p > 1e-6)
#    - autosome only
# ==============================================================================

if [ -f .env ]; then
    set -a; source .env; set +a
else
    echo "Error: .env file not found"
    exit 1
fi

: "${MERGED_FOLDER:?}"

mkdir -p ${MERGED_FOLDER}/all_chrs
IN="${MERGED_FOLDER}/bed"
PREFIX="${MERGED_FOLDER}/all_chrs/all_chrs"
MERGE_LIST="${MERGED_FOLDER}/merge_list.txt"

rm -f ${MERGE_LIST}
for i in {2..22}; do
  echo "${IN}/chr${i}.bed ${IN}/chr${i}.bim ${IN}/chr${i}.fam" >> ${MERGE_LIST}
done

# Merge all chromosomes
plink \
  --bfile ${IN}/chr1 \
  --merge-list ${MERGE_LIST} \
  --make-bed \
  --out ${PREFIX}

# Run QC
plink \
  --bfile ${PREFIX} \
  --mind 0.05 \
  --geno 0.05 \
  --maf 0.01 \
  --hwe 1e-6 \
  --autosome \
  --make-bed \
  --out ${PREFIX}_clean
