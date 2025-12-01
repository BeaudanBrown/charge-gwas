#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# STAGE 4: GWAS ANALYSIS (PARALLELIZED)
# ==============================================================================

# Load environment
if [ -f .env ]; then
    set -a; source .env; set +a
else
    echo "Error: .env file not found"
    exit 1
fi

: "${MERGED_FOLDER:?}"

PREFIX="${MERGED_FOLDER}/all_chrs/all_chrs"
META_DIR="${MERGED_FOLDER}/meta"
FINAL_COV1="${META_DIR}/covars_with_batch.txt"
FINAL_COV2="${META_DIR}/covars2_with_batch.txt"
FINAL_PHENO="${META_DIR}/pheno.txt"

plink \
  --bfile ${PREFIX}_clean \
  --covar ${FINAL_COV1} \
  --pheno ${FINAL_PHENO} \
  --linear no-x-sex hide-covar \
  --all-pheno \
  --ci 0.95 \
  --xchr-model 0 \
  --allow-no-sex \
  --out ${PREFIX}_model1_results &
PID1=$!

plink \
  --bfile ${PREFIX}_clean \
  --covar ${FINAL_COV2} \
  --pheno ${FINAL_PHENO} \
  --linear no-x-sex hide-covar \
  --all-pheno \
  --ci 0.95 \
  --xchr-model 0 \
  --allow-no-sex \
  --out ${PREFIX}_model2_results &
PID2=$!

plink \
  --bfile ${PREFIX}_clean \
  --freq \
  --out ${PREFIX}_allele_freq &
PID3=$!

plink \
  --bfile ${PREFIX}_clean \
  --hardy \
  --out ${PREFIX}_hwe &
PID4=$!

plink \
  --bfile ${PREFIX}_clean \
  --missing \
  --out ${PREFIX}_callrate &
PID5=$!

wait $PID1 || true
wait $PID2 || true
wait $PID3 || true
wait $PID4 || true
wait $PID5 || true
