#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# STAGE 2: COVARIATE AND BATCH MAPPING PREPARATION
# ==============================================================================
# This stage:
# 1. Extracts sample lists from VCF files
# 2. Creates batch identifier mapping (1=Batch1, 2=Batch2)
# 3. Merges batch identifiers into covariate files
# ==============================================================================

if [ -f .env ]; then
    set -a; source .env; set +a
else
    echo "Error: .env file not found"
    exit 1
fi

: "${B1_FOLDER:?}" "${B2_FOLDER:?}" "${MERGED_FOLDER:?}"
: "${B1_PREFIX:=ILGSA24-20123}"
: "${B2_PREFIX:=ILGSA24-00261}"
: "${COV1:?}" "${COV2:?}"

META_DIR="${MERGED_FOLDER}/meta"
mkdir -p ${META_DIR}

# Label samples with batch number
bcftools query -l ${B1_FOLDER}/${B1_PREFIX}_chr22_minimac4.dose.vcf.gz > ${META_DIR}/b1.samples
bcftools query -l ${B2_FOLDER}/${B2_PREFIX}_chr22_minimac4.dose.vcf.gz > ${META_DIR}/b2.samples
awk '{print $1,1}' ${META_DIR}/b1.samples > ${META_DIR}/batch_map.txt
awk '{print $1,2}' ${META_DIR}/b2.samples >> ${META_DIR}/batch_map.txt

# Sort for Joining
sort -k2,2 ${COV1} > ${META_DIR}/cov1.tmp
sort -k2,2 ${COV2} > ${META_DIR}/cov2.tmp
sort -k1,1 ${META_DIR}/batch_map.txt > ${META_DIR}/batch.tmp

# Model 1 covariates (age, sex, batch)
awk 'NR==FNR {gsub(/\r$/, "", $2); batch[$1]=$2; next}
     {gsub(/\r$/, ""); print $1, $2, $3, $4, batch[$2]}' \
     ${META_DIR}/batch.tmp ${META_DIR}/cov1.tmp > ${META_DIR}/covars_with_batch.txt

# Model 2 covariates (age, sex, education, batch)
awk 'NR==FNR {gsub(/\r$/, "", $2); batch[$1]=$2; next}
     {gsub(/\r$/, ""); print $1, $2, $3, $4, $5, batch[$2]}' \
     ${META_DIR}/batch.tmp ${META_DIR}/cov2.tmp > ${META_DIR}/covars2_with_batch.txt
