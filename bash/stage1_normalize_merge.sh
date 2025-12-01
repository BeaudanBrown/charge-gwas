#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# STAGE 1: VCF NORMALIZATION AND BATCH MERGING
# ==============================================================================
# This stage:
# 1. Filters VCFs by imputation quality (R2 >= 0.5)
# 2. Normalizes variants using bcftools +fixref and norm
# 3. Merges two batches (B1 and B2)
# 4. Converts to PLINK binary format
# Processes chromosomes 1-22 in parallel using GNU parallel
# ==============================================================================

# Load environment
if [ -f .env ]; then
    set -a; source .env; set +a
else
    echo "Error: .env file not found"
    exit 1
fi

: "${REF:?}" "${B1_FOLDER:?}" "${B2_FOLDER:?}" "${MERGED_FOLDER:?}"
: "${B1_PREFIX:=ILGSA24-20123}"
: "${B2_PREFIX:=ILGSA24-00261}"

mkdir -p ${MERGED_FOLDER}/{b1_norm,b2_norm,merged,bed}

normalise_and_merge() {
  c=$1
  echo "Normalising chr ${c}..."

  in1=${B1_FOLDER}/${B1_PREFIX}_chr${c}_minimac4.dose.vcf.gz
  in2=${B2_FOLDER}/${B2_PREFIX}_chr${c}_minimac4.dose.vcf.gz

  out1=${MERGED_FOLDER}/b1_norm/chr${c}.norm.vcf.gz
  out2=${MERGED_FOLDER}/b2_norm/chr${c}.norm.vcf.gz
  mrg=${MERGED_FOLDER}/merged/chr${c}.vcf.gz

  bcftools view -i 'R2>=0.5' -Oz -o - ${in1} \
  | bcftools +fixref - -- -f ${REF} -m top \
  | bcftools norm -f ${REF} -m -snps -Oz -o ${out1}
  bcftools index -t ${out1}

  bcftools view -i 'R2>=0.5' -Oz -o - ${in2} \
  | bcftools +fixref - -- -f ${REF} -m top \
  | bcftools norm -f ${REF} -m -snps -Oz -o ${out2}
  bcftools index -t ${out2}

  # Merge
  bcftools merge -0 -m none -Oz -o ${mrg} ${out1} ${out2}
  bcftools index -t ${mrg}

  # Convert to bed
  plink \
    --vcf ${mrg} \
    --keep-allele-order \
    --set-missing-var-ids '@:#:$1:$2' \
    --double-id \
    --make-bed \
    --out ${MERGED_FOLDER}/bed/chr${c} \
    --memory 5000
}

export -f normalise_and_merge
export REF B1_FOLDER B2_FOLDER MERGED_FOLDER B1_PREFIX B2_PREFIX

parallel --ungroup normalise_and_merge ::: {1..22}
