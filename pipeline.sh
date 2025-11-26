#!/usr/bin/env bash
set -euo pipefail

# 1. Safer Environment Loading
if [ -f .env ]; then
    set -a; source .env; set +a
else
    echo "Error: .env file not found"
    exit 1
fi

# Ensure required variables are set
: "${REF:?}" "${B1_FOLDER:?}" "${B2_FOLDER:?}" "${OUT_FOLDER:?}"
: "${COV1:?}" "${COV2:?}" "${PHENO:?}"

# Define Prefixes (Adjust these in .env if they differ)
# Example: B1_PREFIX="ILGSA24-20123"
: "${B1_PREFIX:=ILGSA24-20123}"
: "${B2_PREFIX:=ILGSA24-00261}"

# Setup directories
# samtools faidx ${REF}
mkdir -p ${OUT_FOLDER}/{b1_norm,b2_norm,merged,bed,meta}

# Function to process VCFs
normalise_and_merge() {
  c=$1
  echo "Normalising chr ${c}..."
  
  # Use variables for filenames
  in1=${B1_FOLDER}/${B1_PREFIX}_chr${c}_minimac4.dose.vcf.gz
  in2=${B2_FOLDER}/${B2_PREFIX}_chr${c}_minimac4.dose.vcf.gz
  
  out1=${OUT_FOLDER}/b1_norm/chr${c}.norm.vcf.gz
  out2=${OUT_FOLDER}/b2_norm/chr${c}.norm.vcf.gz
  mrg=${OUT_FOLDER}/merged/chr${c}.vcf.gz

  # Batch 1: filter + normalize
  bcftools view -i 'R2>=0.5' -Oz -o - ${in1} \
  | bcftools +fixref - -- -f ${REF} -m top \
  | bcftools norm -f ${REF} -m -snps -Oz -o ${out1}
  bcftools index -t ${out1}

  # Batch 2: filter + normalize
  bcftools view -i 'R2>=0.5' -Oz -o - ${in2} \
  | bcftools +fixref - -- -f ${REF} -m top \
  | bcftools norm -f ${REF} -m -snps -Oz -o ${out2}
  bcftools index -t ${out2}

  # Merge
  bcftools merge -0 -m none -Oz -o ${mrg} ${out1} ${out2}
  bcftools index -t ${mrg}

  plink \
    --vcf ${mrg} \
    --keep-allele-order \
    --set-missing-var-ids '@:#:$1:$2' \
    --double-id \
    --make-bed \
    --out ${OUT_FOLDER}/bed/chr${c} \
    --memory 5000
}

export -f normalise_and_merge
export REF B1_FOLDER B2_FOLDER OUT_FOLDER B1_PREFIX B2_PREFIX

echo "Normalising and merging batches..."
parallel --ungroup normalise_and_merge ::: {1..22}

# ==============================================================================
# INTEGRATION: COVARIATE & BATCH MAPPING
# ==============================================================================
echo "Processing Covariates and Batch Maps..."

META_DIR="${OUT_FOLDER}/meta"

# 1. Extract Samples from a representative chromosome (e.g., chr22) to determine batch
bcftools query -l ${B1_FOLDER}/${B1_PREFIX}_chr22_minimac4.dose.vcf.gz > ${META_DIR}/b1.samples
bcftools query -l ${B2_FOLDER}/${B2_PREFIX}_chr22_minimac4.dose.vcf.gz > ${META_DIR}/b2.samples

# 2. Create Batch Map (1=Batch1, 2=Batch2)
awk '{print $1,1}' ${META_DIR}/b1.samples > ${META_DIR}/batch_map.txt
awk '{print $1,2}' ${META_DIR}/b2.samples >> ${META_DIR}/batch_map.txt

# 3. Format Covariates and Phenotypes
awk '{print $1"_"$2, $1"_"$2, $3, $4}' ${COV1} > ${META_DIR}/covars_fid_iid.txt
awk '{print $1"_"$2, $1"_"$2, $3, $4, $5}' ${COV2} > ${META_DIR}/covars2_fid_iid.txt
awk '{print $1"_"$2, $1"_"$2, $3}' ${PHENO} > ${META_DIR}/pheno_fid_iid.txt

# 4. Sort for Joining
# We sort the processed covariates by IID (column 2)
sort -k2,2 ${META_DIR}/covars_fid_iid.txt > ${META_DIR}/cov1.tmp
sort -k2,2 ${META_DIR}/covars2_fid_iid.txt > ${META_DIR}/cov2.tmp
# We sort the batch map by SampleID (column 1)
sort -k1,1 ${META_DIR}/batch_map.txt > ${META_DIR}/batch.tmp

# 5. Merge Batch ID into Covariates
# Logic: join batch.tmp (key $1) with cov.tmp (key $2)
# Warning: This requires the VCF Sample ID ($1 in batch.tmp) to match the IID constructed in cov1.tmp ($2)
awk 'NR==FNR {gsub(/\r$/, "", $2); batch[$1]=$2; next} 
     {gsub(/\r$/, ""); print $1, $2, $3, $4, batch[$2]}' \
     ${META_DIR}/batch.tmp ${META_DIR}/cov1.tmp > ${META_DIR}/covars_with_batch.txt

awk 'NR==FNR {gsub(/\r$/, "", $2); batch[$1]=$2; next} 
     {gsub(/\r$/, ""); print $1, $2, $3, $4, $5, batch[$2]}' \
     ${META_DIR}/batch.tmp ${META_DIR}/cov2.tmp > ${META_DIR}/covars2_with_batch.txt

# Set new variable paths for final analysis
FINAL_COV1="${META_DIR}/covars_with_batch.txt"
FINAL_COV2="${META_DIR}/covars2_with_batch.txt"
FINAL_PHENO="${META_DIR}/pheno_fid_iid.txt"

# ==============================================================================
# MERGING CHROMOSOMES & QC
# ==============================================================================

mkdir -p ${OUT_FOLDER}/all_chrs
IN="${OUT_FOLDER}/bed"
PREFIX="${OUT_FOLDER}/all_chrs/all_chrs"
MERGE_LIST="${OUT_FOLDER}/merge_list.txt"

rm -f ${MERGE_LIST}
for i in {2..22}; do
  echo "${IN}/chr${i}.bed ${IN}/chr${i}.bim ${IN}/chr${i}.fam" >> ${MERGE_LIST}
done

echo "Merging all chromosomes..."
plink \
  --bfile ${IN}/chr1 \
  --merge-list ${MERGE_LIST} \
  --make-bed \
  --out ${PREFIX}

echo "Running QC..."
# QC
plink \
  --bfile ${PREFIX} \
  --mind 0.05 \
  --geno 0.05 \
  --maf 0.01 \
  --hwe 1e-6 \
  --autosome \
  --make-bed \
  --out ${PREFIX}_clean

# ==============================================================================
# FINAL ANALYSES
# ==============================================================================

echo "Running model 1..."
plink \
  --bfile ${PREFIX}_clean \
  --covar ${FINAL_COV1} \
  --pheno ${FINAL_PHENO} \
  --linear no-x-sex hide-covar \
  --all-pheno \
  --ci 0.95 \
  --xchr-model 0 \
  --allow-no-sex \
  --out ${PREFIX}_model1_results

echo "Running model 2..."
plink \
  --bfile ${PREFIX}_clean \
  --covar ${FINAL_COV2} \
  --pheno ${FINAL_PHENO} \
  --linear no-x-sex hide-covar \
  --all-pheno \
  --ci 0.95 \
  --xchr-model 0 \
  --allow-no-sex \
  --out ${PREFIX}_model2_results

plink \
  --bfile ${PREFIX}_clean \
  --freq \
  --out ${PREFIX}_allele_freq

plink \
  --bfile ${PREFIX}_clean \
  --hardy \
  --out ${PREFIX}_hwe

plink \
  --bfile ${PREFIX}_clean \
  --missing \
  --out ${PREFIX}_callrate

echo "Done."
