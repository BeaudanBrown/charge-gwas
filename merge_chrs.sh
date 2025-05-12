#!/usr/bin/env bash
set -euo pipefail

OUT="qc/all_chrs"

# Step 1: initialize OUT with chr1
echo "Initializing ${OUT} from chr1..."
cp qc/chr1.bed   ${OUT}.bed
cp qc/chr1.bim   ${OUT}.bim
cp qc/chr1.fam   ${OUT}.fam

# Step 2: merge chromosomes 2..22 one at a time
for i in {2..22}; do
  CHR="qc/chr${i}"
  echo "Merging ${CHR} into ${OUT}..."
  plink \
    --bfile ${OUT} \
    --bmerge ${CHR} \
    --memory 25000 \
    --out ${OUT}_tmp

  # on success, replace OUT with the merged version
  mv ${OUT}_tmp.bed ${OUT}.bed
  mv ${OUT}_tmp.bim ${OUT}.bim
  mv ${OUT}_tmp.fam ${OUT}.fam
done

echo "All chromosomes merged into ${OUT}.bed/.bim/.fam"

