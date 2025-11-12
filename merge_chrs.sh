#!/usr/bin/env bash
set -euo pipefail

OUT="bed/all_chrs"

# Step 1: initialize OUT with chr1
echo "Initializing ${OUT} from chr1..."
cp bed/chr1.bed   ${OUT}.bed
cp bed/chr1.bim   ${OUT}.bim
cp bed/chr1.fam   ${OUT}.fam

CHRS=$(echo {2..22})
# Step 2: merge chromosomes 2..22 one at a time
for i in {2..22}; do
  CHR="bed/chr${i}"
  echo "Merging ${CHR} into ${OUT}..."
  plink \
    --bfile ${OUT} \
    --bmerge ${CHR} \
    --out ${OUT}_tmp

  # on success, replace OUT with the merged version
  mv ${OUT}_tmp.bed ${OUT}.bed
  mv ${OUT}_tmp.bim ${OUT}.bim
  mv ${OUT}_tmp.fam ${OUT}.fam
done

echo "All chromosomes merged into ${OUT}.bed/.bim/.fam"

