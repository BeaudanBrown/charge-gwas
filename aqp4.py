from redcap import project
from dotenv import load_dotenv
import os
import pandas as pd
from io import StringIO


load_dotenv()

api_url = 'https://redcap.helix.monash.edu/api/'
api_key = os.environ['REDCAP_API_KEY']
proj = project.Project(api_url, api_key)

# bcftools view -R apoe_snps.txt ~/AGRF_CAGRF25020097_customGenotyping/ILGSA24-00261/ILGSA24-00261_chr19_filtered.vcf.gz -Ov -o apoe.vcf

# bcftools view -R aqp4_snps.txt ~/AGRF_CAGRF25020097_customGenotyping/ILGSA24-00261/ILGSA24-00261_chr18_filtered.vcf.gz -Ov -o aqp4.vcf
def get_results(filename):
    header = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"):
                header = line.split("\t")

    def clean_header(col):
        col = col.strip("\n")
        col = col.split("_", 1)[1] if "_" in col else col
        col = col.strip("#")
        return col

    headers = list((clean_header(col) for col in header))

    with open(filename) as f:
        lines = [line for line in f if not line.startswith("#")]
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", header=None)
    df.columns = headers

    vcf_columns = [
        "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"
    ]
    sample_cols = df.columns[len(vcf_columns):]
    df[sample_cols]

    melted = df.melt(id_vars=["ID"], value_vars=sample_cols,
                     var_name="sample", value_name="genotype")
    pivoted = melted.pivot(index="sample", columns="ID", values="genotype")

    pivoted_df = pivoted.reset_index().rename_axis(None, axis=1)

    long = pivoted_df.melt(id_vars="sample",
                   var_name="snp_id",
                   value_name="gt_prob")

    long[["gt", "dosage"]] = long["gt_prob"].str.split(":", expand=True)
    long["dosage"] = long["dosage"].astype(float)
    long["ref"] = long["snp_id"].str.split(":", expand=True)[2]
    long["alt"] = long["snp_id"].str.split(":", expand=True)[3]

    # def genotype_to_alleles(gt, ref, alt):
    #     return ''.join([ref if allele == '0' else alt for allele in gt.split('|')])

    def genotype_to_alleles(gt):
        alleles = gt.split('|')
        if alleles == ['0', '0']:
            return "1"
        elif alleles == ['2', '2']:
            return "3"
        else:
            return "2"

    long['alleles'] = long.apply(lambda row: genotype_to_alleles(row['gt']), axis=1)
    long["alleles"] = long["alleles"]

    long["pos"] = long["snp_id"].str.split(":", expand=True)[1]

    wide = (
        long
        .pivot(index="sample",
               columns="pos",
               values=["alleles", "dosage"])
    )

    wide.columns = [
        f"{pos}_{field}"
        for field, pos in wide.columns
    ]
    wide = wide.reset_index()

    return wide

aqp4 = get_results("aqp4.vcf").set_index("sample")
apoe = get_results("apoe.vcf").set_index("sample")
all_snps = aqp4.join(apoe)

all_snps = all_snps.rename(columns={
    '24439072_alleles': 'aqp4_allele1',
    '24439072_dosage': 'aqp4_dosage1',
    '24435587_alleles': 'aqp4_allele2',
    '24435587_dosage': 'aqp4_dosage2',
    '24431689_alleles': 'aqp4_allele3',
    '24431689_dosage': 'aqp4_dosage3',
    '45412079_alleles': 'apoe_allele1',
    '45412079_dosage': 'apoe_dosage1',
    '45411941_alleles': 'apoe_allele2',
    '45411941_dosage': 'apoe_dosage2',
})

records = []
for sample_id, row in all_snps.iterrows():
    row_dict = {
        'idno': f"{str(sample_id).strip("BACH")}--1",
        "genomics_complete": "2",
        "redcap_event_name": "baseline_arm_1",
        "aqp4_allele1": row["aqp4_allele1"],
        "aqp4_allele2": row["aqp4_allele2"],
        "aqp4_allele3": row["aqp4_allele3"],
        "aqp4_dosage1": float(row["aqp4_dosage1"]),
        "aqp4_dosage2": float(row["aqp4_dosage2"]),
        "aqp4_dosage3": float(row["aqp4_dosage3"]),

        "apoe_allele1": row["apoe_allele1"],
        "apoe_allele2": row["apoe_allele2"],
        "apoe_dosage1": float(row["apoe_dosage1"]),
        "apoe_dosage2": float(row["apoe_dosage2"])
    }
    records.append(row_dict)

proj.import_records(records)
