from redcap import project
from dotenv import load_dotenv
import os
import csv

load_dotenv()

api_url = 'https://redcap.helix.monash.edu/api/'
api_key = os.environ['REDCAP_API_KEY']
gwas_folder = os.environ['GWAS_FOLDER']

proj = project.Project(api_url, api_key)
records = proj.export_records(fields=["idno", "sex", "age", "education", "vrii_total_raw"], events=["baseline_arm_1"])
records = [{
    "idno": rec["idno"].split("--")[0],
    "sex": rec["sex"],
    "age": rec["age"],
    "education": rec["education"],
    "vrii_total_raw": rec["vrii_total_raw"],
    } for rec in records if rec["idno"].endswith("--1")]

fam_file = gwas_folder + "/chr1.fam"

with open(fam_file, 'r', newline="") as file:
   # Read the file line by line
   lines = file.readlines()

# Split each line by the tab character
data = [line.strip().split(' ') for line in lines]

final_records = {}
for r in records:
    idno = r['idno']
    # Find the corresponding row in the 2D list
    for row in data:
        if idno in row[1]:
            final_records[idno] = {
                    "sex": r["sex"],
                    "age": r["age"],
                    "education": r["education"],
                    "vrii_total_raw": r["vrii_total_raw"],
                    "FID": row[0],
                    "IID": row[1],
                    }
            break

pheno_filename = "pheno.txt"

with open(pheno_filename, 'w', newline="") as file:
    writer = csv.writer(file, delimiter=' ')
    pheno_data = [
            [r["FID"], r["IID"], r["vrii_total_raw"]] for _, r in final_records.items()
            ]
    pheno_data = sorted(pheno_data, key=lambda x: x[0])
    writer.writerows(pheno_data)

covars_filename = "covars.txt"

with open(covars_filename, 'w', newline="") as file:
    writer = csv.writer(file, delimiter=' ')
    covars_data = [
            [r["FID"], r["IID"], r["sex"], r["age"]] for _, r in final_records.items()
            ]
    covars_data = sorted(covars_data, key=lambda x: x[0])
    writer.writerows(covars_data)

covars2_filename = "covars2.txt"

with open(covars2_filename, 'w', newline="") as file:
    writer = csv.writer(file, delimiter=' ')
    covars_data = [
            [r["FID"], r["IID"], r["sex"], r["age"], r["education"]] for _, r in final_records.items()
            ]
    covars_data = sorted(covars_data, key=lambda x: x[0])
    writer.writerows(covars_data)
