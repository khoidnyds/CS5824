from pyfaidx import Fasta
import pandas as pd
import subprocess
from pathlib import Path
from glob import glob

for i in range(1, 5):
    subprocess.run(
        ["wget", f"https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.{i}.1.genomic.fna.gz"])
    subprocess.run(['gunzip', '-d', f'viral.{i}.1.genomic.fna.gz'])

subprocess.run(["cat *.fna > full_viral.fna"], shell=True)
path = Path("raw_data")
Path.mkdir(path, parents=True, exist_ok=True)
subprocess.run(['mv', 'full_viral.fna', str(path)])
subprocess.run(["rm *.fna"], shell=True)

accession = []
description = []
v1 = Fasta("raw_data/full_viral.fna")
for v in v1:
    l = v.long_name.split(" ")
    a = l[0]
    l = l[:-2]
    d = " ".join(l)[:-1]
    accession.append(a)
    description.append(d)
metadata = pd.DataFrame({"Accession": accession, "Description": description})
metadata.to_csv("raw_data/annotation.txt", index=False)
