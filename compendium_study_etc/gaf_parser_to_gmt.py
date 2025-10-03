import pandas as pd
import os 

gaf = pd.read_csv("functionome_release.gaf", sep="\t", comment="!", header=None)
go_terms = gaf.groupby(4)[2].apply(list)  

current_dir = os.path.dirname(__file__)
file_path = os.path.join(current_dir, "go_terms.gmt")

with open(file_path, "w") as f:
    for term, genes in go_terms.items():
        f.write(f"{term}\t" + "\t".join(genes) + "\n")