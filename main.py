import scanpy as sc
import muon as mu
import time


def load_data():
    load_start_time = time.time()
    adata = mu.read("Satija.h5mu")
    print("----Data successfully loaded in %d seconds----" % (time.time() - load_start_time))
    return adata


full_data = load_data()
