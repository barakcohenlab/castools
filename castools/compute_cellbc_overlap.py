import gzip
import sys

"""
 python3 compute_cellbc_overlap.py ../spike-in-090821-10xtranscriptome/SL5-1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz SL4-1_cells_umis_tripbc.tsv
"""

tenx_cells = sys.argv[1]
cas_cells_f = sys.argv[2]

tenx_cellbcs = {}
with gzip.open(tenx_cells, "rb") as fh:
    for line in fh:
        cellbc = line.decode().rstrip("-1\n")
        tenx_cellbcs[cellbc]= 1


cas_cellbcs = {}
n_overlap = 0
with open(cas_cells_f) as fh:
    for line in fh:
        cas_cell = line.rstrip("\n").split()[0]
        if cas_cell not in cas_cellbcs:
            cas_cellbcs[cas_cell] = 1
            if cas_cell in tenx_cellbcs:
                n_overlap += 1
        if cas_cell in tenx_cellbcs:
            print(line.rstrip("\n"))

print(len(tenx_cellbcs), len(cas_cellbcs), n_overlap, file = sys.stderr)
