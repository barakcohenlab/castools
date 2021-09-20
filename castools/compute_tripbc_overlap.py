import sys
"""
python3 compute_tripbc_overlap.py CAS_LP1_pool2_DNA_barcodes.txt SL4-1_cells_umis_tripbc.tsv
"""

bulk_trips = sys.argv[1]
cas_trips_f = sys.argv[2]

bulk_tripbcs = {}
with open(bulk_trips) as fh:
    for line in fh:
        tripbc = line.rstrip("\n")
        bulk_tripbcs[tripbc]= 1


cas_tripbcs = {}
n_overlap = 0
with open(cas_trips_f) as fh:
    for line in fh:
        cas_trip = line.rstrip("\n").split()[2]
        if cas_trip not in cas_tripbcs: #ignore duplicate trip barcodes
            cas_tripbcs[cas_trip] = 1
            if cas_trip in bulk_tripbcs:
                n_overlap += 1
        # Print all the lines
        if cas_trip in bulk_tripbcs:
            print(line.rstrip("\n"), "overlap")
        else:
            print(line.rstrip("\n"), "no-overlap")

print(len(bulk_tripbcs), len(cas_tripbcs), n_overlap, file = sys.stderr)
