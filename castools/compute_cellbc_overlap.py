import gzip
import sys
import pandas as pd
"""
 python3 compute_cellbc_overlap.py ../spike-in-090821-10xtranscriptome/SL5-1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz SL4-1_cells_umis_tripbc.tsv
"""

tenx_cells = sys.argv[1]
cellBC = pd.read_csv(sys.argv[2], sep = '\t' , header = ['cellBC', 'umi', 'tBC', 'count'])
tenx_cellbcs = {}
with gzip.open(tenx_cells, "rb") as fh:
    for line in fh:
        cellbc = line.decode().rstrip("-1\n")
        tenx_cellbcs[cellbc]= 1

def hammingDist(str1, str2):
    '''
    Calculating hamming distance of two strings
    https://www.geeksforgeeks.org/hamming-distance-two-strings/
    '''
    i = 0
    count = 0
    while(i < len(str1)):
        if(str1[i] != str2[i]):
            count += 1
        i += 1
    return count

def ec_cellBC(cellBC, cr_cellBC_list):
    pop_list = []
    for cr_cellBC in cr_cellBC_list:
        hamming = hammingDist(cellBC, cr_cellBC)
        if hamming <= 3:
            pop_list.append([cr_cellBC, hamming])
    if len(pop_list) == 1:
        return pop_list[0][0], pop_list[0][1]
    else:
        #print(f'{cellBC} is a bad cell!')
        return 0, 16

def cell_bc_ec(quad, cr_cellBC_list):
    '''
    Input: quad: tsv file with columns: cellBC, UMI, pBC, rBC, and read count
        cellBC: a tsv file with cellranger corrected cell barcodes.
    Output: cellranger error corrected quad file with the same structure.
    '''
    counter = 0
    barcode_mapping_dict = {}
    # First extract the quad_cellBC set from quad file
    quad_cellBC = list(set(quad['cellBC'].values))
    # Then error correct the quad_cellBCs if it is in the cellBC list or within
    # If the cellBC is the exact match to the cr_cellBC, then the cellBC is kept
    for cellBC in quad_cellBC:
        if cellBC in cr_cellBC_list:
            barcode_mapping_dict[cellBC] = cellBC
        else:
            close_mapping, hamming = ec_cellBC(cellBC, cr_cellBC_list)
            if hamming <= 3:
                if close_mapping != 0:
                    barcode_mapping_dict[cellBC] = close_mapping
    # Then we go through all the list of cellBCs in the quad file and only record 
    # the ones whose cellBC is in the barcode_mapping_dict!
    pop_list = []
    for _, line in quad.iterrows(): 
        cellBC = line['cellBC']
        if cellBC in barcode_mapping_dict:
            # Here I just printed out the error corrected files
            counter += 1
            print(barcode_mapping_dict[cellBC], line['umi'], line['tBC'], line['count'])
    return counter
    
num_cells = cell_bc_ec(cellBC, list(tenx_cellbcs.keys()))


print(len(tenx_cellbcs), num_cells, file = sys.stderr)
