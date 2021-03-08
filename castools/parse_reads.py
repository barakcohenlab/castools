import gzip
import sys

"""
Print the cell-barcode, UMI and trip barcode from the nextseq run. Use the
10x whitelist to whitelist the cell-barcodes. Only print the white-listed
cell-barcodes.
"""

## From https://codereview.stackexchange.com/questions/151329/reverse-complement-of-a-dna-string
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

def hamdist(str1, str2): # From http://code.activestate.com/recipes/499304-hamming-distance/
    """Count the # of differences between equal length strings str1 and str2"""
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

def complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna])

#Compute minimum hamming distance between a given barcodes and all barcodes in the whitelist
def min_hamming(v3_whitelist, bc1):
    min_hamming = 200
    for whitelist_bc in v3_whitelist:
        hamming = hamdist(bc1, whitelist_bc)
        if hamming < min_hamming:
            min_hamming = hamming
    return(min_hamming)

#Read in the 10X v3 whitelist
v3_whitelist = {}
with open("../dat/10xv3_whitelist.txt") as whitelist_fh:
    for bc in whitelist_fh:
        v3_whitelist[bc.rstrip("\n")] = 1

translate_table = {}
with open("../dat/3M-february-2018-translatetable.txt") as translate_table_f:
    for line in translate_table_f:
        fields = line.split()
        translate_table[fields[1]] = fields[0]

#Prints the TRIP barcode from each read
usable_reads_n = 0
not_usable_reads_n = 0
total_reads_n = 0
line_num = 0
polya_n = 0
captureseq_n = 0
polya_notvalid_n = 0
other_notvalid_n = 0
cellumitrip = {} #key is cellbc+umi+tripbc, value is number of reads
with gzip.open(sys.argv[1], 'rt') as r1:
    with gzip.open(sys.argv[2], 'rt') as r2:
        for line1,line2 in zip(r1,r2):
            line1 = line1.rstrip("\n")
            line2 = line2.rstrip("\n")
            line_num += 1
            if line_num % 4 == 2: #Only the sequence lines
                total_reads_n += 1
                #Check for end of BFP sequence, this is the sequence after the TRIP barcode
                #if "TTGCTAGGAC" in line and "TCTAGACTCGAAGCGAGAGCT" in line:
                capture_seq = "TTGCTAGGAC"
                end_of_bfp = "TCGCTTCGAGTCTAGA"
                before_trip = "CCGGCCACAACTCGAG"
                if end_of_bfp in line2: #Try to get tripBC using Read2, look for end of BFP in read2
                    beg_pos = line2.find(end_of_bfp)
                    # Get umi and cell bc from read 1
                    cell_bc = line1[0:16]
                    umi = line1[16:28]
                    trip_bc = line2[beg_pos + 16 :beg_pos + 32]
                    if cell_bc in v3_whitelist:
                        uid = cell_bc + " " + umi + " " + trip_bc
                        if uid not in cellumitrip:
                            cellumitrip[uid] = 0
                        cellumitrip[uid] += 1
                        #print(cell_bc, umi, trip_bc)
                        polya_n += 1
                        usable_reads_n += 1
                    else:
                        polya_notvalid_n += 1
                        print("Not in whitelist", cell_bc, line1, file = sys.stderr)
                        #if reverse_complement(cell_bc) in v3_whitelist:
                        #    print("RC in whitelist - polyA", cell_bc, line1, file = sys.stderr)
                        #if cell_bc in translate_table:
                        #    print("RC in translate_table - polyA", cell_bc, line1, file = sys.stderr)
                        #else:
                        #    print("RC not in translate_table - polyA", cell_bc, line1, file = sys.stderr)
                elif before_trip in line1: #Try to get tripBC using just read1, look for capture sequence
                    beg_pos = line1.find(before_trip)
                    if len(line1) - beg_pos >= 31:
                        cell_bc = line1[0:16]
                        umi = line1[16:28]
                        trip_bc = line1[(beg_pos) + 16 :beg_pos + 32]
                        if cell_bc in v3_whitelist:
                            uid = cell_bc + " " + umi + " "  + trip_bc
                            if uid not in cellumitrip:
                                cellumitrip[uid] = 0
                            cellumitrip[uid] += 1
                            #print(cell_bc, umi, trip_bc)
                            captureseq_n += 1
                            usable_reads_n += 1
                        else:
                            other_notvalid_n += 1
                            print("not in whitelist", cell_bc, file = sys.stderr)
                else:
                    not_usable_reads_n += 1
                    #print("not-usable-read", line1, line2, file = sys.stderr)
total_notusable = not_usable_reads_n + polya_notvalid_n + other_notvalid_n
for uid in cellumitrip:
    print(uid + " " + str(cellumitrip[uid]))
print(sys.argv[1], file = sys.stderr)
print("Total, usable, usable-fraction, polya-reads, captureseq-reads, total_not_usable, polya_notusable, withcapture_notusable, other_notusable \n", total_reads_n, usable_reads_n, usable_reads_n/(total_reads_n), polya_n, captureseq_n,  total_notusable, not_usable_reads_n, polya_notvalid_n, other_notvalid_n, file = sys.stderr)
