import gzip
import sys
import argparse
from fuzzysearch import find_near_matches
#from fuzzywuzzy import fuzz
#from fuzzywuzzy import process

def parse_arguments():
    parser = argparse.ArgumentParser(description='Print the cell-barcode, UMI and trip barcode from the nextseq run. Use the 10x whitelist to whitelist the cell-barcodes. Only print the white-listed cell-barcodes.')
    parser.add_argument("R1",
                        help = "Read1 of the FASTQ pair")
    parser.add_argument("R2",
                        help = "Read2 of the FASTQ pair")
    args = parser.parse_args()
    return args

## From https://codereview.stackexchange.com/questions/151329/reverse-complement-of-a-dna-string
def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in dna[::-1]])

 # From http://code.activestate.com/recipes/499304-hamming-distance/
def hamdist(str1, str2):
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

def read_whitelist():
    #Read in the 10X v3 whitelist
    v3_whitelist = {}
    with gzip.open("../dat/10xv3_whitelist.txt.gz", "rt") as whitelist_fh:
        for bc in whitelist_fh:
            v3_whitelist[bc.rstrip("\n")] = 1
    return v3_whitelist

def translate_10x_whitelist():
    translate_table = {}
    with open("../dat/3M-february-2018-translatetable.txt") as translate_table_f:
        for line in translate_table_f:
            fields = line.split()
            translate_table[fields[1]] = fields[0]
    return translate_table

def parse_fastq(args, v3_whitelist):
    cellumitrip = {} #key is cellbc+umi+tripbc, value is number of reads
    #Prints the TRIP barcode from each read
    usable_reads_n = 0
    not_usable_reads_n = 0
    total_reads_n = 0
    line_num = 0
    polya_n = 0
    captureseq_n = 0
    polya_notvalid_n = 0
    other_notvalid_n = 0
    fuzzy_cutoff_ratio = 0.85 # 0.85 * 16 = 13.6, so 14, 15, 16
    with gzip.open(args.R1, 'rt') as r1:
        with gzip.open(args.R2, 'rt') as r2:
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
                    after_trip_bottom = "TCTAGA"
                    after_trip_top = "CTCGAG"
                    #print(fuzz.partial_ratio(line2, end_of_bfp))
                    #print(fuzz.partial_ratio(line1, before_trip))
                    #print(find_near_matches(end_of_bfp, line2, max_l_dist=2))
                    m1 = find_near_matches(end_of_bfp, line2, max_l_dist=2)
                    m2 = find_near_matches(before_trip, line1, max_l_dist=2)
                    #print(process.extractOne(end_of_bfp, line2, scorer=fuzz.token_sort_ratio))
                    if end_of_bfp in line2 or m1: #Try to get tripBC using Read2, look for end of BFP in read2
                        #print("m1", m1[0].start, m1, end_of_bfp, line2, file = s)
                        beg_pos = m1[0].start
                        print(str(len(line2) - beg_pos) + " pA" + " " + line1[(beg_pos) + 16 :beg_pos + 32] + " " + line2[beg_pos + 32: beg_pos + 38], file = sys.stderr)
                        # Get umi and cell bc from read 1
                        if len(line2) - beg_pos >= 38 and line2[beg_pos + 32: beg_pos + 38] == after_trip_top:
                            cell_bc = line1[0:16]
                            umi = line1[16:28]
                            trip_bc = line2[beg_pos + 16 :beg_pos + 32]
                            # Should we do the error correction first before filtering through the whitelist?
                            #if cell_bc in v3_whitelist:
                            uid = cell_bc + " " + umi + " " + trip_bc
                            if uid not in cellumitrip:
                                cellumitrip[uid] = 0
                            cellumitrip[uid] += 1
                            polya_n += 1
                            usable_reads_n += 1
                            #else:
                            #    polya_notvalid_n += 1
                    elif before_trip in line1 or m2: #Try to get tripBC using just read1, look for capture sequence
                        beg_pos = m2[0].start
                        print(str(len(line1) - beg_pos) + "\tCS" + " " + line1[(beg_pos) + 16 :beg_pos + 32] +  " " + line1[beg_pos + 32: beg_pos + 38], file = sys.stderr)
                        if len(line1) - beg_pos >= 38 and line1[beg_pos + 32: beg_pos + 38] == after_trip_bottom:
                            cell_bc = line1[0:16]
                            umi = line1[16:28]
                            trip_bc = line1[(beg_pos) + 16 :beg_pos + 32]
                            # Should we do the error correction first before filtering through the whitelist?
                            #if cell_bc in v3_whitelist:
                            uid = cell_bc + " " + umi + " "  + trip_bc
                            if uid not in cellumitrip:
                                cellumitrip[uid] = 0
                            cellumitrip[uid] += 1
                            captureseq_n += 1
                            usable_reads_n += 1
                            #else:
                            #    other_notvalid_n += 1
                            #    print("not in whitelist", cell_bc, file = sys.stderr)
                    else:
                        not_usable_reads_n += 1
                        #print("not-usable-read", line1, line2, file = sys.stderr)
    total_notusable = not_usable_reads_n + polya_notvalid_n + other_notvalid_n
    print(args.R1, file = sys.stderr)
    print("Total, usable, usable-fraction, polya-reads, captureseq-reads, total_not_usable, polya_notusable, withcapture_notusable, other_notusable \n", total_reads_n, usable_reads_n, usable_reads_n/(total_reads_n), polya_n, captureseq_n,  total_notusable, not_usable_reads_n, polya_notvalid_n, other_notvalid_n, file = sys.stderr)
    return cellumitrip

def main():
    v3_whitelist = read_whitelist()
    #translate_table = translate_10x_whitelist()
    ########
    # How did you translate capture sequence bcs to polyA BCs?
    #######
    args = parse_arguments()
    cellumitrip = parse_fastq(args, v3_whitelist)
    for uid in cellumitrip:
        print(uid + " " + str(cellumitrip[uid]))

if __name__ == "__main__":
    main()
