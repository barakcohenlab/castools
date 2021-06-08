import gzip
import sys
import argparse
import csv

def parse_arguments():
    parser = argparse.ArgumentParser(description="collapse cell-barcodes using translate table")
    parser.add_argument(
        'trio',
        help="Path to the trip barcode corrected TRIOs (with header)"
    )
    args = parser.parse_args()
    return args

def read_translate_table():
    translate_table = {}
    with gzip.open("../dat/3M-february-2018-translatetable.txt.gz", "rt") as translate_table_f:
        for line in translate_table_f:
            fields = line.split()
            translate_table[fields[1]] = fields[0]
    return translate_table

def translate_cellbarcodes(trio_file, translate_table):
    with open(trio_file) as trio_fh:
        reader = csv.DictReader(trio_fh, delimiter = "\t")
        print("\t".join(reader.fieldnames))
        for line in reader:
            try:
                if line['cellBC'] in translate_table:
                    line['cellBC'] = min(line['cellBC'], translate_table[line['cellBC']])
            except KeyError:
                print("Please make sure trios file has the right header.", file = sys.stderr)
                sys.exit(1)
            print("\t".join(line.values()))

def main():
    args = parse_arguments()
    translate_table = read_translate_table()
    translate_cellbarcodes(args.trio, translate_table)

if __name__ == "__main__":
    main()
