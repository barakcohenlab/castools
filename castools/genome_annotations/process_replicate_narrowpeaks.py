#!/usr/bin/python3

import sys
import argparse
import os
import gzip as gz

parser = argparse.ArgumentParser(description = 'Remove duplicated narrowpeaks')
parser.add_argument('narrowpeak_file', help = 'Narrowpeaks file downloaded from ENCODE')
args = parser.parse_args()

basename = os.path.basename(args.narrowpeak_file)

non_duplicated_lines = {}

with gz.open(args.narrowpeak_file, 'rt') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        loc = '-'.join(line[:3])
        signal_value = float(line[6])
        if loc in non_duplicated_lines.keys():
            if signal_value > float(non_duplicated_lines[loc][6]):
                non_duplicated_lines[loc] = line
        else:
            non_duplicated_lines[loc] = line

with gz.open(basename + '_no_dup.bed.gz', 'wt') as f:
    for name, line in non_duplicated_lines.items():
        f.write('\t'.join(line) + '\n')