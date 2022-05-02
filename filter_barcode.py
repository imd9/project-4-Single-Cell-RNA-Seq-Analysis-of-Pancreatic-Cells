# Author: Monica Roberts
# Purpose: Takes barcode list text file, filters for barcodes that have > 1000 reads, and outputs list of filtered barcodes in txt file.

import sys

text_file = sys.argv[1]
sample_name = text_file.split('.')[0].split('/')[-1]

filtered_barcode = []

with open(text_file, 'r') as f:
    for line in f:
        barcode = line.split(':')[0]
        read_count = line.split(':')[1]
        if int(read_count) > 1000:
            filtered_barcode.append(barcode)

with open('filtered_barcode_list_' + sample_name, 'w') as f:
    for item in filtered_barcode:
        f.write(item + '\n')
