# Author: Monica Roberts (data curator)
# Purpose: Takes zipped fasta file as input and outputs txt file that contains unique barcodes with number of reads for each.

import sys
import gzip

barcode_dict = {}

file_input = sys.argv[1]
file_name = file_input.split('.')[0]
doc_title = file_name.split('/')[-1]

with gzip.open(file_input, 'rt') as f:
    for line in f:
        if '=' in line:
            barcode1 = line.split('=')[1].split(' ')[0]
            barcode2 = line.split('=')[2].split(' ')[0]
            barcode = barcode1+barcode2
            if barcode in barcode_dict:
                barcode_dict[barcode] += 1
            else:
                barcode_dict[barcode] = 1

barcode_freq = list(barcode_dict.values())

with open('barcode_list_'+ doc_title + '.txt', 'w') as g:
    for key, value in barcode_dict.items(): 
        g.write('%s:%s\n' % (key, value))


