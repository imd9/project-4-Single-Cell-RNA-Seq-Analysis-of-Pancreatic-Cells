# Author: Monica Roberts
# Purpose: Takes text files of filtered barcodes from each sample and combines them into single list in text file for salmon alevin whitelist input.

import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]

files = [file1, file2, file3]
unique_barcodes = []

for file in files:
    with open(file, 'r') as f:
        for line in f:
            barcode = line
            if barcode not in unique_barcodes:
                unique_barcodes.append(barcode)

with open('barcode_whitelist.txt', 'w') as f:
    for item in unique_barcodes:
        f.write(item)


