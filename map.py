# Author: Monica Roberts
# Purpose: Creates mapping file with transcript ID and gene ID for salmon alevin tool.

import sys
import gzip

reference = sys.argv[1]

mapping = {}

with gzip.open(reference, 'rt') as f:
    for line in f:
        if line[0] == '>':
            transcriptID = line.split('|')[0][1:]
            geneID = line.split('|')[1]
            mapping[transcriptID] = geneID

with open('mapping.tsv', 'w') as g:
    for key, value in mapping.items(): 
        g.write('%s\t%s\n' % (key, value))


    