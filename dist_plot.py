import matplotlib.pyplot as plt
plt.switch_backend('agg')
import sys

file_name = sys.argv[1]
sample_name = file_name.split('.')[0].split('/')[-1]

counts = []
with open(file_name, 'r') as f:
    for line in f:
        read_count = line.split(':')[1]
        counts.append(read_count)

counts.sort()

plt.yscale('log')
plt.hist(counts, cumulative=True, color='green')
plt.savefig(sample_name + '.png')