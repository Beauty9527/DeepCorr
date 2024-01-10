"""
plt.plot, draw the legend with max-min normalization
show the results in the same scale
Proovread can not be run in 3 and 6, expected in above two figures
"""
import numpy as np
import matplotlib.pyplot as plt
methods = ['LoRDEC-29', 'Jabba', 'ColorMap', 'Proovread', 'Hercules', 'DeepCorr']

#  Ecoli pacbio
total_bases = np.array([720464739, 611947598, 730726602, 537183316, 742217998, 748417338])
aligned_bases = np.array([702509360, 611947598, 715441895, 537132179, 723994674, 748015166])
alignment_identity = np.array([0.9750, 1.0000, 0.9790, 0.9999, 0.9754, 0.9995])
max_length = np.array([44078, 41342, 44113, 39836, 44113, 44113])
average_length = np.array([8435, 7880, 8529, 4971, 8691, 8768])
n50 = np.array([13494, 12352, 13641, 9435, 13887, 14020])
genome_fraction = np.array([100.00, 99.25, 100.00, 100.00, 100.00, 100.00])

# min-max normalization
normalized_total_bases = (total_bases - np.min(total_bases)) / (np.max(total_bases) - np.min(total_bases))
normalized_aligned_bases = (aligned_bases - np.min(aligned_bases)) / (np.max(aligned_bases) - np.min(aligned_bases))
normalized_alignment_identity = (alignment_identity - np.min(alignment_identity)) / (np.max(alignment_identity) - np.min(alignment_identity))
normalized_max_length = (max_length - np.min(max_length)) / (np.max(max_length) - np.min(max_length))
normalized_average_length = (average_length - np.min(average_length)) / (np.max(average_length) - np.min(average_length))
normalized_n50 = (n50 - np.min(n50)) / (np.max(n50) - np.min(n50))
normalized_genome_fraction = (genome_fraction - np.min(genome_fraction)) / (np.max(genome_fraction) - np.min(genome_fraction))

# create a figure
plt.figure(figsize=(9, 6))

plt.plot(methods, normalized_total_bases, label='Total Bases', marker='o', markersize=8,  linestyle=':', color='blue')
plt.plot(methods, normalized_aligned_bases, label='Aligned Bases', marker='^', markersize=8, linestyle=':', color='green')
plt.plot(methods, normalized_alignment_identity, label='Alignment Identity', marker='s',  markersize=8, linestyle=':', color='orange')
plt.plot(methods, normalized_max_length, label='Max Length', marker='D', markersize=8, linestyle=':', color='red')
plt.plot(methods, normalized_average_length, label='Average Length', marker='P', markersize=8, linestyle=':', color='purple')
plt.plot(methods, normalized_n50, label='N50', marker='*', markersize=8, linestyle=':', color='brown')
plt.plot(methods, normalized_genome_fraction, label='Genome Fraction', marker='X', markersize=8, linestyle=':', color='pink')

# title and markers
plt.title('Normalized Mapping Statistics for Different Methods on Ecoli PacBio', fontsize=16)
plt.xlabel('Methods', fontsize=18)
plt.ylabel('Normalized Values', fontsize=16)
plt.xticks(rotation=20, ha='right')  # 使x轴标签倾斜才放得下

plt.legend(fontsize=10, loc='lower right')
plt.xticks(fontsize=14)
plt.grid(True)

plt.tight_layout()
plt.show()
