import matplotlib.pyplot as plt
import gzip

relative_positions = []
fragment_lengths = []
input_file = 'mapped.bed.gz'

with gzip.open(input_file, 'rt') as file:
    for line in file:
        if line.startswith('#'):
            continue
        cols = line.split()
        if len(cols) < 12:
            continue
        start = int(cols[2])
        end = int(cols[3])
        frag_start = int(cols[8])
        frag_end = int(cols[9])
        fragment_length = int(cols[11])
        fragment_midpoint = (frag_start + frag_end) / 2
        genomic_midpoint = (start + end) / 2
        relative_position = fragment_midpoint - genomic_midpoint
        if -500 <= relative_position <= 500:
            relative_positions.append(relative_position)
            fragment_lengths.append(fragment_length)

plt.figure(figsize=(50, 30))
plt.hexbin(relative_positions, fragment_lengths, gridsize=200, cmap='Blues', mincnt=1)

plt.xlim(-500, 500)
plt.ylim(0, 500)
plt.xlabel('Position relative to midpoint', fontsize=50)
plt.ylabel('Fragment Length (bp)', fontsize=50)
plt.title('V-Plot: Gene Fragments', fontsize=70)

cb = plt.colorbar()
cb.set_label('Frequency', fontsize=50)
cb.ax.tick_params(labelsize=40)

plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.grid(False)

plt.savefig('v_plot_output.png')
plt.show()
