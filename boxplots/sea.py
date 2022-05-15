import json

import matplotlib.pyplot as plt
import seaborn as sns


def read_data(path):
    with open(path, 'r') as f:
        data = json.loads(f.readlines()[0])
    return data

def draw_boxplot(ax, title, data, tool_blacklist, data_label, colors):
    ax.set_title(title)
    ax.tick_params(labelrotation=60, axis='x')
    all_plot_data = []
    all_labels = []
    for i in range(len(data)):
        plot_data = []
        labels = []
        for tool in data[i]:
            if tool in tool_blacklist:
                continue
            labels.append(tool)
            plot_data.append(data[i][tool][data_label])
        all_plot_data.append(plot_data)
        all_labels.append(labels)
    final_plot_data = []
    final_labels = []
    final_colors = []
    final_vertical_lines = []

    total_i = 0
    for plot_data, label_data, color in zip(all_plot_data, all_labels, colors):
        for val, label, i in zip(plot_data, label_data, range(len(label_data))):
            total_i += 1
            final_plot_data.append(val)
            final_labels.append(label)
            final_colors.append(sns.xkcd_rgb[color])
            if i == len(label_data)-1:
                final_vertical_lines.append(total_i-0.5)
    g = sns.boxplot(ax=ax,
                    data=final_plot_data,
                    palette=final_colors,
                    showmeans=False,
                    meanprops={"marker":"s",
                               "markerfacecolor":"white",
                               "markeredgecolor":"blue"})
    # remove last line to prevent thick line on the right
    final_vertical_lines.pop()
    for index in final_vertical_lines:
        plt.axvline(x=index, color="black")
    g.set_xticklabels(final_labels)
    plt.setp(ax.xaxis.get_majorticklabels(), ha='right')

# TP, FP, TN, FN, P, N, TPR, TNR, PPV, ACC, BA, runtime
# Accuracy = ACC
# Precision = PPV
# Sensitivity = TPR
# Runtime = runtime
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10,5), constrained_layout=True)
data = 'oral'
hiseq = read_data(f'hiseq_{data}_microbiome_paired_parsed.txt')
miseq = read_data(f'miseq_{data}_microbiome_paired_parsed.txt')
nanopore = read_data(f'nanopore_{data}_microbiome_unpaired_parsed.txt')
draw_boxplot(axes,
             'Simulated oral microbiome sensitivity',
             [hiseq, miseq, nanopore],
             [],#['BBMap_default', 'BWA_MEM2', 'Bowtie2_end_to_end_un_conc', 'Minimap2_nanopore', 'Kraken2Minimap2_nanopore'],
             'TPR',
             ['yellow', 'cyan', 'pale red'])
plt.savefig('fig.png')
plt.show()
