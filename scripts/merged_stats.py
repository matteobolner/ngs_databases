import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import os
import json
import numpy as np

genome_sizes = {'9823':'2501912388', '9986':'2737490501', '9913':'2715853792'}
stats = ['general_stats', 'wgs_illumina_stats', 'library_source', 'library_selection', 'library_strategy']

for stat in stats:
    merge_stats(stat)

def merge_stats(stat):
    path_to_folders = '/home/pelmo/data_and_pipelines/ena_databases/data/updated_databases'
    taxa = ['9823', '9913', '9986']
    names = ['Pig', 'Cattle', 'Rabbit']
    filenames = []
    for taxon in taxa:
        path = os.path.join(path_to_folders, str(taxon) + "_data")
        path_to_folder = os.path.join(path, str(taxon) + "_stats")
        filenames.append(os.path.join(path_to_folder, stat + "_" + str(taxon) + ".csv"))
        stats = pd.concat([pd.read_csv(file, index_col = 0) for file in filenames ], 1)
    stats.columns = names
    stats.to_csv(os.path.join(path_to_folders, stat +'_livestock.csv'))
    return(stats)

def plot_stats(stat, path_to_folders):
    stats_df = merge_stats(stat)
    normalized_stats = stats_df.apply(lambda x: (x/x.sum()*100), axis=0)
    normalized_stats = normalized_stats.reset_index().melt(id_vars='index').sort_values(by=['value'], ascending = False)
    normalized_stats.columns = [library, 'Species','Percentage']
    ##Discard all values lower than 1% for plot readability
    normalized_stats = normalized_stats[normalized_stats['Percentage']>=1]
    print(normalized_stats)
    normalized_stats_plot = sns.catplot(data=normalized_stats, kind = 'bar', x=stat, y='Percentage', hue='Species', height=4, aspect=10/8, hue_order = ['Pig','Cattle','Rabbit'], palette = ["C0", "C1", "C2"] )
    for axes in normalized_stats_plot.axes.flat:
        _ = axes.set_xticklabels(axes.get_xticklabels(), rotation=70)
    sns.set(rc={'figure.figsize':(15,10)})
    plt.savefig(os.path.join(path_to_folders, stat +'.jpg'), bbox_inches = "tight",dpi = 300)

    return(normalized_stats_plot)



for library in ['library_source', 'library_selection', 'library_strategy']:
    plot_stats(library, '/home/pelmo/data_and_pipelines/ena_databases/data/updated_databases')
