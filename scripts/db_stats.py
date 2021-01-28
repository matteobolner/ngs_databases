import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import os
import json
import numpy as np

get_stats(9822, '/home/pelmo/db_final/9822_data/')
get_stats(9984, '/home/pelmo/db_final/9984_data/')
get_stats(9913, '/home/pelmo/db_final/9913_data/')



def get_stats(animal_taxon, path_to_folder):
    #animal_taxon = 9822
    genome_sizes = {'9822':'2501912388', '9984':'2737490501', '9913':'2715853792'}
    #path_to_folder = '/home/pelmo/db_final/9822_data/'
    path_to_db = os.path.join(path_to_folder, str(animal_taxon)+'_database.csv')
    path_to_db
    df = pd.read_csv(path_to_db)
    df = df.drop(columns = ['Unnamed: 0'])
    df = df.replace('', 'NaN')

    samples_number = len(df['sample_accession'].unique().tolist())
    samples_number
    projects_number = len(df['study_accession'].unique().tolist())
    projects_number
    samples_with_sex = len(df[['sample_accession','sex']].drop_duplicates()['sex'].dropna())
    #df[['sample_accession','sex']].drop_duplicates()['sex'].dropna().value_counts()
    instrument_platform = len(df[['sample_accession','instrument_platform']].drop_duplicates()['sample_accession'].dropna())
    instrument_platform
    end_type_samples = (df[['sample_accession','library_layout']].drop_duplicates()['library_layout'].dropna())
    paired_end_samples = end_type_samples.value_counts()['PAIRED']
    paired_end_samples
    single_end_samples = end_type_samples.value_counts()['SINGLE']
    single_end_samples
    library_strategy = len(df['library_strategy'].dropna())
    library_strategy
    tissue_biosamples = len(df[['sample_accession','tissue_biosamples']].drop_duplicates()['tissue_biosamples'].dropna())
    tissue_biosamples
    breed_samples = len(df[['sample_accession','breed']].drop_duplicates()['breed'].dropna())
    #tissue_ena = len(df[df['tissue_type_ENA']!='NaN'])
    general_stats = {'samples_number': samples_number, 'projects_number':projects_number, 'samples_with_sex':samples_with_sex, 'samples_with_breed': breed_samples, 'instrument_platform':instrument_platform, 'paired_end_samples':paired_end_samples, 'single_end_samples':single_end_samples, 'library_strategy':library_strategy, 'tissue_biosamples':tissue_biosamples}
    general_stats
    df['base_count']=df['base_count'].fillna(0)

    df['depth'] = df['base_count'].astype(int).div(int(genome_sizes[str(animal_taxon)])).round(2)
    deep_df = df[df['depth'] >=5]

    df['base_count']=df['base_count'].replace(0, 'NONE')
    df['depth']=df['depth'].replace(0, 'NONE')

    #select only wgs, illumina with depth >=5

    deep_df = deep_df[deep_df['instrument_platform']=='ILLUMINA']
    deep_df = deep_df[deep_df['library_source']=='GENOMIC']
    samples_number_w = len(deep_df['sample_accession'].unique().tolist())
    samples_number_w
    projects_number_w = len(deep_df['study_accession'].unique().tolist())
    projects_number_w
    breeds_w = len(deep_df[['sample_accession','breed']].drop_duplicates()['breed'].dropna())
    tissue_w = len(deep_df[['sample_accession','tissue_biosamples']].drop_duplicates()['tissue_biosamples'].dropna())

    breeds_w
    median_depth = deep_df['depth'].median()
    median_depth
    mean_depth = deep_df['depth'].mean()
    mean_depth
    depth_sd = deep_df['depth'].std()
    depth_sd
    min_depth = deep_df['depth'].min()
    min_depth
    max_depth = deep_df['depth'].max()
    max_depth

    wgs_illumina_stats = {'samples_number' : samples_number_w, 'projects_number': projects_number_w, 'breeds':breeds_w, 'tissue': tissue_w,'median_depth':median_depth, 'mean_depth':mean_depth, 'depth_sd':depth_sd, 'min_depth': min_depth, 'max_depth':max_depth }
    general_stats
    wgs_illumina_stats

    library_source_count = df['library_source'].value_counts().to_dict()
    library_source_count
    library_selection_count = df['library_selection'].value_counts().to_dict()
    library_selection_count
    library_strategy_count = df['library_strategy'].value_counts().to_dict()
    library_strategy_count
    library_stats = {'library_source':library_source_count, 'library_selection':library_selection_count, 'library_strategy':library_strategy_count}
    library_stats

    breed_stats = df['breed'].value_counts().to_dict()

    breed_stats_df = pd.DataFrame.from_records(breed_stats, index =[animal_taxon]).transpose()
    general_stats_df = pd.DataFrame.from_records(general_stats, index =[animal_taxon]).transpose()
    wgs_illumina_stats_df = pd.DataFrame.from_records(wgs_illumina_stats, index =[animal_taxon]).transpose()




    library_source_df = pd.DataFrame.from_records(library_source_count, index =[animal_taxon]).transpose()
    library_strategy_df = pd.DataFrame.from_records(library_strategy_count, index =[animal_taxon]).transpose()
    library_selection_df = pd.DataFrame.from_records(library_selection_count, index =[animal_taxon]).transpose()

    stats_to_save = general_stats_df, wgs_illumina_stats_df, breed_stats_df, library_source_df, library_strategy_df, library_selection_df, breed_stats_df
    stats_to_save_names = ['general_stats_df','wgs_illumina_stats_df','breed_stats_df','library_source_df','library_strategy_df','library_selection_df', 'breed_stats_df' ]

    directory = str(animal_taxon) + '_stats'
    stats_path = os.path.join(path_to_folder, directory)
    os.mkdir(stats_path)
    for name, stat in zip(stats_to_save_names, stats_to_save):
        name = name.rstrip('_df') + '_' + str(animal_taxon) + ".csv"
        df_path = os.path.join(stats_path, name)
        stat.to_csv(df_path)

        #df.to_csv(os.path.join(path, str(animal_taxon)+'_database.csv'))
    return()
'''
multiple_samples = df['sample_accession'].value_counts()
multiple_samples
multiple_samples = df[df['sample_accession'].isin(multiple_samples.index[multiple_samples.gt(1)])][['study_accession','sample_accession']]
multiple_samples
multiple_samples['count'] = multiple_samples.groupby(['sample_accession'])['study_accession'].transform('count')
multiple_samples
multiple_samples = multiple_samples.drop_duplicates().sort_values(by= 'count', ascending=False).reset_index(drop=True)
multiple_samples



multiple_runs = df['run_accession'].value_counts()
multiple_runs
multiple_runs = multiple_runs[multiple_runs.gt(1)]
multiple_runs

multiple_samples = df['sample_accession'].value_counts()
multiple_samples
multiple_samples = df[df['sample_accession'].isin(multiple_samples.index[multiple_samples.gt(1)])][['study_accession','sample_accession']]
multiple_samples
multiple_samples['count'] = multiple_samples.groupby(['sample_accession'])['study_accession'].transform('count')
multiple_samples
multiple_samples = multiple_samples.drop_duplicates().sort_values(by= 'count', ascending=False).reset_index(drop=True)

multiple_samples
subspecies = df['scientific_name'].value_counts()

#sample_plot = sns.barplot(data = multiple_samples, x = 'sample_accession', y = 'count', palette = 'muted')
#plt.xticks(rotation=90)

#plt.show()

grouped_df = df.groupby(df.library_source)
library_sources = [name for name,unused_df in grouped_df]
for source in library_sources:
    grouped_df.get_group(source).to_csv(os.path.join(path, str(animal_taxon) + "_" + str(source) + "_database.csv"))
    pd.value_counts(grouped_df.get_group(source)['library_strategy']).to_csv(os.path.join(path,  str(animal_taxon) + "_" + str(source) + "_library_strategies.csv"))

df.to_csv(os.path.join(path, str(animal_taxon)+'_database.csv'))
#multiple_runs.to_csv(os.path.join(path, str(animal_taxon) + '_multiple_runs.csv'))
multiple_samples.to_csv(os.path.join(path, str(animal_taxon) + '_multiple_samples.csv'))
subspecies.to_csv(os.path.join(path, str(animal_taxon) + '_subspecies.csv'))

#plot the library information

sns.countplot(x = df['library_strategy'], palette= 'muted')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(path, str(animal_taxon)+'_library_strategy.png'))

sns.countplot(x = df['library_source'], palette= 'muted')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(os.path.join(path, str(animal_taxon)+'_library_source.png'))
'''
