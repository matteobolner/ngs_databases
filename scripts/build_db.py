import csv
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import sys
import os
import json
import pickle
import numpy as np


def get_data(animal_taxon):
    #animal_taxon = 9822
    genome_sizes = {'9822':'2501912388', '9984':'2737490501', '9913':'2715853792'}
    headers = {
      'Content-Type': 'application/x-www-form-urlencoded',
    }
    data = {
    'result': 'read_run',
    'query': 'tax_tree('+ str(animal_taxon)+')',
    'fields': 'scientific_name,tax_id,study_accession,study_alias,study_title,sample_accession,sample_alias,accession,secondary_sample_accession,run_accession,base_count,read_count,description,experiment_title,instrument_platform,instrument_model,library_source,library_selection,library_layout,library_name,library_strategy,sequencing_method,sample_description,sample_title,sex,cell_type,tissue_type,first_public,last_updated,fastq_galaxy,fastq_ftp,sra_ftp,sra_galaxy,submitted_ftp,submitted_galaxy,project_name',
    'format': 'json'
    }

    r = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', headers=headers, data=data)
    r_json =r.json()
    df = pd.DataFrame.from_dict(r_json)

    ids = df['sample_accession'].unique().tolist()
    json_list = []
    biosamples_server = 'https://www.ebi.ac.uk/biosamples/samples/'
    i=0
    for id in ids:
        i+=1
        ext = id
        r = requests.get(biosamples_server+ext, headers={ "Content-Type" : "application/json"})

        if r.ok:
            decoded = r.json()
            json_list.append(decoded)
        else:
            continue
        print(str(i) + '/' + str(len(ids)))

    #with open('/home/pelmo/data_and_pipelines/json_list', 'rb') as pickled_list:
    #    json_list =pickle.load(pickled_list)
    ####BIOSAMPLES#####

    biosamples_df = pd.DataFrame.from_records(json_list)
    biosamples_df_chars = biosamples_df['characteristics'].apply(pd.Series)
    #biosamples_df_chars.columns.tolist()
    interesting_fields = []
    interesting_fields = ['sex', 'breed', 'tissue']
    biosamples_ids = biosamples_df_chars['External Id'].apply(pd.Series)[0].apply(pd.Series)['text']
    biosamples_ids.name = 'biosamples_id'
    fields_df = pd.DataFrame(index=np.arange(len(biosamples_df)), columns=interesting_fields)
    for field in interesting_fields:
        fields_df[field] = biosamples_df_chars[field].apply(pd.Series)[0].apply(pd.Series)['text']

    df = df.drop(columns = ['sex'])
    biosamples_id_fields = pd.concat([biosamples_ids, fields_df], axis=1)

    df_with_biosamples = df.merge(biosamples_id_fields, left_on='sample_accession', right_on='biosamples_id')
    df_with_biosamples = df_with_biosamples.rename(columns={'tissue_type':'tissue_type_ENA', 'tissue':'tissue_biosamples'})

    #####XREF#####
    project_ids = df_with_biosamples['study_accession'].unique().tolist()
    headers={ "Content-Type" : "application/json"}
    json_list = []
    xrefs_df = pd.DataFrame(project_ids, columns = ['study_accession'])
    xrefs_df['xrefs'] = ''
    xrefs_list = []
    for id in project_ids:
        r = requests.get('https://www.ebi.ac.uk/ena/xref/rest/json/search?accession=' + id, headers=headers)
        if r.ok:
            decoded = r.json()
            xrefs_list.append([d['Source Primary Accession'] for d in decoded])
    xrefs_df['xrefs'] = xrefs_list
    df_complete = df_with_biosamples.merge(xrefs_df, left_on='study_accession', right_on='study_accession')

    return(df_complete)

def get_stats(df_complete):
    df = df_complete
    df = df.replace('', 'NaN')

    samples_number = len(df['sample_accession'].unique().tolist())
    samples_number
    projects_number = len(df['study_accession'].unique().tolist())
    projects_number
    samples_with_sex = len(df['sex'].dropna())
    samples_with_sex
    instrument_platform = len(df['instrument_platform'].dropna())
    instrument_platform
    paired_end_samples = len(df[df['library_layout'] == 'PAIRED'])
    paired_end_samples
    single_end_samples = len(df[df['library_layout'] == 'SINGLE'])
    single_end_samples
    library_strategy = len(df['library_strategy'].dropna())
    library_strategy
    tissue_biosamples = len(df['tissue_biosamples'].dropna())
    tissue_biosamples
    tissue_ena = len(df[df['tissue_type_ENA']!='NaN'])

    general_stats = {'samples_number': samples_number, 'projects_number':projects_number,'samples_with_sex':samples_with_sex, 'instrument_platform':instrument_platform, 'paired_end_samples':paired_end_samples, 'single_end_samples':single_end_samples, 'library_strategy':library_strategy, 'tissue_biosamples':tissue_biosamples}

    df['base_count']=df['base_count'].replace('NaN', 0)
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
    breeds_w = len(deep_df['breed'].dropna())
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

    wgs_illumina_stats = {'samples_number' : samples_number_w, 'projects_number': projects_number_w, 'breeds':breeds_w, 'median_depth':median_depth, 'mean_depth':mean_depth, 'depth_sd':depth_sd, 'min_depth': min_depth, 'max_depth':max_depth }
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
    breed_stats
    return(df)

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


def save_files(df, path):
    #create directory to save files
    directory = str(animal_taxon) + '_data'
    path = os.path.join(directory_path, directory)
    os.mkdir(path)
    df.to_csv(os.path.join(path, str(animal_taxon)+'_database.csv'))
    return()

if __name__ == "__main__":
    animal_taxon = sys.argv[1]
    directory_path = sys.argv[2]
    save_files(get_data(animal_taxon), directory_path)
