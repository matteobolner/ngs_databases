import csv
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import sys
import os
import json
import pickle


def get_data(animal_taxon, directory_path):
    animal_taxon = 9822

    #create directory to save files
    directory = str(animal_taxon) + '_data'
    path = os.path.join(directory_path, directory)
    os.mkdir(path)

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
        print(i)


    '''with open('/home/pelmo/data_and_pipelines/json_list', 'rb') as pickled_list:
        json_list =pickle.load(pickled_list)'''

    biosamples_df = pd.DataFrame.from_records(json_list)
    biosamples_df_chars = biosamples_df['characteristics'].apply(pd.Series)
    biosamples_ids = biosamples_df_chars['External Id'].apply(pd.Series)[0].apply(pd.Series)['text']
    biosamples_breed = biosamples_df_chars['breed'].apply(pd.Series)[0].apply(pd.Series)['text']
    biosamples_id_breed = pd.concat([biosamples_ids, biosamples_breed], axis=1)

    biosamples_id_breed.columns = ['biosamples_id', 'breed']

    df = biosamples_id_breed.merge(df, left_on='biosamples_id', right_on='sample_accession')

    #return some information about the libraries
    print(pd.value_counts(df['library_source']))
    print(pd.value_counts(df['library_strategy']))
    print(pd.value_counts(df['breed']))



    multiple_runs = df['run_accession'].value_counts()
    multiple_runs
    multiple_runs = multiple_runs[multiple_runs.gt(1)]
    multiple_runs
    multiple_samples = df['sample_accession'].value_counts()
    multiple_samples
    multiple_samples = df[df['sample_accession'].isin(multiple_samples.index[multiple_samples.gt(1)])][['study_accession','sample_accession']]
    multiple_samples['count'] = multiple_samples.groupby(['sample_accession'])['study_accession'].transform('count')
    multiple_samples = multiple_samples.drop_duplicates().sort_values(by= 'count', ascending=False).reset_index(drop=True)

    multiple_samples
    subspecies = df['scientific_name'].value_counts()
    print(subspecies)

    #sample_plot = sns.barplot(data = multiple_samples, x = 'sample_accession', y = 'count', palette = 'muted')
    #plt.xticks(rotation=90)

    #plt.show()

    grouped_df = df.groupby(df.library_source)
    library_sources = [name for name,unused_df in grouped_df]
    for source in library_sources:
        grouped_df.get_group(source).to_csv(os.path.join(path, str(animal_taxon) + "_" + str(source) + "_database.csv"))
        pd.value_counts(grouped_df.get_group(source)['library_strategy']).to_csv(os.path.join(path,  str(animal_taxon) + "_" + str(source) + "_library_strategies.csv"))

    df.to_csv(os.path.join(path, str(animal_taxon)+'_database.csv'))
    multiple_runs.to_csv(os.path.join(path, str(animal_taxon) + '_multiple_runs.csv'))
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
    return(df)

if __name__ == "__main__":
    animal_taxon = sys.argv[1]
    directory_path = sys.argv[2]
    get_data(animal_taxon, directory_path)
