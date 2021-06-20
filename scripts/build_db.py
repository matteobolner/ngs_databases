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
    genome_sizes = {'9822':'2501912388', '9984':'2737490501', '9913':'2715853792', '7460':'225,250,884'}
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
