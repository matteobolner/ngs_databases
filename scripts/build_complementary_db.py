
#this script is outdated since the problem was that the original db consisted only of the specific species name and not e.g Sus

import pandas as pd
import sys
import requests
from io import StringIO


def get_complementary_samples(input_csv):
    reads_df = pd.read_csv(input_csv)
    study_accession_list = reads_df['study_accession'].unique()
    taxon_id_list = reads_df['tax_id'].unique()
    print(study_accession_list)
    complementary_df = pd.DataFrame(columns=list(reads_df.columns))
    print(complementary_df)
    i= 0
    for id in study_accession_list: 
        i +=1
        print(str(i) +"/" + str(len(study_accession_list)))
        headers = {
            'Content-Type': 'application/x-www-form-urlencoded',
        }

        data = {
        'result': 'read_run',
        'query': 'study_accession="'+ str(id)+'"',
        'fields': 'scientific_name,tax_id,study_accession,study_alias,study_title,sample_accession,sample_alias,accession,secondary_sample_accession,run_accession,base_count,read_count,description,experiment_title,instrument_platform,instrument_model,library_source,library_selection,library_layout,library_name,library_strategy,sequencing_method,sample_description,sample_title,sex,cell_type,tissue_type,first_public,last_updated,fastq_galaxy,fastq_ftp,sra_ftp,sra_galaxy,submitted_ftp,submitted_galaxy,project_name',
        'format': 'tsv'
        }

        response = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', headers=headers, data=data)

        complementary_data = StringIO(response.text)
        temp_df = pd.read_csv(complementary_data, sep = '\t')
        temp_df = temp_df.loc[temp_df['tax_id'].isin(taxon_id_list)]
        complementary_df = pd.concat([complementary_df, temp_df]).reset_index(drop=True)
    print(complementary_df.info())
    complementary_df.to_csv("complete_prova_9986.csv")
    return(reads_df)

if __name__ == "__main__":
    input_csv = sys.argv[1]
    get_complementary_samples(input_csv)
