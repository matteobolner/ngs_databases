import csv
from io import StringIO
import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
import sys
import os

def get_data(animal_taxon, directory_path):
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
    'format': 'tsv'
  }

  response = requests.post('https://www.ebi.ac.uk/ena/portal/api/search', headers=headers, data=data)

  data = StringIO(response.text)
  df = pd.read_csv(data, sep = '\t')

  #return some information about the libraries
  print(pd.value_counts(df['library_source']))
  print(pd.value_counts(df['library_strategy']))




  multiple_runs = df['run_accession'].value_counts()
  multiple_runs = multiple_runs[multiple_runs.gt(1)]

  multiple_samples = df['sample_accession'].value_counts()
  multiple_samples = df[df['sample_accession'].isin(multiple_samples.index[multiple_samples.gt(1)])][['study_accession','sample_accession']]
  multiple_samples['count'] = multiple_samples.groupby(['sample_accession'])['study_accession'].transform('count')
  multiple_samples = multiple_samples.drop_duplicates().sort_values(by= 'count', ascending=False).reset_index(drop=True)

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
