import csv
import pandas as pd
import requests
import sys 
import os
from io import StringIO

def enrich_df(gvm_df, curated_df):
    gvm_df = pd.read_csv(gvm_df)
    curated_df = pd.read_csv(curated_df)
    complete_df = gvm_df.merge(curated_df, left_on='secondary_sample_accession', right_on='secondary_sample_accession', how = 'outer')
    complete_df.to_csv("complete_db.csv")
    print(complete_df.info())
    return()


if __name__ == "__main__":
    df_to_enrich = sys.argv[1]
    enrichment_df = sys.argv[2]
    enrich_df(df_to_enrich, enrichment_df)