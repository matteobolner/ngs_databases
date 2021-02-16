import sys
import pandas as pd

def split_df(input_df):
    df = pd.read_csv(input_df, index_col=0)
    print(df.info())
    grouped_df = df.groupby(df.library_source)
    #for key, item in grouped_df:
    #    print(grouped_df.get_group(key))
    #print(grouped_df['GENOMIC'])
    #print(grouped_df.keys)
    library_sources = [name for name,unused_df in grouped_df]
    for source in library_sources:
        grouped_df.get_group(source).to_csv(str(source)+"_database.csv")
        pd.value_counts(grouped_df.get_group(source)['library_strategy']).to_csv(str(source)+"_library_strategies.csv")
    return()



if __name__ == "__main__":
    input_df = sys.argv[1]
    split_df(input_df)