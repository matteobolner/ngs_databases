#also not need anymore but may still be useful

import pandas as pd
import sys

def merge_dataframes(df_1, df_2):
    df_1 = pd.read_csv(df_1)
    df_2 = pd.read_csv(df_2)
    df_1 = df_1.drop(df_1.columns[[0]],axis = 1)
    df_2 = df_2.drop(df_2.columns[[0,1]],axis = 1)
    complete_df = pd.concat([df_1, df_2]).reset_index(drop=True)
    complete_df.to_csv("complete_db.csv")
    return(complete_df)


if __name__ == "__main__":
    df_1 = sys.argv[1]
    df_2 = sys.argv[2]
    merge_dataframes(df_1, df_2)