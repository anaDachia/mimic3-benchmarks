from __future__ import absolute_import
from __future__ import print_function

import pandas as pd


def dataframe_from_csv(path, header=0, index_col=0):
    return pd.read_csv(path, header=header, index_col=index_col)

#
# def create_symptom_file():
#     old_edges_df = pd.read_csv("~/Desktop/Edge.csv", header=None, names=["NODE1", "NODE2", "COUNT", "TYPE"], delimiter="\t")
#     old_edges_df = old_edges_df[old_edges_df["TYPE"] == "pati.symp"]
#     old_edges_df.to_csv("~/Desktop/symptoms.csv", index_label='NODE1', mode="w", header=False, index=False )
