import pandas as pd
import numpy as np
from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt
pd.set_option('display.max_columns', 50)


def draw_venn2(list1, list2):
    venn2([set(list1), set(list2)], set_labels=('EZH2', 'ADAR'))
    plt.title('EZH2 & ADAR regulated genes\n')
    plt.show()
    pass


def get_getlist(df, cutoff1, cutoff2):
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]>cutoff1) | (df["log2(fold_change)"]<cutoff2)]
    df_list = df_c["gene_id"].tolist()
    return df_list


def up(df, cutoff1):
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]>cutoff1)]
    df_list = df_c["gene_id"].tolist()
    return df_list

def down(df, cutoff2):
    df = df[~df["log2(fold_change)"].isin(["#NAME?"])]
    df["log2(fold_change)"] = df["log2(fold_change)"].apply(lambda x: np.float(x))
    df_c = df [(df["log2(fold_change)"]<cutoff2)]
    df_list = df_c["gene_id"].tolist()
    return df_list



Ezh = pd.read_csv("EZH2_gene_exp.diff", sep="\t")
ADAR = pd.read_csv("ADAR_gene_exp.diff", sep="\t")

#coregulated
# Ezh_list = get_getlist(Ezh, 1, -1)
# ADAR_list = get_getlist(ADAR, 1, -1)

# co_upregulated
# Ezh_list = up(Ezh, 1)
# ADAR_list = up(ADAR, 1)

# co_downregulated
Ezh_list = down(Ezh, -1)
ADAR_list = down(ADAR, -1)


draw_venn2(Ezh_list, ADAR_list)

