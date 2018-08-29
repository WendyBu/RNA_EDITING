"""
Compare A6, A7, A8
EZH2 knockdown
output : "WTvsEZH2sh/A7vsCtrl_final.xls"
Venn diagram tutorial https://www.badgrammargoodsyntax.com/compbio/2017/10/29/compbio-012-making-venn-diagrams-the-right-way-using-python

"""


import pandas as pd
import numpy as np
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
from matplotlib import pyplot as plt

pd.set_option('display.max_columns', 50)


def get_list(df):
    chrom_pos = df["chr_pos"]
    chrom_pos.drop_duplicates(inplace=True)
    l_position = chrom_pos.tolist()
    return l_position


def get_gene_df(df,genelist):
    dfs = df[df["chr_pos"].isin(genelist)]
    dfs = dfs[["chr_pos", "gene", "annot1", "alu?",	"non_alu_repetitive?"]]
    dfs.drop_duplicates(inplace=True)
    return dfs


def draw_venn3(lst1,lst2,lst3):
    venn3([set(lst1), set(lst2), set(lst3)], set_labels=('Control', 'EZH2_sh1', 'EZH2_sh2'))
    plt.title('Comparison RNA edited genes in EZH2 knockdown\n')
    plt.show()
    pass

def draw_venn2(list1,list2):
    venn2([set(list1), set(list2)], set_labels=('Control', 'EZH2 shRNA'))
    plt.title('RNA editing genes in control and EZH2 knockdown\n')
    plt.show()
    pass


def get_gene_list(df):
    gene_list = list(set(df["gene"].tolist()))
    return gene_list


def genelist_count(file, exclude_list):
    A = pd.read_csv(file, sep="\t")
    A = A[~A["gene"].isin(exclude_list)]
    A_df = A.groupby("gene").size()
    A_df.sort_values(inplace=True, ascending=False)
    A_df.drop("intergenic", axis=0, inplace=True)
    return A_df


def add_annot(df, comb_table):
    df2 = df.groupby(["gene", "annot1"]).size()
    UTR3 = df2[:, "3UTR"]
    UTR5 = df2[:, "5UTR"]
    INTRON = df2[:, "intronic"]
    df3 = pd.concat([comb_table,UTR3, UTR5, INTRON], axis=1, join="outer", sort=True)
    df3.rename({0: "total_edits", 1: "3UTR", 2:"5UTR", 3:"intronic"}, inplace=True)
    df3.replace(" ", np.nan, inplace=True)
    df3.fillna(0, inplace=True)
    return df3


def main():
    A6_df=pd.read_csv("EZH2Rediting/A6all1.xls", sep="\t")
    A7_df=pd.read_csv("EZH2Rediting/A7all1.xls", sep="\t")
    A8_df=pd.read_csv("EZH2Rediting/A8all1.xls", sep="\t")
    lst1=get_list(A6_df)
    lst2=get_list(A7_df)
    lst3=get_list(A8_df)
    draw_venn3(lst1,lst2,lst3)

    # for three groups
    gained = list((set(lst3).union(lst2)).difference(lst1))
    lost =  list(set(lst1).difference(set(lst3).union(lst2)))


    A7uniq_df = get_gene_df(A7_df, gained)
    A8uniq_df = get_gene_df(A8_df, gained)
    # A7uniq_df.to_csv("WTvsEZH2sh/A7gained.xls", sep="\t")
    # A8uniq_df.to_csv("WTvsEZH2sh/A8gained.xls", sep="\t")
    A6uniq_df = get_gene_df(A6_df,lost)
    # A6uniq_df.to_csv("WTvsEZH2sh/A6gained.xls", sep="\t")

    A7_gained_list = get_gene_list(A7uniq_df)
    A8_gained_list = get_gene_list(A8uniq_df)
    A6_lost_list = get_gene_list(A6uniq_df)
    draw_venn3(A6_lost_list,A7_gained_list,A8_gained_list)


    A7_table = A7uniq_df.groupby("gene").size()
    A8_table = A8uniq_df.groupby("gene").size()
    A6_table = A6uniq_df.groupby("gene").size()


    A6_distrib = A6_df.groupby("annot1").size()
    A7_distrib = A7_df.groupby("annot1").size()
    A8_distrib = A8_df.groupby("annot1").size()

    distrib = pd.DataFrame()
    distrib["control"] = A6_distrib
    distrib["EZH2_sh1"]= A7_distrib
    distrib["Ezh2_sh2"]= A8_distrib
    # distrib.to_csv("WTvsEZH2sh/distribution.xls", sep="\t")

    A7_final = add_annot(A7uniq_df, A7_table)
    A8_final = add_annot(A8uniq_df, A8_table)
    A6_final = add_annot(A6uniq_df, A6_table)

    A7_final.to_csv("WTvsEZH2sh/A7vsCtrl_final.xls", sep="\t")
    A8_final.to_csv("WTvsEZH2sh/A8vsCtrl_final.xls", sep="\t")
    A6_final.to_csv("WTvsEZH2sh/CtrvsA7A8_final.xls", sep="\t")   # add 3utr, 5utr, intron info for each gene
    pass


if __name__ =="__main__":
    main()