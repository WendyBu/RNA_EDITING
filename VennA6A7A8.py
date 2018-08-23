"""
Compare A6, A7, A8
EZH2 knockdown
"""


import pandas as pd
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


def get_list(input_f):
    df = pd.read_csv(input_f, sep="\t")
    chrom_pos = df["chr_pos"]
    chrom_pos.drop_duplicates(inplace=True)
    l_position = chrom_pos.tolist()
    return l_position


def get_gene_list(input_f):
    df = pd.read_csv(input_f, sep="\t")
    dfs = df[df["gene"]!="intergenic"]
    df_gene = dfs["gene"]
    df_gene.drop_duplicates(inplace=True)
    return df_gene



def draw_venn3(lst1,lst2,lst3):
    venn3([set(lst1), set(lst2), set(lst3)], set_labels=('Control', 'EZH2_sh1', 'EZH2_sh2'))
    plt.title('Comparison RNA edited genes in EZH2 knockdown\n')
    plt.show()
    pass


def call_venn():
    # RNA editing sites
    lst1=get_list("EZH2Rediting/A6all1.xls")
    lst2=get_list("EZH2Rediting/A7all1.xls")
    lst3=get_list("EZH2Rediting/A8all1.xls")
    # draw_venn3(lst1,lst2,lst3)
    print "editing sites:", len(lst1), len(lst2), len(lst3)

    # RNA editing genes
    genelst1=get_gene_list("EZH2Rediting/A6all1.xls")
    genelst2=get_gene_list("EZH2Rediting/A7all1.xls")
    genelst3=get_gene_list("EZH2Rediting/A8all1.xls")
    draw_venn3(genelst1,genelst2,genelst3)
    print "editing genes:", len(genelst1), len(genelst2), len(genelst3)
    return lst1,lst2,lst3,genelst1,genelst2,genelst3



def diff_genes():





    pass

def main():
    lst1,lst2,lst3, glst1,glst2,glst3= call_venn()

    pass


if __name__ =="__main__":
    main()