### combine editing sites and mRNA expression
### in cells with EZH2 knockdown

import pandas as pd
pd.set_option('display.max_columns', 50)

df = pd.read_table("gene_exp.diff", sep="\t")
# print df.head(), df.shape
dfs = df[["gene_id", "log2(fold_change)", "p_value"]]
dfs.set_index("gene_id", inplace=True)   ### set_index, not reset_index


lost = pd.read_table("../WTvsEZH2sh/CtrvsA7A8_final.xls",sep="\t", index_col=0)
lost.sort_values("0", axis=0, inplace=True, ascending=False)
df_exp_loss = pd.concat([lost,dfs], axis=1, join='inner')
df_exp_loss.to_csv("lost_sites_exp.xls", sep="\t")


gain = pd.read_table("../WTvsEZH2sh/A8vsCtrl_final.xls",sep="\t", index_col=0)
gain.sort_values("0", axis=0, inplace=True, ascending=False)
df_exp_gain = pd.concat([gain,dfs], axis=1, join='inner')
df_exp_gain.to_csv("gained_exp.xls", sep="\t")


