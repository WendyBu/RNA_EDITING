import pandas as pd
import glob
pd.set_option('display.max_columns', 100)
import sys
### run in the server, it is slow


def get_df2(file):
    df = pd.read_csv(file, sep="\t", header=None, low_memory=False)
    df2 = pd.DataFrame()
    df2[["chr", "start", "end", "paired"]] = df.iloc[:, 0:4]
    return df2


def editing_input(input_file):
    A8 = pd.read_csv("WTvsEZH2sh/A8gained.xls", sep="\t", index_col=0)
    A8_df = pd.DataFrame()
    A8_df[["chr-pos", "gene", "annot1"]] = A8.iloc[:, 0:3]
    # print A8["annot1"].unique()
    A8_df_utr = A8_df[A8_df["annot1"]=="3UTR"]
    A8_df_utr[["chr", "pos"]] = A8_df_utr["chr-pos"].str.split("_", expand=True)
    A8_df_utr["pos"] = pd.to_numeric(A8_df_utr["pos"])
    return input_df


def findmatch(position, df2):
    df_match = df2[(df2["start"]<position) & (df2["end"]>position)]
    if not df_match.empty:
        return df_match



def main(arg1, arg2):
    result = pd.DataFrame(columns=["chr", "start", "end", "paired"])
    RNA_editing_input = editing_input(arg1)
    for index, row in RNA_editing_input.iterrows():
        for file in glob.glob("miR_all_location/*"):
            df2 = get_df2(file)
            match = findmatch(row["pos"], df2)
            if match is not None:
                result = pd.concat([result, match], axis=0)
    result.to_csv(arg2, sep="\t")
    pass


if __name__ =="__main__":
    main(sys.argv[1], sys.argv[2])
