import argparse
import sys
import pandas as pd
import ast
import numpy as np
import os

def main():
    parser = argparse.ArgumentParser(description="This script parses the stats file generated by methphaser and"
                                                 "reports more detailed information for each reported phase block pair")
    parser.add_argument("input_stats",nargs="+")
    parser.add_argument("output_stats")
    parser.add_argument("--min_reads",default=10,type=int)
    parser.add_argument("--min_percentage",default=0.5,type=float)
    parser.add_argument("--min_bias", default=0.05,type=float)
    args = parser.parse_args()

    combined_out = pd.DataFrame()
    for f in args.input_stats:
        try:
            stats = pd.read_csv(f, index_col=0)[1:].reset_index()
            if "myth_phasing_relationship" not in stats.columns:
                continue

            out = pd.DataFrame()
            stats["snp_phased_block_1"] = stats["snp_phased_block_1"].apply(ast.literal_eval)
            stats["snp_phased_block_2"] = stats["snp_phased_block_2"].apply(ast.literal_eval)
            out["Contig"] = [os.path.basename(os.path.dirname(f))] * stats.shape[0]
            out["Blockset"] = [os.path.splitext(os.path.basename(f))[0]] * stats.shape[0]
            out["Block1"] = stats.apply(get_start, block=1, axis=1)
            out["Block2"] = stats.apply(get_start, block=2, axis=1)
            out["Distance"] = stats.apply(calc_distance, axis=1)
            out["Total"] = stats.apply(calc_total, axis=1)
            out["Same"] = stats.apply(get_count, name="same_hap_num", axis=1)
            out["Flip"] = stats.apply(get_count, name="diff_hap_num", axis=1)
            out["Diff"] = stats.apply(calc_diff, axis=1)
            out["B1_CpG1"] = stats.apply(get_count, name="olp_b1_hp1_CpG_num", axis=1)
            out["B1_CpG2"] = stats.apply(get_count, name="olp_b1_hp2_CpG_num", axis=1)
            out["B1_Bias"] = stats.apply(calc_bias, block=1, axis=1)
            out["B2_CpG1"] = stats.apply(get_count, name="olp_b2_hp1_CpG_num", axis=1)
            out["B2_CpG2"] = stats.apply(get_count, name="olp_b2_hp2_CpG_num", axis=1)
            out["B2_Bias"] = stats.apply(calc_bias, block=2, axis=1)
            out.insert(5, "Call", out.apply(classify_pair, min_reads=args.min_reads, min_diff=args.min_percentage,
                                              min_bias=args.min_bias, axis=1))
            combined_out = pd.concat([combined_out, out], axis=0)
        except pd.errors.EmptyDataError:
            pass
    combined_out.to_csv(args.output_stats, sep="\t", index=False)


def classify_pair(row, min_reads, min_diff, min_bias):
    if row['Total'] == 0:
        result = "no_overlap"
    elif row['Total'] <= min_reads:
        result = "low_overlap"
    elif row["Diff"] <= min_diff:
        result = "low_difference"
    elif row["B1_Bias"] < min_bias or row["B2_Bias"] < min_bias:
        result = "biased_haplotyping"
    elif row["Same"] > row["Flip"]:
        result = "same"
    else:
        result = "flip"
    return result


def get_count(row, name):
    if pd.isna(row[name]):
        return 0
    return int(row[name])


def get_start(row, block):
    return row["snp_phased_block_{0}".format(block)][0]


def calc_distance(row):
    return row['snp_phased_block_2'][0] - row['snp_phased_block_1'][1]


def calc_total(row):
    if pd.isna(row["same_hap_num"]):
        return 0
    return int(row['same_hap_num'] + row['diff_hap_num'])


def calc_diff(row):
    if pd.isna(row["same_hap_num"]):
        return "NA"
    same_p = row['same_hap_num'] / (row['same_hap_num'] + row['diff_hap_num'])
    flip_p = row['diff_hap_num'] / (row['same_hap_num'] + row['diff_hap_num'])
    return round(max([flip_p, same_p]) - min([flip_p, same_p]), 3)


def calc_bias(row, block):
    if pd.isna(row["same_hap_num"]):
        return "NA"
    h1 = row["olp_b{0}_hp1_CpG_num".format(block)]
    h2 = row["olp_b{0}_hp2_CpG_num".format(block)]
    return round(min([h1, h2]) / sum([h1, h2]), 3)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("User interrupted, exiting")
        sys.exit(1)
