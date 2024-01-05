import os
import argparse
import concurrent.futures
import pandas as pd
from collections import Counter
from shlex import quote
import time


parser = argparse.ArgumentParser(description='Combine filter results')
parser.add_argument('--df_guides', type=str, required=True, default = "")
parser.add_argument('--kraken2_res_folder', type=str, required=True,
                    help='The path to the folder containing kragen2 results (output of run_kraken2.py)',
                    default = "")
parser.add_argument('-n', '--num_workers', type=int, default=32, help="number of workers for parallel processing of kraken2 result parsing")
parser.add_argument('--RNAfold_res_path', type=str, required=True, default = "")
parser.add_argument('--bowtie_hg_res_path', type=str, required=True, default = "")
parser.add_argument('--bowtie_infbcd_res_path', type=str, required=False, default = "")

args = parser.parse_args()


def chunk_list(lst, chunk_size):
    """
    Breaks a list into chunks of the specified size.
    The last chunk may have fewer items if the list size is not divisible by the chunk size.
    """
    return [lst[i:i+chunk_size] for i in range(0, len(lst), chunk_size)]


def worker_readcsv(filepath, white_list = [0, 131567, 1, 488241, 1173138,335341]):
    """Read Kraken2 results from a file and return the guide name and a dataframe
       Filters based on a white list of taxids
    Args:
        filepath (_type_): Kraken2 result file
    Returns:
        guide_name (str): name of the guide
        dataframe: taxid (not in while list), hammingDistance
    """
    df = pd.read_csv(filepath, sep='\t', header=None)
    df.columns = ["classification", "sequence_ID", "taxonomy_ID", "sequence_length", "LCA"]
    df["taxid"] = df["LCA"].apply(lambda x: x.split(":")[0])
    df["taxid"] = df["taxid"].astype(int)
    df["hd"] = df["sequence_ID"].apply(lambda x: x.split("_")[1])
    # subset df to only include taxids that not are in the white list
    df = df[~df["taxid"].isin(white_list)]
    #get guide name
    filename = os.path.basename(filepath)
    guide_name = filename.split("_hd.4.txt.kraken2.out")[0]
    return [guide_name,df[["taxid","hd"]]]

def read_kraken2_res(kraken2_res_folder):
    starttime = time.time()
    filelist = os.listdir(kraken2_res_folder)
    filepathlist = [os.path.join(kraken2_res_folder, filename) for filename in filelist]
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        results = executor.map(worker_readcsv, filepathlist)
    endtime = time.time()
    print(f"read_kraken2_res took {endtime - starttime:.2f} seconds, read {len(filepathlist)} files)")
    return results

def main():

    ############################
    # read the guides dataframe#
    ############################
    if os.path.exists(args.df_guides):
        print(f"Reading df_guides dataframe")
        df_guides = pd.read_csv(args.df_guides, sep='\t')
        print(f"df_guides shape: {df_guides.shape}")
    else:
        exit(f"df_guides file does not exist\n path={args.df_guides}")

    #################################
    # read and parse kraken2 results#
    #################################
    if os.path.exists(args.kraken2_res_folder):
        print(f"Reading kraken2 results using {args.num_workers} workers")
        kraken2_res = read_kraken2_res(args.kraken2_res_folder) # takes 30 seconds with 64 cores
        kraken2_res = {guide_name:df for guide_name, df in kraken2_res} # convert to dictionary
        # parse through each guide and determine cross-reactivity with common microbes, takes ~15 seconds
        print(f"Calculating cross-reactivity")
        start = time.time()
        guide_crossreact_df = pd.DataFrame(columns=["hd1", "hd2", "hd3", "hd4"])
        for guide_name, df in kraken2_res.items(): # go through each guide
            match_counts = {"hd1":0, "hd2":0, "hd3":0, "hd4":0}
            for hd in df["hd"].unique(): # go each hamming distance
                df_subset = df[df["hd"] == hd]
                taxid_counts = sum(Counter(df_subset["taxid"]))
                match_counts[hd] = taxid_counts
            guide_crossreact_df.loc[guide_name] = pd.Series({"hd1": match_counts["hd1"], "hd2": match_counts["hd2"], "hd3":  match_counts["hd3"], "hd4": match_counts["hd4"]})
        end = time.time()
        print(f"Calculating cross-reactivity took {end - start:.2f} seconds")
    else:
        exit(f"kraken2_res_folder does not exist\n path={args.kraken2_res_folder}")

    #################################################################
    # merge the cross-reactivity dataframe with the guides dataframe#
    #################################################################
    df_guides["guideID"] = df_guides["sequenceID (bvbrc)"] + "_" + df_guides["start"].apply(str)
    df_guides = df_guides.merge(guide_crossreact_df, left_on="guideID", right_index=True)

    # read RNAfold results
    if os.path.exists(args.RNAfold_res_path):
        print(f"Reading RNAfold results")
        df_RNAfold = pd.read_csv(args.RNAfold_res_path, sep='\t')
        print(f"df_RNAfold shape: {df_RNAfold.shape}")
        df_RNAfold["guideID"] = df_RNAfold["segment"] + "_" + df_RNAfold["start"].apply(str)
        df_RNAfold = df_RNAfold[["guideID", "MFE", "has_hairpin", "spacer_basepairs"]]
        df_guides = df_guides.merge(df_RNAfold, on="guideID")
    else:
        exit(f"RNAfold_res_path does not exist\n path={args.RNAfold_res_path}")

    ############################################
    # read bowtie cross-reference check results#
    ############################################
    if os.path.exists(args.bowtie_hg_res_path):
        print(f"Reading bowtie cross-reference check results")
        df_bowtie = pd.read_csv(args.bowtie_hg_res_path, sep='\t')
        print(f"df_bowtie shape: {df_bowtie.shape}")
        df_bowtie["guideID"] = df_bowtie["segment"] + "_" + df_bowtie["start"].apply(str)
        df_bowtie["bt_hg"] = df_bowtie["ct"]
        df_bowtie = df_bowtie[["guideID", "bt_hg"]]
        df_guides = df_guides.merge(df_bowtie, on="guideID")
    else:
        exit(f"bowtie_hg_res_path does not exist\n path={args.bowtie_hg_res_path}")

    if os.path.exists(args.bowtie_hg_res_path) and os.path.exists(args.bowtie_infbcd_res_path) and args.bowtie_infbcd_res_path != "":
        print(f"Reading bowtie cross-reference check results")
        df_bowtie = pd.read_csv(args.bowtie_infbcd_res_path, sep='\t')
        print(f"df_bowtie shape: {df_bowtie.shape}")
        df_bowtie["guideID"] = df_bowtie["segment"] + "_" + df_bowtie["start"].apply(str)
        df_bowtie["bt_bcd"] = df_bowtie["ct"]
        df_bowtie = df_bowtie[["guideID", "bt_bcd"]]
        df_guides = df_guides.merge(df_bowtie, on="guideID")

    ############################################
    # write the combined dataframe to a file   #
    ############################################
    print(f"Writing combined dataframe to a file")
    outfile = os.path.splitext(args.df_guides)[0] + "_filters_res_combined.tsv"
    df_guides.to_csv(outfile, sep='\t', index=False)
    print("Done...wrote result to ", outfile)

if __name__ == "__main__":
    main()