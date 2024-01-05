import argparse
from Bio import SeqIO
import pandas as pd
import shutil
import os
import numpy as np

parser = argparse.ArgumentParser(description='Get top guides and produce four outut files: targets.fa windows.txt targets.txt and spacers.txt')
parser.add_argument("-i", "--input", default="3.guides/guides_table.tsv",
                    help="input tab delimited file, contains guides, windows, coverage, produced by get_guides_table.py") 
# parser.add_argument("-f", "--cov", default=0.3,
#                     help="minimum coverage of segment variants, default is 0.5")
parser.add_argument("-m", "--metric", default="strain_targeted_count", help = "name of the column used to rank guides, possible values are 'segment:targeted_variant_count','subtype_targeted_count', 'strain_targeted_count'")
parser.add_argument("-c", "--cutoff", default=1000, help = "guides with metric value less than cutoff will be removed, default is 1000")
parser.add_argument("-o", "--outputdir", default="top_guides_out",
                    help="output directory, default is top_guides_out")
                    
args = parser.parse_args()

# create output directory
# if os.path.exists(args.output):
#     shutil.rmtree(args.output)
if not os.path.exists(args.outputdir):
    os.mkdir(args.outputdir)


# load the guide frequency file
df_guides = pd.read_csv(args.input, sep='\t', low_memory=False)

# set column names
df_guides.columns = ["segment (from sequence desc.)", 'guide', 'segment:targeted_variant_count',"subtype_targeted_count", "strain_targeted_count", "sequenceID (bvbrc)","start","target","spacer","strand","GC_content","A_content"]

#remove rows where guides are invalid
df_guides.drop(df_guides[(df_guides['guide'] == "nnnnnnnnnnnnnnnnnnnn") | (df_guides['guide'] == "aaaaaaaaaaaaaaaaaaaa") | (df_guides['guide'] == "target")].index, inplace=True)

df_guides["segment (from sequence desc.)"] = df_guides["segment (from sequence desc.)"].fillna(0).astype(int).astype(str)
df_guides["segment (from sequence desc.)"].replace("0", np.nan, inplace=True)

df_guides["segment:targeted_variant_count"] = df_guides["segment:targeted_variant_count"].fillna(pd.NA)

#remove rows where segment is invalid
#df_guides.dropna(subset=["segment"], inplace=True)

#reindex
df_guides.reset_index(drop=True, inplace=True)


# check if selected metric is in the dataframe
if args.metric not in df_guides.columns:
    print(f"Error: selected metric:{args.metric} is not in the dataframe")
    exit(1)

# filter by multi-targeting count and subtypes count
mask = df_guides[args.metric]  >= float(args.cutoff)
df_guides = df_guides[mask]

print(df_guides.shape[0], " guides were extracted")

# remove duplicates guides
df_guides.drop_duplicates(subset=['guide'], keep='first', inplace=True)

print(df_guides.shape[0], " guides remained after removing duplicates")

# sort by metric
df_guides.sort_values(by=[args.metric], ascending=False, inplace=True)

# write guides dataframe
df_guides.to_csv(os.path.join(args.outputdir,"df_top_guides.tsv"), sep='\t', index=False)

# write targets.fa
with open(os.path.join(args.outputdir,"targets.fa"), "w") as f:
    for index, row in df_guides.iterrows():
        f.write(f">{row['sequenceID (bvbrc)']}_{row['start']}\n")
        f.write(f"{row['target'].upper()}\n")
# write windows.txt
with open(os.path.join(args.outputdir,"windows.txt"), "w") as f:
    f.write("segment\tstart\ttarget\tspacer\tstrand\tGC_content\tA_content\n")
    for index, row in df_guides.iterrows():
        f.write(f"{row['sequenceID (bvbrc)']}\t{row['start']}\t{row['target'].upper()}\t{row['spacer'].upper()}\t{row['strand']}\t{row['GC_content']}\t{row['A_content']}\n")
# write targets.txt
with open(os.path.join(args.outputdir,"targets.txt"), "w") as f:
    for index, row in df_guides.iterrows():
        f.write(f"{row['target'].upper()}\n")
# write spacers.txt
with open(os.path.join(args.outputdir,"spacers.txt"), "w") as f:
    for index, row in df_guides.iterrows():
        f.write(f"{row['spacer'].upper()}\n")