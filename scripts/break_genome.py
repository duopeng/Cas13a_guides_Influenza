import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default="ref_data/NC_045512v2.fa",
                    help="genome .fa file name")
parser.add_argument("-d", "--scriptdir", default=".",
                    help="dummy argument for compatibility with break_genomes.py")
parser.add_argument("-w", "--window", type=int, default=20,
                    help="window size")
parser.add_argument("-e", "--enzyme", type=str, help="Cas enzyme type")  # "Cas13a" or "Cas12"
parser.add_argument("-p", "--pfs_length", type=int, default=4,
                    help="length of PFS/PAM to evaluate")
parser.add_argument("-g", "--genome_strand", type=str, default="+",
                    help="strandedness of viral genome")
parser.add_argument("-s", "--strand", type=str, default="+",
                    help="strand to target")
parser.add_argument("-o", "--out", type=str, default="break_genome_output",
                    help="output directory")
args = parser.parse_args()

if args.enzyme is None:
    print("\nNO ENZYME SPECIFIED\n")
    exit(1)
if not os.path.exists(args.out):
    os.makedirs(args.out)


genome_seq = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))

windows = []
for seq_id, seq in genome_seq.items():
    seq_str = str(seq.seq)
    num_windows = len(seq_str) - args.window - 3
    for i in range(num_windows):
        windows.append({
            'segment': seq_id.split(" ")[0],
            'start': i + 1,
            'target': seq_str[i: i + args.window]
        })

df = pd.DataFrame(windows)
num_rows, num_columns = df.shape
if num_rows == 0:
    print(f"No windows found in genome {args.input}. Exiting...")
    exit(1)
if args.enzyme == "Cas13a":
    if args.strand == args.genome_strand:
        df['spacer'] = df['target'].apply(lambda x: Seq(x).reverse_complement().transcribe())
        df['strand'] = args.genome_strand
    else:
        df['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
        df['target'] = df['target'].apply(lambda x: Seq(x).reverse_complement())
        df['strand'] = args.strand
else:
    df_minusStrand = df.copy()
    df_minusStrand['target'] = df['spacer'].apply(lambda x: Seq(x).transcribe())
    df_minusStrand['spacer'] = df['target'].apply(lambda x: Seq(x).transcribe())
    df_minusStrand['strand'] = '-'
    df = pd.concat([df, df_minusStrand], ignore_index=True)

df['GC_content'] = df['spacer'].apply(lambda x: round(gc_fraction(Seq(x)),2))
df['A_content'] = df['spacer'].apply(lambda x: x.upper().count('A')/args.window)

df.to_csv(f"{args.out}/windows.txt", sep='\t', index=False)

df['target'].to_csv(f"{args.out}/targets.txt", index=False)
df['spacer'].to_csv(f"{args.out}/spacers.txt", index=False)

with open(f"{args.out}/targets.fa", 'w') as f:
    for idx, row in df.iterrows():
        f.write(f">{row['segment']}_{row['start']}\n")
        f.write(f"{row['target']}\n")
