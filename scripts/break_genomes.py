import argparse
from Bio import SeqIO
import shutil
import os
import concurrent.futures
import subprocess

parser = argparse.ArgumentParser(description='Wrapper for break_genome.py. Break each genome in the input file with break_genome.py.')
parser.add_argument("-d", "--scriptdir", default=".",
                    help="directory containing break_genome.py")
parser.add_argument("-i", "--input", default="SARS_CoV2_genomes/Coronaviridae.fna.SARS-CoV-2.fa",
                    help="genomes .fa file name")
parser.add_argument("-w", "--window", type=int, default=20,
                    help="window size")
parser.add_argument("-e", "--enzyme", type=str, help="Cas enzyme type")  # "Cas13a" or "Cas12"
parser.add_argument("-p", "--pfs_length", type=int, default=4,
                    help="length of PFS/PAM to evaluate")
parser.add_argument("-g", "--genome_strand", type=str, default="+",
                    help="strandedness of viral genome")
parser.add_argument("-s", "--strand", type=str, default="+",
                    help="strand to target")
parser.add_argument("-o", "--out", type=str, default="output",
                    help="output directory")
parser.add_argument('-n', '--num_workers', type=int, default=8, 
                    help='Number of concurrent tasks.')
args = parser.parse_args()

def run_script(args_string):
    args_list = args_string.split(" ")
    subprocess.run(['python', f'{args.scriptdir}/break_genome.py'] + args_list)


def main():
    if args.enzyme is None:
        print("\nNO ENZYME SPECIFIED\n")
        exit(1)
    # make output directory
    if not os.path.exists(args.out):
        os.makedirs(args.out)
    else:
        print(f"Purging output directory {args.out}")
        shutil.rmtree(args.out)
        os.makedirs(args.out)

    #read in genomes and write to individual fasta files
    print("Preparing genomes for break_genome.py...", flush=True)
    list_of_genomes = []
    genome_count = 0
    for record in SeqIO.parse(args.input, "fasta"):
        if not os.path.exists(f"{args.out}/{record.id}"):
            os.makedirs(f"{args.out}/{record.id}")
        with open(f"{args.out}/{record.id}/{record.id}.fa", "w") as handle:
            handle.write(f">{record.description}\n{record.seq}\n")
        list_of_genomes.append(record.id)
        genome_count += 1
        print(f"{genome_count} genomes prepared", flush=True) if genome_count % 10000 == 0 else None
    print(f"Done. {genome_count} genomes prepared")


    # prepare arguments for break_genome.py
    args_dict = vars(args)
    exclude_args = ["input", "out", "num_workers"]
    shared_args = [f"--{k} {v}" for k, v in args_dict.items() if k not in exclude_args]
    shared_args = " ".join(shared_args)

    list_of_args_string = [f"{shared_args} --input {args.out}/{genome}/{genome}.fa --out {args.out}/{genome}/" for genome in list_of_genomes]

    #print(shared_args)
    #print(args_list)

    #run break_genome.py on each genome
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        executor.map(run_script, list_of_args_string)

if __name__ == '__main__':
    main()

