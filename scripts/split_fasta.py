from Bio import SeqIO
import math
import argparse

def main():
    parser = argparse.ArgumentParser(description="Break a fasta file into chunks.")
    parser.add_argument('-f', '--file', required=True, help="Input fasta file.")
    parser.add_argument('-n', '--num_chunk', required=True, help="number of chunks to split into.")
    args = parser.parse_args()
    num_chunk = int(args.num_chunk)

    records = list(SeqIO.parse(args.file, "fasta"))
    num_records = len(records)
    chunk_size = math.ceil(num_records / num_chunk)  # determine chunk size

    for i in range(num_chunk):
        print(f"Writing chunk {i+1} of {num_chunk}, {chunk_size} records. Index {i * chunk_size} to {(i + 1) * chunk_size}")
        chunk_records = records[i * chunk_size:(i + 1) * chunk_size]
        with open(f"{args.file}_chunk_{i+1}.fasta", "w") as f:
            SeqIO.write(chunk_records, f, "fasta")

if __name__ == "__main__":
    main()
