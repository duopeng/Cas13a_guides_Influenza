import argparse
from Bio import SeqIO

def search_keyword_in_fasta(fasta_file, keyword):
    with open(fasta_file + f".{keyword}.fa", "w") as handle, open(fasta_file + f".{keyword}.txt", "w") as handle2:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if keyword in record.description:
                # print(f"Keyword found in record {record.id}")
                # print(f"Description: {record.description}")
                # print(f"Sequence: {record.seq}\n")
                handle.write(f">{record.description}\n{record.seq}\n")
                handle2.write(f"{record.id}\t{len(record.seq)}\t{record.description}\n")

def main():
    parser = argparse.ArgumentParser(description='Search for keyword in FASTA file.')
    parser.add_argument('-f', '--file', required=True, help='Path to the FASTA file.')
    parser.add_argument('-k', '--keyword', required=True, help='Keyword to search for.')
    args = parser.parse_args()

    search_keyword_in_fasta(args.file, args.keyword)

if __name__ == "__main__":
    main()