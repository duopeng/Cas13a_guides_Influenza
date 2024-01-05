#!/usr/bin/env python3

import os
import argparse
import glob
from collections import defaultdict
from math import ceil

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import screed
from utils import InputStream, OutputStream, cat_files, fasta_iter, command_output, multiprocessing_map, command


M_SCHEMA = {
    "target_name": str,
    "target_seq": str,
    "kmer_seq": str,
    "hamming_dist": int,
    "hits_counts": int,
    "rc": int
}


def reverse_complement(dna):
    """ Find the reverse complement of a DNA string """
    dnadict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    reverseDna = [dnadict[c] for c in dna]
    return reverseDna[::-1]


def hamming_distance(str1, str2):
    """ Compute the Hamming distance between two strings """
    hd = 0
    if len(str1) != len(str2):
        print("the two strings different length ERROR")
        return None
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hd += 1
    return hd


def read_candidate_spacers(candidate_filename):
    candidates = defaultdict()
    with open(candidate_filename) as stream:
        # skip header
        next(stream)
        for line in stream:
            line = line.strip("\n").split("\t")
            start_pos = str(line[0])
            target_kmer = line[1]
            strand = 1 if line[3] == "+" else 0
            target_kname = f"sp.{start_pos}_strand.{strand}_hd.4"
            candidates[target_kname] = target_kmer
    return candidates


def read_matched_neighbors(matched_file):
    global candidates
    # extract target kmer name and look up target kmer sequence
    matched_filename = os.path.basename(matched_file)
    kmer_name = matched_filename.split('.txt')[0]
    target_kmer = candidates[kmer_name]

    # read in the matched kmer counts
    matched_neighbors = list()
    with open(matched_file) as stream:
        for line in stream:
            line = line.strip('\n').split(' ')
            hamming_dist = hamming_distance(line[0], target_kmer)
            rc = 0
            if hamming_dist > 4:
                hamming_dist = hamming_distance(reverse_complement(line[0]), target_kmer)
                rc = 1
            matched_neighbors.append((target_kmer, line[0], hamming_dist, line[1], rc))
    return matched_neighbors


def parse_matched_kmer(args):
    matched_kmer_files = glob.glob(f"{args.jellyfish_ctdir}/*.txt")
    matched_kmernames = [os.path.basename(x).split(".txt")[0] for x in matched_kmer_files]
    list_of_matched_recrods = multiprocessing_map(read_matched_neighbors, matched_kmer_files, args.num_cores)
    with open(args.output_file, "w") as ostream:
        ostream.write("\t".join(M_SCHEMA.keys()) + "\n")
        for mi, mr_list in enumerate(list_of_matched_recrods):
            target_kmer_name = matched_kmernames[mi]
            for mrow in mr_list:
                ostream.write(f"{target_kmer_name}\t")
                ostream.write("\t".join(map(str, mrow)) + "\n")


def target_to_kmers(args, merge_jf_kmer_dir="step3_merge_jf_kmer"):
    sample_name = args.sample_name
    filename = f"{merge_jf_kmer_dir}/{sample_name}.tsv"
    assert os.path.exist(filename), f"{filename} doesn't exist"

    mytable = defaultdict(list)
    with open(filename) as stream:
        next(stream)
        current_target = None
        current_target_hits = None
        for line in stream:
            line = line.strip('\n').split('\t')
            target, kmer_seq = line[0], line[2]
            if target != current_target:
                if current_target_hits is not None:
                    mytable[current_target] = current_target_hits
                current_target = target
                current_target_hits = defaultdict(list)
            # caution: there are kmers forward and reversecomplement are the same
            current_target_hits[kmer_seq] = []
            mytable[kmer_seq][""]
        # for the last target sequences
        if current_target_hits is not None:
            mytable[current_target] = current_target_hits


def generate_target_to_cakmer(args, merge_jf_kmer_dir="step3_merge_jf_kmer"):
    sample_name = args.sample_name
    filename = f"{merge_jf_kmer_dir}/{sample_name}.tsv"
    assert os.path.exist(filename), f"{filename} doesn't exist"

    matched_canonical_kmers_meta = defaultdict(list)
    with open(filename) as stream:
        # skip the header
        next(stream)
        for line in stream:
            line = line.strip('\n').split('\t')
            target, kmer_seq = line[0], line[2]
            matched_canonical_kmers_meta[kmer_seq].append(target)


def kmerize_seq(sequence, k=20):
    # when we search for the reads,
    # it might need to be both forward and reverse complement since we used canonicalized k-mers for the reads
    for start in range(0, len(sequence) - k + 1):
        kmer = sequence[start:start+k]
        # canonicalize the k-mer
        revcomp = screed.rc(kmer)
        if kmer < revcomp:
            yield kmer
        else:
            yield revcomp


def search_canonical_kmers_against_reads(packed_args):
    chunk_id, id_start, id_end, input_fasta_file, flag = packed_args
    global matched_canonical_kmers
    global chunk_tsvs_path

    out_tsvfile = chunk_tsvs_path[chunk_id]
    line_start, line_end = 2*id_start+1, 2*id_end
    cat_cmd = f"cat {input_fasta_file} | awk 'NR >= {line_start} && NR <={line_end}' "
    chunk_in_memory = command_output(cat_cmd, quiet=False).strip('\n').split('\n')
    num_of_reads = 2 * (id_end - id_start)
    assert len(chunk_in_memory) == num_of_reads, f"read in chunk_id {chunk_id} failed"

    with OutputStream(out_tsvfile) as ostream:
        for lineno in range(0, num_of_reads, 2):
            read_header = chunk_in_memory[lineno].strip('\n')
            read_seq = chunk_in_memory[lineno+1].strip('\n')
            read_ca_kmer_set = set(kmerize_seq(read_seq, k=20))
            matched_cakmers = list(matched_canonical_kmers.intersection(read_ca_kmer_set))
            for matched_ca_kmer in matched_cakmers:
                ostream.write(f"{matched_ca_kmer}\t{read_header}\n")
                # one read could have multiple mathched canonical kmers
    return "it worked"


def query_matched_reads(args):
    sample_name = args.sample_name
    output_dir = args.output_dir

    # Read in merged matched neighbor kmers from jellyfish
    global matched_canonical_kmers
    matched_canonical_kmers = []
    with open(args.merged_kmers_file) as stream:
        next(stream)
        for line in stream:
            line = line.strip('\n').split('\t')
            target, kmer_seq = line[0], line[2]
            matched_canonical_kmers.append(kmer_seq)
    matched_canonical_kmers = set(matched_canonical_kmers)

    # Now I need to scan the reads: one read may include more than one candidate kmers
    sub_reads_file = args.subset_reads_fasta
    total_read_counts = int(command_output(f"grep -c '@' {sub_reads_file}", quiet=False))

    global chunk_tsvs_path
    chunk_tsvs_path = []
    command(f"mkdir -p {output_dir}/temp/{sample_name}")

    # multiprocessing
    chunk_size = args.chunk_size
    chunk_id = 0
    number_of_chunks = ceil(total_read_counts/chunk_size) - 1
    arguments_list = []
    for ni, ci in enumerate(range(0, total_read_counts, chunk_size)):
        headerless_path = f"{output_dir}/temp/{sample_name}/chunkid_{chunk_id}.tsv"
        if ni == number_of_chunks:
            slice_args = (chunk_id, ci, total_read_counts, sub_reads_file, True)
        else:
            slice_args = (chunk_id, ci, ci+chunk_size, sub_reads_file, False)
        arguments_list.append(slice_args)
        chunk_tsvs_path.append(headerless_path)
        chunk_id += 1

    results = multiprocessing_map(search_canonical_kmers_against_reads, arguments_list, args.num_cores)
    assert all(s == "it worked" for s in results)

    # Write the neighbor-reads search results to file
    out_file = f"{output_dir}/{sample_name}.tsv"
    with OutputStream(out_file) as stream:
        stream.write(f"ca_kmer\treads\n")
    cat_files(chunk_tsvs_path, out_file, 20)


def compute_primer_hamming_dist(args):
    # Read in sc2 primers fasta
    primers = {}
    for rec in Bio.SeqIO.parse(args.primer_fasta_file, 'fasta'):
        primers[rec.id] = str(rec.seq)

    with open(args.output_file, "w") as ostream:
        ostream.write("\t".join(["reads", "hd", "primer", "position"]) + "\n")
        for header, seq, qual in FastqGeneralIterator(args.reads_fastq):
            # for each read
            close_primer = (None, None)
            close_dist = 10000
            for k, v in primers.items():
                if len(seq) >= len(v):
                    fwd_primer_dist = hamming_distance(seq[:len(v)], v)
                    rev_primer_dist = hamming_distance(seq[-len(v):], v)
                    if fwd_primer_dist < close_dist:
                        close_dist = fwd_primer_dist
                        close_primer = (k, "fwd")
                    if rev_primer_dist < close_dist:
                        close_dist = rev_primer_dist
                        close_primer = (k, "rev")
            ostream.write("\t".join([header, str(close_dist), close_primer[0], close_primer[1]]) + "\n")


def main():
    p = argparse.ArgumentParser(prog="python isscr.py", description='insilico screening helper scripts.')
    p.add_argument(
        "--sample_name", type=str,
        help="sample name") #required=True,
    p.add_argument(
        "--output_file",
        help="Output pileup file")
    p.add_argument(
        "--jellyfish_ctdir", type=str,
        help="path to jellyfish kmer count results")
    p.add_argument(
        '--parse_matched_kmers', action='store_true', default=False,
        help=f"parse_matched_kmers.")
    p.add_argument(
        '--query_matched_reads', action='store_true', default=False,
        help=f"parse_matched_kmers.")
    p.add_argument(
        '--compute_primer_hamming_dist', action='store_true', default=False,
        help=f"compute_primer_hamming_dist.")
    p.add_argument(
        '--num_cores',
        dest='num_cores',
        type=int, default=8,
        help=f"Number of physical cores to use (8)")
    p.add_argument(
        '--chunk_size', dest='chunk_size', type=int, default=50000,
        help=f"per chunk size (50000) when reading in fasta reads files")
    p.add_argument(
        "--candidate_spacers_file",
        type=str,
        help="path to candidate_spacers.txt or guides.fasta")
    p.add_argument(
        "--merged_kmers_file",
        type=str,
        help="path to merged jellyfish search results per sample")
    p.add_argument(
        "--subset_reads_fasta",
        type=str,
        help="path to subseted fasta file")
    p.add_argument(
        "--output_dir",
        type=str,
        help="path to output directory")
    p.add_argument(
        "--reads_fastq", #required=True,
        type=str,
        help="path to preprocessed reads as input")
    args = p.parse_args()

    global candidates
    if args.candidate_spacers_file.endswith(".txt"):
        candidates = read_candidate_spacers(args.candidate_spacers_file)
    else:
        candidates = {f"{record.id}_hd.4":str(record.seq) for record in SeqIO.parse(args.candidate_spacers_file, "fasta")}

    if args.parse_matched_kmers:
        parse_matched_kmer(args)

    if args.query_matched_reads:
        query_matched_reads(args)

    if args.compute_primer_hamming_dist:
        compute_primer_hamming_dist(args)

main()
