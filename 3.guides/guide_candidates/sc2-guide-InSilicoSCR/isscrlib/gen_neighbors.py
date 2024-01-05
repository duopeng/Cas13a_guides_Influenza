#!/usr/bin/env python3

import os, sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(dna):
    """ Find the reverse complement of a DNA string """
    dnadict = {'A':'T','C':'G','G':'C','T':'A'}
    reverseDna = [ dnadict[c] for c in dna ]
    return reverseDna[::-1]


def hamming_distance(str1,str2):
    """ Compute the Hamming distance between two strings """
    hd = 0
    if len(str1) != len(str2):
        print("the two strings different length ERROR")
        return None
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            hd += 1
    return hd


def immediateNeighbor(pattern):
    """ 1 edit distant away neighbor from the pattern """
    neighbor = [pattern]
    for i in range(len(pattern)):
        symbol = pattern[i]
        for x in ['A', 'C', 'G', 'T']:
            if x != symbol:
                neighbor.append(pattern[:i] + x + pattern[i+1])
    return neighbor


def generate_neighbors(pattern, d):
    """ Generate the d-neighborhood, the set of all k-mers with hamming dist not exceed d """
    if d == 0:
        return pattern

    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighbors = []
    suffixNeighbors = generate_neighbors(pattern[1:], d)
    for text in suffixNeighbors:
        if hamming_distance(pattern[1:], text) < d:
            for x in ['A', 'C', 'G', 'T']:
                neighbors.append(x + text)
        else:
            neighbors.append(pattern[:1] + text)
    neighbros = list(set(neighbors))
    return neighbors


def convert_guide_to_target(guides_fasta_fp):
    """ convert RNA guide sequences to cDNA target sequences"""
    for record in SeqIO.parse(guides_fasta_fp, "fasta"):
        print(">" + record.id)
        my_rna = Seq(str(record.seq))
        #print(my_rna)
        my_target = my_rna.back_transcribe().reverse_complement()
        #print(my_rna.back_transcribe())
        print(my_target)


def main():
    """ candidate_spacers 6969 """
    # Usage: grep -v "start" candidate_spacers.txt | xargs -Ixx -P 48 bash -c "python3.7 gen_neighbors.py xx"
    if len(sys.argv) == 5:
        # 6969 candidate_spacers
        start_pos = sys.argv[1]
        input_kmer = sys.argv[2]
        strand = sys.argv[4]
        sl = 1 if strand == "+" else 0

        d_neighbors = generate_neighbors(input_kmer, 4)

        outfile = f"kmer_targets/sp.{start_pos}_strand.{sl}_hd.4.txt"
        with open(outfile, "w") as stream:
            for line, kmer in enumerate(d_neighbors):
                hd = hamming_distance(kmer, input_kmer)
                #stream.write("\n".join(d_neighbors) + "\n")
                stream.write(f">{line}_hd{hd}\n{kmer}\n")

    """ given target.fasta """
    if len(sys.argv) == 3:
        target_fastafile = sys.argv[1]
        outdir = sys.argv[2]
        for record in SeqIO.parse(target_fastafile, "fasta"):
            seq_name = record.id
            input_kmer = str(record.seq)
            d_neighbors = generate_neighbors(input_kmer, 4)
            outfile = f"{outdir}/{seq_name}_hd.4.txt"
            with open(outfile, "w") as stream:
                for line, kmer in enumerate(d_neighbors):
                    hd = hamming_distance(kmer, input_kmer)
                    stream.write(f">{line}_hd{hd}\n{kmer}\n")

main()
