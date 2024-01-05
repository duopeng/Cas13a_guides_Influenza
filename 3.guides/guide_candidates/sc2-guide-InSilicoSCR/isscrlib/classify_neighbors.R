## Chunyu Zhao 2020-06-04

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(readr)
library(magrittr)
library(data.table)
library(seqinr)
library(stringdist)



########### Input arguments
guide_name <- args[1]
kq_clas_reads_fp <- args[2]
output_file <- args[3]
kq_db_dir <- args[4]
guides_fasta <- args[5]


########### Read in KrakenUniq results
# classified_kmers from krakenuniq
classified_kmers <- read_delim(kq_clas_reads_fp, delim="\t", col_names = c("label", "kmer_name", "taxid", "length", "classification")) %>% select(kmer_name, taxid) %>% mutate(taxid = as.character(taxid))

# Read in corresponding k-mer sequence
input_kq_dir <- dirname(kq_clas_reads_fp)
kmer_list <- read.fasta(file.path(input_kq_dir, "classified_out.fasta"), as.string = TRUE, forceDNAtolower=F)
kseqs = unlist(rbind(lapply(1:length(kmer_list), function(x) getSequence(kmer_list[[x]], as.string = T)[[1]])))
knames = unlist(rbind(lapply(1:length(kmer_list), function(x) getAnnot(kmer_list[[x]], as.string = T))))
kmer_fasta <- data.frame(kmer_name = knames, kmer_seq = kseqs) %>% mutate(kmer_name = gsub(">", "", kmer_name))
# Add k-mer sequence
classified_kmers %<>% left_join(kmer_fasta, by = c("kmer_name"))


########### Read in KrakenUniq seqid2taxid mapping
# Add seq_id (aka contig accession name) (from fasta header)
seqid2taxid <- read_delim(file.path(kq_db_dir, "seqid2taxid.map"), delim="\t", col_names = c("seqid", "taxid"), col_types = list(col_character(), col_character()))
classified_kmers %<>% left_join(seqid2taxid, by=c("taxid"))

# Add taxName
taxDB <- read_delim(file.path(kq_db_dir, "taxDB"), delim="\t", col_names = c("taxid", " parent_taxid", "taxName", "rank")) %>% mutate(taxid = as.character(taxid))
classified_kmers %<>% left_join(taxDB, by=c("taxid" = "taxid"))


########### Compuate hamming distance
if (FALSE) {
  # This works for the 6969 guides list
  candidate_targets <- read_delim(file.path(project_dir, "candidate_spacers.txt"), delim = "\t") %>%
    mutate(strand_label = ifelse(strand == "+", 1, 0)) %>%
    mutate(kmer_fname = paste("sp.", start, "_strand.", strand_label, "_hd.4", sep=""))
  target_kmer_seq <- candidate_targets %>% filter(kmer_fname %in% guide_name) %>% .$target
}


guides <- read.fasta(guides_fasta, as.string = TRUE, forceDNAtolower=F)
kseqs = unlist(rbind(lapply(1:length(guides), function(x) getSequence(guides[[x]], as.string = T)[[1]])))
knames = unlist(rbind(lapply(1:length(guides), function(x) getAnnot(guides[[x]], as.string = T))))
curr_guide_seq <- data.frame(kmer_name = knames, target_kmer = kseqs) %>% mutate(kmer_name = gsub(">", "", kmer_name)) %>% mutate(target_kmer = as.character(target_kmer)) %>% filter(kmer_name %in% guide_name) %>% .$target_kmer

classified_kmers %<>% mutate(hamming_dist = stringdist(curr_guide_seq, kmer_seq, method="hamming")) %>%
  select(kmer_name, kmer_seq, hamming_dist, taxid, everything())


########### Write classified_kmers to file
classified_kmers %>% arrange(hamming_dist) %>% write.table(output_file, quote=F, sep="\t", row.names = F)
