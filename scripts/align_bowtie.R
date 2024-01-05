##################################################
### count bowtie alignments

library(optparse)
library(here)

option_list <- list(make_option(c("-g", "--genome"), type="character", default=NULL,
                                help="bowtie index prefix", metavar="character"),
                    make_option(c("-m", "--mismatch"), type="character", default=1,
                                help="number of mismatches allowed", metavar="integer"),
                    make_option(c("-e", "--enzyme"), type="character", default=NULL, # "Cas13a" or "Cas12"
                                help="Cas enzyme type", metavar="character"),
                    make_option(c("-v", "--omit"), type="character", default=NULL,
                                help="transcripts to omit", metavar="character"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"),
                    make_option(c("-b", "--bowtie"), type="character", default=NULL,
                                help="path to bowtie", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

bowtie_path <- system("which bowtie", intern=T)
if(length(bowtie_path)==0 & is.null(opt$bowtie)) {
  cat("ERROR: need to supply path to bowtie")
  q(save="no")
}

if(is.null(opt$genome)) {
  cat("\nERROR: no reference genome specified")
  q(save="no")
} else {
  if(!("enzyme" %in% names(opt))) {
    cat("\nERROR: no enzyme specified\n")
    q(save="no")
  } else {
    cat(paste("\nAligning against off-target:", opt$genome, "(bowtie) \n"))
  }
}

# align windows with bowtie
cts_fname <- file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam"))
if(!file.exists(cts_fname)) {
  cat(paste("- aligning windows to", opt$genome, "\n"))
  system(paste(bowtie_path,
               ifelse(opt$enzyme=="Cas13a", "--norc", ""), # for Cas13a, do not align to reverse complement
               "-k 50", # report up to 50 alignments
               "-v", opt$mismatch, # up to opt$mismatch mismatches allowed
               "-S", # output as .sam alignment file
               "--un", file.path(opt$out, paste0("bowtie_", opt$genome, "_unmapped.fa")), # fasta file of unmapped windows
               "-f", file.path(here(), "ref_data", opt$genome),  # path to bowtie index
               file.path(opt$out, "targets.fa"), # fname of windows fasta file
               ">", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam")), # sam alignment file of mapped windows
               "2>", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.bowtiestats")))) # fname of bowtie output
}

# read in alignments
unmapped_fname <- file.path(opt$out, paste0("bowtie_", opt$genome, "_unmapped.fa"))
if(file.exists(unmapped_fname)) {
  unmapped <- read.table(unmapped_fname, stringsAsFactors=F)$V1
  unmapped <- unmapped[!grepl(">", unmapped)]
} else {
  unmapped <- c()
}
mapped <- system(paste("grep -v ^@", file.path(opt$out, paste0("bowtie_", opt$genome, "_mapped.sam")),
                          "| cut -f1,2,3,10"), intern=T)
mapped <- data.frame(matrix(unlist(strsplit(mapped, split="\t")), ncol=4, byrow=T), stringsAsFactors=F)
colnames(mapped) <- c("window", "flag", "aligned_to", "template")
mapped <- subset(mapped, mapped$aligned_to != "*")
if(opt$enzyme == "Cas13a") {
  mapped <- subset(mapped, flag != 16) # do not report alignments to minus strand
}
if(!is.null(opt$omit)) {
  omit_transcripts <- readLines(opt$omit)
  mapped <- subset(mapped, !(aligned_to %in% omit_transcripts))
}

# format table of alignment counts
cat("- counting alignments\n")
windows <- read.table(file.path(opt$out, "windows.txt"), header=T, stringsAsFactors=F)
window_cts <- data.frame(seq=windows$target,
                         segment=windows$segment,
                         start=windows$start,
                         strand=windows$strand,
                         ct=sapply(windows$target,
                                   function(x) {
                                     if(x %in% unmapped) {
                                       return(0)
                                     } else {
                                       tmp <- subset(mapped, mapped$template==x)
                                       return(length(unique(tmp$aligned_to)))
                                     }
                                   }))
cat("- number of alignments:\n")
print(data.frame("number of alignments"=0:5,
                 count=sapply(0:5, function(x) sum(window_cts$ct==x))))

# write alignment counts to output
write.table(window_cts,
            file=file.path(opt$out, paste0("alignment_cts_", opt$genome, ".txt")),
            quote=F, sep="\t", row.names=F)
