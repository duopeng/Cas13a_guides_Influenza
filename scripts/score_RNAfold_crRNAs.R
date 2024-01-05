##################################################
### compute crRNA folding structures

library(optparse)

option_list <- list(make_option(c("-e", "--enzyme"), type="character", default=NULL, # "Cas13a" or "Cas12"
                                help="Cas enzyme type", metavar="character"),
                    make_option(c("-c", "--cas_repeat"), type="character",
                                default="uagaccaccccaaaaaugaaggggacuaaaac",
                                help="crRNA repeat sequence", metavar="character"),
                    make_option(c("-s", "--spacer"), type="character", default=NULL,
                                help="example `good` spacer", metavar="character"),
                    make_option(c("-o", "--out"), type="character", default=".",
                                help="output directory", metavar="character"),
                    make_option(c("-r", "--rnafold"), type="character", default=NULL,
                                help="path to rnafold", metavar="character"))
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

cat("\nRNAfold score: crRNA structure\n")

rnafold_path <- system("which RNAfold", intern=T)
if(is.null(opt$rnafold)) {
  opt$rnafold <- rnafold_path
} else {
  if(length(rnafold_path)==0) {
    cat("ERROR: need to supply path to rnafold")
    q(save="no")
  }
}

# establish repeat and guide sequences
if(is.null(opt$spacer)) {
  if(opt$enzyme == "Cas13a") {
    # cas_repeat <- "UAG ACC AGC CCA AAA AUG AAG GGC ACU AAA AC" # from Gavin
    # cas_repeat <- gsub(" ", "", cas_repeat)
    # cas_repeat <- toupper("uagaccaccccaaaaaugaaggggacuaaaac") # from Gavin
    cas_repeat <- gsub(" ", "", toupper(opt$cas_repeat))
    good_spacer <- "GCA GCG CCU CUU GCA ACG AU" # from Gavin
    good_spacer <- gsub(" ", "", good_spacer)
  } else {
    if(opt$enzyme == "Cas12") {
      # cas_repeat <- toupper("aauuucuacuaaguguagau") # from Gavin
      cas_repeat <- gsub(" ", "", toupper(opt$cas_repeat))
      good_spacer <- "UCGAUGGGGAAACCUUACCCUCCAG" # from Creutzburg et al. NAR (2020) - Figure 1A
    }
  }
}

# fold good crRNA
good_guide_structure <- system(paste("echo", paste0(cas_repeat, good_spacer), "| ", rnafold_path, "--noPS"), intern=T)[2]
good_guide_structure <- strsplit(good_guide_structure, split=" ")[[1]][1]
hairpin_structure <- substr(good_guide_structure, start=1, stop=nchar(cas_repeat))

# fold windows
cat("- computing folding energies (crRNAs)\n")
spacers <- readLines(file.path(opt$out, "spacers.txt"))
crRNA_seq <- paste0(cas_repeat, spacers)
writeLines(crRNA_seq, con=file.path(opt$out, "crRNAs.txt"))
system(paste0(rnafold_path, " -i ", file.path(opt$out, "crRNAs.txt"),
              " --noPS --outfile=crRNAs_RNAfold.txt"))

# read in RNAfold output: windows
RNAfold <- data.frame(matrix(readLines(file.path("crRNAs_RNAfold.txt")),
                             ncol=2, byrow=T), stringsAsFactors=F)
colnames(RNAfold) <- c("sequence", "structure")
RNAfold$MFE <- as.numeric(sub(" ", "", sub("\\)", "", sub(".*\\(", "", RNAfold$structure))))
RNAfold$structure <- sub(" .*", "", RNAfold$structure)

# flag if hairpin has been interfered
RNAfold$has_hairpin <- (substr(RNAfold$structure, start=1, stop=nchar(cas_repeat)) == hairpin_structure)

# flag how many base-pairs in spacer
RNAfold$spacer_basepairs <- sapply(RNAfold$structure,
                                   function(x) {
                                     spacer <- substr(x, start=nchar(cas_repeat)+1, stop=nchar(x))
                                     spacer <- strsplit(spacer, split="")[[1]]
                                     return(sum(spacer %in% c("(", ")")))
                                   })

# add start index and strand
RNAfold$segment <- read.table(file.path(opt$out, "windows.txt"), header=T)$segment
RNAfold$start <- read.table(file.path(opt$out, "windows.txt"), header=T)$start
RNAfold$strand <- read.table(file.path(opt$out, "windows.txt"), header=T)$strand

# output RNAfold scores
write.table(RNAfold,
            file=file.path(opt$out, "score_RNAfold_crRNAs.txt"),
            quote=F, sep="\t", row.names=F)
