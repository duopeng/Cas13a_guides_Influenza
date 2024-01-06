# Cas13a guides for Influenza detection

This workflow is used to predict Cas13a guides for Influenza detection. 
- Analyzes ***all*** 20nt windows in ***all*** influenza sequences (downloaded from bbvrc.org in this workflow).  
- Compute the number of segment, strain and subtype each guide can target.
- Select guides that target the most segments, strains or subtypes.
- Analyze crRNA folding structures using RNAfold.
- Check cross-reactivity against human transcriptome and common microbes.
- The computation pipline is divided into three stages:  
    1. download Influenza genomes
    2. preprocess Influenza genomes  
    3. analyze Cas13a guides  

<br>

# Setup

### clone repository
```shell
git clone https://github.com/duopeng/Cas13a_guides_Influenza
```

### create conda environment
```shell
conda env create -f environment.yml
```

### install kraken2
Follow  [instructions](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#installation) to install Kraken2.   
Add executables `kraken2` and `kraken2-build` to the PATH environment variable.  
(Kraken2 is used to check guide's cross-reactivity with microbes) 

<br>

# Guide prediction pipeline stage 1: download and process Influenza genomes

### download Orthomyxoviridae genomes/sequences
note that for the downloads, we use genomes and sequences interchangeably.

```shell
cd 1.download_influenza_genomes
wget -N "ftp://ftp.bvbrc.org/viruses/Orthomyxoviridae.fna"
```

### 1.1 extract influenza genomes

```shell
module load anaconda # hpc-specific command to activate conda env
conda activate Cas13 # hpc-specific command to activate conda env
# in folder "1.download_influenza_genomes"
python ../scripts/parse_fasta.py -f Orthomyxoviridae.fna -k Influenza
```

### 1.2 split the genome file into 10 chunks (for parallelization)

```shell
n=10
# in folder "1.download_influenza_genomes"
python ../scripts/split_fasta.py -f Orthomyxoviridae.fna.Influenza.fa -n $n
```

<br>

# Guide prediction pipeline  Stage 2: preprocess Influenza genomes
### 2.1 break genomes into tiling windows of 20nt (guides length)

(submitting 10 slurm jobs via sbatch, takes ~1hr on hpc)  
Note that ***all*** 20nt windows in ***all*** sequences are analyzed and saved to `tiled_genomes` folder

```shell
# in folder "1.download_influenza_genomes"
for ((i=1;i<=n;i++))
do
    sbatch ../scripts/break_genome_hpc.sh "$i" "../scripts"
    echo "submitted job $i"
done
```

### 2.2 parse strains, subtypes and segments
From description of sequences, extract strain, subtype and segment information.  
produce output file: `Influenza_sequence_description_parsed.txt` in directory `output`  
**This step is performed in Jupyter notebook**: `analyze_influenza_genomes.ipynb` in directory `2.analyze_influenza_genomes`
    
    

<br>

# Guide prediction pipeline stage 3: analyze Cas13a guides
### 3.1 get guides table
This step (3.1):
- analyzes ***all*** 20nt windows in ***all*** sequences (enumerated in previous steps) 
- Compute segment, strain and subtype coverage for each guide
- Result table is saved to `guides_table.tsv` in directory `3.guides`
- takes ~3hrs with multiprocessing, num_workers=10
```shell
# in repository root directory
folder_path_prefix="1.download_influenza_genomes/tiled_genomes/Cas13a_Influenza_chunk_"
genome_df="2.analyze_influenza_genomes/output/Influenza_sequence_description_parsed.txt"
python scripts/get_guides_table.py --folder_path_prefix "${folder_path_prefix}" --genome_df "${genome_df}" --outputdir "3.guides"
```
### 3.2 get top guides
get top guides and produce four outut files (for subsequent analyses): targets.fa windows.txt targets.txt and spacers.txt

```shell
# in repository root directory
python scripts/top_guides.py --input "3.guides/guides_table.tsv" --outputdir "3.guides/guide_candidates" --metric "strain_targeted_count" --cutoff 2000
```

### 3.3 calculate crRNA folding structures 
script `score_RNAfold_crRNAs.R` taken from https://github.com/lareaulab/cas13a_guide_design  
output files are in folder `3.guides/guide_candidates/crRNAs_RNAfold`

```shell
module load R/4.3
cd 3.guides/guide_candidates # execute R script in folder "3.guides/guide_candidates"
mkdir -p crRNAs_RNAfold 
rm crRNAs_RNAfold/crRNAs_RNAfold.txt 
Rscript ../../scripts/score_RNAfold_crRNAs.R --enzyme Cas13a --rnafold RNAfold --out .
mv crRNAs.txt crRNAs_RNAfold.txt score_RNAfold_crRNAs.txt crRNAs_RNAfold/
```

### 3.4 check cross-reactivity against hg38 transcriptome

get human latest transcriptome and build bowtie index
```shell
cd ..
# in 3.guides/guide_candidates
mkdir -p check_alignment_to_hg38 && cp targets.fa windows.txt check_alignment_to_hg38/
mkdir -p refdata && cd refdata
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz
gunzip GRCh38_latest_rna.fna.gz
bowtie-build GRCh38_latest_rna.fna GRCh38_latest_rna
```

align guides against human transcriptome
```shell
# in 3.guides/guide_candidates/check_alignment_to_hg38
rm bowtie_GRCh38_latest_rna_mapped.sam
Rscript ../../../scripts/align_bowtie.R --genome GRCh38_latest_rna --enzyme Cas13a --out .
```

### 3.5 check cross-reactivity against common microbes
prepare microbial genomes database (kraken2 db), takes a ~22 hours  
logs will be written to file `kraken2_std_db_build.log`
```shell
module load gcc #required by kraken2
kraken2-build --standard --threads 64 --kmer-len 20 --minimizer-len 20 --minimizer-spaces 5 --db kraken2_standard_db &> kraken2_std_db_build.log
```

generate guide neighbors (sequences with 1, 2, 3, 4 mismatches), takes a few hours
```shell
conda activate /hpc/mydata/duo.peng/anaconda/2023.03/x86_64/envs/isscr
cd ../sc2-guide-InSilicoSCR && cp ../targets.fa ./
# in folder sc2-guide-InSilicoSCR
snakemake --configfile=config.yml --cores 32 -p generate_neighbors
```

run kraken2 on guide neighbors (~10min)
```shell
cd ..
DBNAME="kraken2_db/kraken2_standard_db/"
# in folder 3.guides/guide_candidates
python ../../scripts/run_kraken2.py --num_workers 32 --kraken2_db $DBNAME --folder_path sc2-guide-InSilicoSCR/1_neighbors/
# parse kraken2 results, and check LCAs of all 4-neighbors, generate hd_LCA.counts.txt
python ../../scripts/parse_Kraken2_res.py --num_workers 64 --kraken2_res_folder sc2-guide-InSilicoSCR/1_neighbors_kraken2.out/
```

### 3.6. combine cross-reactivity and fold results

```shell
# in folder 3.guides/guide_candidates
python ../../scripts/combine_filter_results.py \
--num_workers 64 --df_guides df_top_guides.tsv \
--kraken2_res_folder sc2-guide-InSilicoSCR/1_neighbors_kraken2.out \
--RNAfold_res_path crRNAs_RNAfold/score_RNAfold_crRNAs.txt \
--bowtie_hg_res_path check_alignment_to_hg38/alignment_cts_GRCh38_latest_rna.txt
```

## Final guide table
contains cross-reactivity and folding results  
`df_top_guides_filters_res_combined.tsv` in directory guide_candidates