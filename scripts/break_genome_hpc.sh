#!/bin/bash

#SBATCH --job-name=Cas13a_gRNA
#SBATCH --time=14-00:00:00
#SBATCH --partition cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=256G
#SBATCH --cpus-per-task=32
#SBATCH -e slurm-%A_%a.err
#SBATCH -o slurm-%A_%a.out

enzyme=Cas13a
input_prefix="Orthomyxoviridae.fna.Influenza.fa_chunk_"
input_suffix=".fasta"
outputdir="tiled_genomes"
script_path="${2}"

python "${script_path}/break_genomes.py" --enzyme "${enzyme}" \
    --scriptdir "${script_path}" \
    --input "${input_prefix}${1}${input_suffix}" \
    --out "${outputdir}/${enzyme}_Influenza_chunk_${1}" \
    --num_workers 32