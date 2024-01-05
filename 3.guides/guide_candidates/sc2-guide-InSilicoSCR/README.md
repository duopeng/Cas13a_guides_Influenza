# InSilico SCReening for SC2 guides

Please refer to [wiki](https://github.com/czbiohub/sc2-guide-InSilicoSCR/wiki/InSilicoSCR) for details of the pipeline.


## Create the conda environment `isscr`

```
conda env update --name isscr -f env.yml
```

## Compuate exact match via grep

```
paste - - < guides.fa  | awk '{print $2, substr($1,2)}' | xargs -Ixx -P 8 bash -c "bash grep_kmer_counts.sh xx"
```
