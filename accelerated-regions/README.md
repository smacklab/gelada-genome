# Gelada accelerated region pipeline

This pipeline works off a modified version of [maffilter](https://github.com/jydu/maffilter). maffilter includes many algorithms, including two—`AlnFilter` and `EntropyFilter`—that work by applying filtering criteria to sliding windows. Normally, windows that fail the criteria are removed and blocks are split, with blocks encompassing failed coordinates optionally being written to a "trash" file. The modified maffilter used in this pipeline instead writes coordinates of windows that pass filters to the "trash" file. The filtering criteria across windows thus emphasizes keeping alignment tracts that pass filters in _any_ sliding window rather than keeping alignment tracts that pass filters in _all_ sliding windows.

## Configure

To configure this pipeline, edit the file in `config/settings.cfg`

## Make

In order to build the modified version, run the following:

```
mkdir -p out
sbatch slurm/make-maffilter.sh
```

## Download maf files

```
scripts/list-maf.sh
sbatch slurm/download-maf.sh
```

## Process MAFs

```

sbatch slurm/filter-species.sh
sbatch --array=1-$(wc -l data/maf-pass.txt | xargs | cut -d ' ' -f 1) slurm/preprocess-maf-single.sh
sbatch --mem=16gb --time=24:00:00 --array=1-$(wc -l data/maf-pass.txt | xargs | cut -d ' ' -f 1) slurm/filter-gaps-single.sh 50 90

sbatch --exclusive slurm/filter-entropy.sh 50 90 90

sbatch --array=1-$(wc -l data/maf-pass.txt | xargs | cut -d ' ' -f 1) slurm/extract-blocks-single.sh 50 90 90

sbatch slurm/postprocess-maf.sh 50 90 90

sbatch slurm/split-single-blocks.sh 50 90 90
```

## Fit phylogenetic models
```
sbatch slurm/fit-model.sh 50 90 90
sbatch slurm/fit-cons.sh 50 90 90

sbatch slurm/fit-phylop.sh 50 90 90

sbatch slurm/block-stats.sh

scripts/summarize-cons.R 50 90 90
scripts/summarize-stats.R 50 90 90
scripts/summarize-phylop.R 50 90 90
scripts/summarize-all.R 50 90 90
```
