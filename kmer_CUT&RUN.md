## Reference free kmer analysis for CUT&RUN data

The purpose of this analysis is to process CUT&RUN data for identifying kmers enriched in CENPA binding sites.

### 0. Data Location and Set Up

Celine's cut&run data:
http://public.gi.ucsc.edu/~khmiga/celine/

Project IDs:
http://public.gi.ucsc.edu/~khmiga/celine/Demultiplexing_Stats.csv

snakemake pipeline for kmer analysis
https://github.com/straightlab/xenla-cen-dna-paper

Working directory on courtyard
```
/public/groups/migalab/projects/hprc_censat_kmer_cenpa
```

Commands I used to set up conda environment
```
source /public/groups/migalab/programs/miniconda/etc/profile.d/conda.sh

conda create -n snake_kmer
conda activate snake_kmer

conda install -c conda-forge mamba
conda install -c bioconda kmc bbmap
conda install pandas
```
Location of snakemake installation:
```
/public/groups/migalab/programs/miniconda/bin/snakemake
```

Use these commands every time you run the pipeline, to activate the correct environment
```
conda activate snake_kmer
export PATH="$PATH:/public/groups/migalab/programs/miniconda/envs/snakemake/bin"

```

### 1. Workflow overview

- Run fastqc on input data to check quality
- Then, run the k-mer snakemake pipeline. The steps in the pipeline are:
  - Takes unprocessed, zipped fastqs as input. Removes PCR duplicates using clumpify.sh and trims adapters using trimmomatic.

  - Processed reads are then used as input for KMC which generates a k-mer database for each fastq. K-mer databases are exported to text files as well as fasta files.

  - The abundance of each k-mer in each k-mer database is then counted and normalized to the number of basepairs in the processed fastq file.

  - An enrichment ratio is then calculated for each k-mer by dividing the number of times it is found in dataset1 by the number of times it is found in dataset2 (where dataset1 and dataset2 are specified in the config file as 'pairing').

  - K-mers are then filtered based on enrichment value cutoff (specified in config). Cutoffs are based on the number of median absolute deviations from the median enrichment value for the given pairing. K-mers with enrichment values above the specified cutoff are written to a fasta file.

  - Filtered k-mers are then used to select reads from a specified fastq (in config file: 'your_fav_seqdata'). Reads containing at least one k-mer from input k-mer file are outputted to a fasta file. This file is then converted into a k-mer table file where each row is a read and each column is a k-mer ID number. Each line contains either a 1 or 0 denoting the presence or absence of the corresponding k-mer in that read.

config file:
```
/public/groups/migalab/programs/Trimmomatic-0.39/trimmomatic-0.39.jar
```

Questions:
- I'm unclear from looking at the README and the paper what the input to `your_fav_seqdata` should be in the config file. Is this an assembly fasta file for the entire genome? the paper mentions the input genome assembly is first separated into 50kb segments, so do I need to do this outside of the snakemake, and then pass into `your_fav_seqdata` the 50kb "chunked" fa file?
- I'm also unclear about the differences between the two sets of config files and snakefiles in the github repo. It looks like they are almost interchangeable except for different example input files and additional values passed in to `KMER_LENS` `CIVAL` and `NUM_MADS` but I wanted to be sure I was using the correct one since there are two.
