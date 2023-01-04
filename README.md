# CDR_detect
Detecting CDRs in long reads with methylation tags

```
|__ scripts
    |__ CDR_detect_reads.py returns readnames that contain potential CDR (older version)
    |__ CDR_detect.py returns readnames and coordinates of CDR in the reads
    |__ evaluate_CDR_detect.py calculates P,R,F1 against input truth CDR reads
    |__ CDR_detect.wdl wdl workflow for hprc assemblies
    |__ Dockerfile used for docker in CDR_detect.wdl
|__ conference_summary.md summary of project July 2022
|__ CDR_detect_notes.md
```

## How to run CDR_detect.py

Python libraries needed:
```
- pysam
- numpy
- Bio.seq
```

Command:
```
python3 CDR_detect.py -b bamfile -o outfile -w <windowsize> -x <methylation threshold> -n <minimum CDR size>
```

Inputs:
- `-b` (required) bamfile of reads from the alpha satellite HOR regions. These reads have methylation tags giving the percent likelihood of methylation at each CpG site on the read
- `-o` (required) output text file name
- `-w` (default = 1000 bp) Window size of CDRs to consider in bp
- `-x` (default = 0.3) methylation frequency threshold for calling a CDR in a window
- `-n` (default = 3000 bp) minimum CDR size in bp


Description of method:

- For each read, move forward in a sliding window of size `w` by 1 bp
- Calculate methylation frequency in that window
- If the methylation frequency drops below `x`, record current window start coordinate as start coordinate of a CDR
- Once the methylation frequency rises above `x`, record current window end coordinate as end coordinate of CDR
- Return readnames and the coordinates of predicted CDRs inside them. Only return CDRs > `n` (size in bp)
