# CDR_detect
Contains methods for detecting CDRs in methyl-tagged long read data

## How to run CDR_detect.py

Requirements:
```
- pysam
- numpy
- Bio.seq
- argparse
```

Command:
```

```

Description of method:
- takes in bamfile of reads aligning to the HOR regions
- these reads have methylation tags giving the percent likelihood of methylation at each CpG site on the read
- First, we calculate the average methylation probability across the whole read. If it is < 40%, we automatically call that read a CDR and move on
- For the other reads, we move in a sliding window across each read, by 1 bp. Optimal window size right now is 3000bp
- for each window, we take the average of the methylation probabilities in that window. If any window drops below 30%, we exit and call that read a CDR
- 
