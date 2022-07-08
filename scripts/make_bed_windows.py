# prints bedfile with windows of binsize s
# python3 make_bed_windows.py -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.nonCDR.CDR.new.pb_cpg.combined.reference.bed -s 1000

import argparse
import numpy as np
import time

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='smooth_methyl.py',
        description="""smooth methylation data into bins of size x""")

    parser.add_argument("-b", "--bed",
                        required=True,
                        metavar="input.bed",
                        help="Bedfile in format output by https://github.com/PacificBiosciences/pb-CpG-tools")
    parser.add_argument("-s", "--binsize",
                        required=True,
                        metavar="binsize",
                        help="binsize to smooth windows")

    return parser.parse_args()

def get_genome(bedPath):
    '''
    retrieves start and end coords for each contig in input bed file
    basically created .fai from a bedfile
    '''
    # make genome file - one line per contig with total genome size
    data = []
    genome = {}
    with open(bedPath, "r") as inFile:
        lines = inFile.readlines()
        for line in lines:
            row = [i for i in line.strip().split("\t")]
            data.append(row)
    data = np.array(data)
    contigs = set(data[:, 0])
    for c in contigs:
        starts = data[data[:, 0] == c, 1]
        ends = data[data[:, 0] == c, 2]

        genome[c] = (int(min(starts)), int(max(ends)))

    return genome

def make_windows(genome, binsize):
    '''
    args.genome=dictionary where {contig: (start,end)}
    return generator of bed coords with windows of binsize from start to end in each contig
    '''
    windows=[]
    for contig in genome.keys():
        for i in range(genome[contig][0], genome[contig][1]-binsize+1, binsize):
            windows.append((contig,i,i+binsize))

    return windows


def main():
    '''

    '''
    # parse command line arguments
    t0 = time.time()
    args = arg_parser()

    # get chromsizes AKA genome file for bedtools makewindows
    genome = get_genome(args.bed)

    with open(args.bed + ".genome", "w") as gFile:
        for key,value in genome.items():
            print(key, value[0],value[1],sep="\t",file=gFile)

    # make windows based on binsize
    windows=make_windows(genome, int(args.binsize))
    with open(args.bed +"_windows"+args.binsize+".bed", "w") as outFile:
        for line in windows:
            print(*line,sep="\t",file=outFile)

if __name__ == '__main__':
    main()