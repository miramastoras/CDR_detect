#!/usr/bin/env python3

'''
Purpose: Determine whether reads contain a potential CDR
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 CDR_detect.py -b bamfile -o outfile_prefix
'''


import argparse
import pysam
import re
import numpy as np
from Bio.Seq import Seq
import time, sys
import bisect

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='CDR_detect.py',
        description="""Separate candidate CDR containing reads from a bamfile""")

    parser.add_argument("-b", "--bam",
                        required=True,
                        metavar="input.bam",
                        help="The aligned BAM file.")
    parser.add_argument("-o", "--out",
                        required=True,
                        metavar="output text file prefix",
                        help="prefix for output list of CDR read names and non CDR read names")
    parser.add_argument("-w", "--windowsize",
                        required=False,
                        metavar="windowsize in bp",
                        help="Windowsize in bp on each read to consider smoothed methylation percentage")
    parser.add_argument("-x", "--window_threshold",
                        required=False,
                        metavar="methylation threshold to call CDR in a window",
                        help="Percentage methylation below which to call a window containing a CDR")
    parser.add_argument("-t", "--read_threshold",
                        required=False,
                        metavar="read length average methylation threshold to call CDR",
                        help="Percentage methylation below which to call a read containing a CDR")
    parser.add_argument("-i", "--step_size",
                        required=False,
                        metavar="read length average methylation threshold to call CDR",
                        help="Percentage methylation below which to call a read containing a CDR")

    return parser.parse_args()

def get_base_indices(query_seq, base, reverse):
    """
    Reference: https://github.com/PacificBiosciences/pb-CpG-tools/blob/main/aligned_bam_to_cpg_scores.py
    Find all occurrences of base in query sequence and make a list of their
    indices. Return the list of indices.
    :param query_seq: The original read sequence (not aligned read sequence). (str)
    :param base: The nucleotide modifications occur on ('C'). (str)
    :param reverse: True/False whether sequence is reversed. (Boolean)
    :return: List of integers, 0-based indices of all bases in query seq. (list)
    """
    if not reverse:
        return [i.start() for i in re.finditer(base, query_seq)]
    # if seq stored in reverse, need reverse complement to get correct indices for base
    # use biopython for this (convert to Seq, get RC, convert to string)
    else:
        return [i.start() for i in re.finditer(base, str(Seq(query_seq).reverse_complement()))]

def get_mod_sequence(integers):
    """
    Reference: https://github.com/PacificBiosciences/pb-CpG-tools/blob/main/aligned_bam_to_cpg_scores.py
    A generator that takes an iterable of integers coding mod bases from the SAM Mm tags, and yields an iterable of
    positions of sequential bases.
    Example: [5, 12, 0] -> [6, 19, 20]
    In above example the 6th C, 19th C, and 20th C are modified
    See this example described in: https://samtools.github.io/hts-specs/SAMtags.pdf; Dec 9 2021
    :param integers: Iterable of integers (parsed from SAM Mm tag). (iter)
    :return mod_sequence: Iterator of integers, 1-based counts of position of modified base in set of bases. (iter)
    """
    base_count = 0
    for i in integers:
        base_count += i + 1
        yield base_count

def parse_tag(read):
    t0=time.time()
    reverse = read.is_reverse
    mmtag = read.get_tag('Mm')
    mltag = read.get_tag('Ml')
    read_seq = read.query_sequence

    # get 0 based locations of all C sites
    base_indeces = get_base_indices(read_seq,"C",reverse)# convert to generator to speed up

    # parse mmtag for mod C+m positions
    modline = next(x[len("C+m") + 1:] for x in mmtag.split(';') if x.startswith("C+m"))

    # get 1-based position of each C from read e.g., [6, 19, 20] = the 6th, 19th, and 20th C bases are modified in the set of Cs
    mod_sequence = get_mod_sequence((int(x) for x in modline.split(',')))
    # get indeces of modified C's from tag relative to the read ( 0 based )
    potModCs = [base_indeces[i - 1] for i in mod_sequence]

    # convert ml tags to probabilities
    ml_probs = [round(x / 256, 3) if x > 0 else 0 for x in mltag]

    # get indeces of modC calls that pass the 50% threshold, relative to their position in the read
    modCs=[]

    modCs=[potModCs[i] for i in range(len(ml_probs)) if ml_probs[i] >=0.5]
    #for i in range(len(ml_probs)):
        #if ml_probs[i] >= 0.5:
            #modCs.append(potModCs[i])

    return base_indeces, modCs, time.time()-t0

def calc_avg_methyl(read,modCs):
    '''
    Returns average methylation frequency across the entire read
    '''
    # count total number of CpG sites in the read
    totalCpG = len(get_base_indices(read.query_sequence, "CG", read.is_reverse))
    # get number of modified CpGs
    mod_count = len(modCs)

    return (mod_count / totalCpG)

def is_read_CDR(read,modCs,w, x, o, t):
    '''
    :param: read: AlignedSegment object representing the read
    :param: allCs: 0-based list of indeces of all C sites in the read
    :param: modCs: 0-based list of indeces of all methylated C sites in the read
    :param: w:  window size as percent of read length
    :j: number of previous windows to compare to
    :x: percent change from current read to number of previous windows to use as threshold
    :return: a boolean determining whether read is a candidate for containing a CDR or not
    '''
    t0=time.time()
    reason=""
    # get average methylation freq across whole read
    avg_methyl=calc_avg_methyl(read, modCs)
    # if its below 40%, call CDR and exit.
    #print(read.query_name, "average methylation ", avg_methyl)
    if avg_methyl < t :
        #print(read.query_name, "\t","less_than_50")
        return "T",time.time() - t0

    else: # else, do sliding window
        # obtain list of frequencies for every sliding window
        #freqs = get_sliding_window_freqs(read, w, allCs, modCs)

        # loop through all windows
        # get index of all CpG sites in read, relative to read
        totalCpGs = get_base_indices(read.query_sequence, "CG", read.is_reverse)
        readlength = len(read.query_sequence)

        # loop through sliding window, in windows with overlap size defined by o
        for i in range(0,readlength - w + 1, o):
            # get all CpG sites in window and all mod CpG sites in window
            allCpGsright=bisect.bisect_left(totalCpGs,i+w)
            allCpGsleft=bisect.bisect_left(totalCpGs,i) # come back to this, should we add 1?

            if allCpGsright == allCpGsleft: # this means there are no CpG sites in the current window
                continue

            modCsright=bisect.bisect_left(modCs,i+w)
            modCsleft=bisect.bisect_left(modCs,i) # come back to this, should we add 1?

            allCs_inside=totalCpGs[allCpGsleft:allCpGsright]
            modCs_inside=totalCpGs[modCsleft:modCsright]

            freq = len(modCs_inside) / len(allCs_inside)

            # if curr methylation is > x, flag as CDR
            if freq < x:
                #print(read.query_name,"\t","methylation_drop")
                return "T", time.time() - t0

        # we hit end of read and methylation never dropped below x, call nonCDR
        return "F" , time.time() - t0

def find_cdr_candidates(bamfile, w,x,o,t ):
    '''

    '''
    nonCDRs=[]
    CDRs=[]

    parse_tag_time=0
    is_read_CDR_time=0
    empty = 0
    for read in bamfile.fetch():  # specify bamfile.fetch(contig=ref, start=pos_start, stop=pos_stop) for specific regions
        #t0=time.time()
        # later, subset reads by quality?
        # get methyl tags
        if read.has_tag('Mm') and read.has_tag('Ml'):
            # if methyl tag is there but empty, skip
            if len(read.get_tag('Mm')) == 0:
                empty+=1
                continue
            # deal with number of methyl tags being super low compared to number of Cs in read
            # output another file for reads with empty methyl tags

            # parse mm and ml tag, get list of indeces of all C bases, and indeces of all C mod bases
            allCs,modCs, parsetime = parse_tag(read)
            parse_tag_time+=parsetime
            # decide if CDR
            CDR,CDRtime=is_read_CDR(read, modCs, w, x,o,t)
            if CDR=="F":
                CDR=False
            else:
                CDR=True
            is_read_CDR_time+=CDRtime

            # append to CDR or non CDR list based on CDR bool
            if CDR:
                CDRs.append(read.query_name)
            else:
                nonCDRs.append(read.query_name)
        else:
            empty+=1
        #print('\n','processed read ',read.query_name,'in ',' %.3f' % (time.time() - t0), file=sys.stderr)
    print("w:", w, "x:", x, "o",o,"t:", t)
    print("time in parsing bam tag: ", parsetime,"seconds")
    print("time in deciding if CDR: ", is_read_CDR_time, "seconds")
    print("number of reads missing or empty Mm/Ml tag: ", empty)

    return CDRs, nonCDRs

def main():
    '''
    args.w windowsize in percent of read length
    args.j number of previous windows to consider in deciding on CDR candidacy
    args.x percent change of current window to previous j windows
    '''
    # parse command line arguments
    t0 = time.time()
    args = arg_parser()

    # open bam file
    bamfile = pysam.AlignmentFile(args.bam, "rb",check_sq=False)
    bamfile.fetch()

    # define default parameters
    w= 3000 # windowsize in bp
    x = 0.3 # threshold for methylation in sliding window
    t = 0.4 # threshold to automatically call a CDR across the whole read
    o = 1 # number bp to move sliding window over by

    # use command line parameters if they are provided
    if args.windowsize is not None:
        w=int(args.windowsize)
    if args.window_threshold is not None:
        x=float(args.window_threshold)
    if args.read_threshold is not None:
        t=float(args.read_threshold)
    if args.step_size is not None:
        o=int(args.step_size)
    # separate cdr candidates from noncdr candidates
    CDRs,nonCDRs =find_cdr_candidates(bamfile, w, x, o, t)

    # write output text files
    with open(args.out+"_CDR.txt",'w') as out:
        print(*CDRs , sep="\n",file=out)

    with open(args.out+"_nonCDR.txt",'w') as out:
        print(*nonCDRs , sep="\n",file=out)

    bamfile.close()
    print('\n', 'total time for the program %.3f' % (time.time() - t0), file=sys.stderr)

if __name__ == '__main__':
    main()
