'''
Purpose: Evaluate performance of CDR detect in classifying CDR vs nonCDR reads
Author: Mira Mastoras, mmastora@ucsc.edu
Usage:
python3 evaluate_CDR_detect.py -c /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.CDR_500.quant.readnames.txt \
-n /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.nonCDR_500.quant.readnames.txt \
-b /Users/miramastoras/Desktop/IGV_files/S3CXH1L.hg002_t2tX.srt.bam \
-o /Users/miramastoras/Desktop/test_CDR500_CDR_detect_coords.txt
'''

import argparse
import pysam
import numpy as np
import time,sys
import CDR_detect

def arg_parser():
    '''
    Parses command line arguments with argparse
    '''
    parser = argparse.ArgumentParser(
        prog='evaluate_CDR_detect.py',
        description="""Evaluate performance of CDR detect in classifying CDR vs nonCDR reads """)

    parser.add_argument("-b", "--bam",
                        required=True,
                        metavar="input.bam",
                        help="The aligned BAM file to clasify.")
    parser.add_argument("-o", "--out",
                        required=True,
                        metavar="output tsv file prefix",
                        help="prefix for output list of parameters and their performance")
    parser.add_argument("-c", "--truthCDR",
                        required=True,
                        help="text file listing read names of truth CDRs")
    parser.add_argument("-n", "--truthNonCDR",
                        required=True,
                        help="text file listing read names of truth nonCDRs")

    return parser.parse_args()

def calc_performance(CDRs,nonCDRs,tCDR,tnonCDR,args):
    '''
    evaluates precision, recall, f1 of cdr classifications
    '''
    print(len(CDRs), len(nonCDRs), len(tCDR), len(tnonCDR))
    # calc true positives (called CDR, actually CDR)
    tp=len(np.intersect1d(CDRs,tCDR)) # yields the elements in CDRs that are NOT in tCDR
    # calc false positives (called CDR, actually nonCDR)
    fp=len(np.setdiff1d(CDRs,tCDR))
    # calc true negatives (called nonCDR, actually non CDR)
    tn=len(np.intersect1d(nonCDRs,tnonCDR))
    # calc false negatives (called nonCDR, actually CDR)
    fn=len(np.setdiff1d(nonCDRs,tnonCDR))

    # calc p,r,f1
    p= tp / (tp+fp)
    r= tp / (tp+fn)
    f1 = 2 * ((p*r) / (p + r))

    FPreads= np.setdiff1d(CDRs,tCDR)
    FNreads= np.setdiff1d(nonCDRs,tnonCDR)
    TPreads= np.intersect1d(CDRs,tCDR)
    TNreads= np.intersect1d(nonCDRs,tnonCDR)

    with open(args.out + "_FP.txt","w") as out:
        print(*FPreads, sep="\n",file=out)
    with open(args.out + "_FN.txt","w") as out:
        print(*FNreads, sep="\n", file=out)

    with open(args.out + "_TP.txt","w") as out:
        print(*TPreads, sep="\n",file=out)
    with open(args.out + "_TN.txt","w") as out:
        print(*TNreads, sep="\n", file=out)

    return p,r,f1,tp,fp,tn,fn,len(tCDR),len(tnonCDR)

def main():
    '''

    '''
    # parse command line arguments
    t0 = time.time()
    args = arg_parser()

    # read in list of truth CDRs and truth non CDRs
    truthNonCDR=open(args.truthNonCDR)
    truthCDR=open(args.truthCDR)

    tCDR=[]
    tnonCDR=[]

    # get lists of truth CDRs and truth nonCDRs
    for line in truthCDR:
        tCDR.append(line.split())
    for line in truthNonCDR:
        tnonCDR.append(line.split())

    truthNonCDR.close()
    truthCDR.close()

    # open bamfile
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    bamfile.fetch()

    # loop through all combinations of parameters,
    for w in [500,1000,2000,3000,4000,5000,6000,10000]:
         for x in [0.25,0.3,0.35,0.4,0.5]:
             for n in [500,1000,2000,3000]:
                 CDR_results = CDR_detect.find_cdr_candidates(bamfile, w, x, n)
                 CDRs = []
                 for c in CDR_results:
                     CDRs.append(c[0])

                 # get nonCDR reads from those not labelled CDR
                 allReads = tCDR + tnonCDR
                 nonCDRs = np.setdiff1d(allReads, CDRs)
                 p, r, f1, tp, fp, tn, fn, totalPosCDR, totalNegnonCDR = calc_performance(CDRs,nonCDRs,tCDR,tnonCDR,args)
                 with open(args.out, "a") as outFile:
                     print(w, x, n, p, r, f1, tp, fp, tn, fn, totalPosCDR, totalNegnonCDR, sep="\t", file=outFile)

    # define parameters (later switch to command line)
    #w = 1000  # windowsize in bp
   #  x = 0.3  # threshold for methylation in window
   #  n = 3000  #
   #  CDR_results = CDR_detect.find_cdr_candidates(bamfile,w,x,n)
   #  CDRs=[]
   #  for c in CDR_results:
   #      CDRs.append(c[0])
   #
   #  # get nonCDR reads from those not labelled CDR
   #  allReads=tCDR + tnonCDR
   #  nonCDRs = np.setdiff1d(allReads, CDRs)
   #  # yields the elements in `list_2` that are NOT in `list_1`
   #  #with open(args.out + "_predictedCDR.txt", "w") as out:
   #      #print(*CDRs, sep="\n", file=out)
   # # with open(args.out + "_predictednonCDR.txt", "w") as out:
   #      #print(*nonCDRs, sep="\n", file=out)
   #  # evaluate performance
   #  p,r,f1,tp,fp,tn,fn,totalPosCDR,totalNegnonCDR = calc_performance(CDRs,nonCDRs,tCDR,tnonCDR,args)
   #          # write to output file
   #  with open(args.out, "a") as outFile:
   #      print(w,x, n,p,r,f1,tp,fp,tn,fn,totalPosCDR,totalNegnonCDR,sep="\t",file=outFile)


    bamfile.close()

    print('\n', 'total time for the program %.3f' % (time.time() - t0), file=sys.stderr)
if __name__ == '__main__':
    main()


#