# CDR detect publication analyses

This markdown contains the analyses done for the CDR publication in benchmarking CDR_detect.py and applying it to the HPRC panel of individuals

Table of Contents:
1. [Quantifying CDRs in CHM13 and HG002-t2tX for benchmarking ](##1-quantifying-cdrs-in-chm13-and-hg002-t2tx-for-benchmarking)


## 1. Quantifying CDRs in CHM13 and HG002-t2tX for benchmarking

In order to know whether our method of detecting a CDR from a single read in a reference-free manner is performing well, we need to establish a quantitative method for defining a CDR across a pileup of CpG mod-tagged reads.

The t2t-CHM13 reference and HG002-t2tX assembly have fully resolved centromeric repeats, and their CDRs have already been presented in previous publications (Altemose et. al t2tcensat, Sup Table S11) (Altemose, et al, dimelo-seq), so we can use these as a truthset to benchmark our method before applying it to assemblies with known errors in the centromere / alpha satellites (ie HPRC).

## 2. Benchmark CDR_detect performance across CHM13 and HG002 t2t-X

## 3. Run CDR_detect across HPRC individuals

## 4. Expand CDR_detect to ONT data 
