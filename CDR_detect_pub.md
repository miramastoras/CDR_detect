# CDR detect publication analyses

This markdown contains the analyses done for the CDR publication in benchmarking CDR_detect.py and applying it to the HPRC panel of individuals

Table of Contents:
1. [Quantifying CDRs in CHM13 and HG002-t2tX for benchmarking ](##1-quantifying-cdrs-in-chm13-and-hg002-t2tx-for-benchmarking)


## 1. Quantifying CDRs in CHM13 and HG002-t2tX for benchmarking

In order to know whether our method of detecting a CDR from a single read in a reference-free manner is performing well, we need to establish a quantitative method for defining a CDR across a pileup of CpG mod-tagged reads.

The t2t-CHM13 reference and HG002-t2tX assembly have fully resolved centromeric repeats, and their CDRs have already been presented in previous publications (Altemose et. al t2tcensat, Sup Table S11) (Altemose, et al, dimelo-seq), so we can use these as a truthset to benchmark our method before applying it to assemblies with known errors in the centromere / alpha satellites (i.e. HPRC assemblies).

### 1.1 Comparing two ways to smooth methylation calls into windows in an aligned bamfile

#### 1.1A Preparing data
Starting with an aligned bamfile containing [CpG mod samtags](https://samtools.github.io/hts-specs/SAMtags.pdf), we use [pb-cpg tools](https://github.com/PacificBiosciences/pb-CpG-tools) to get 5mC modification probabilities at every CpG site in the genome.

```
# model mode used instead of count mode because it was recommended by authors
conda activate cpg
python /Users/miramastoras/Desktop/Miga_lab/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b /Users/miramastoras/Desktop/IGV_files/S3CXH1L.hg002_t2tX.srt.bam -f /Users/miramastoras/Desktop/IGV_files/HG002_t2t_chrX.fa -o /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg -p model -d /Users/miramastoras/Desktop/Miga_lab/pb-CpG-tools/pileup_calling_model/ -m reference
```
Next, we make windows of size 500 and 1000 bp across the regions represented in `S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed`
```
python3 /Users/miramastoras/Desktop/Miga_lab/CDR_detect/scripts/make_bed_windows.py -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed -s 500

python3 /Users/miramastoras/Desktop/Miga_lab/CDR_detect/scripts/make_bed_windows.py -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed -s 1000
```

##### Method 1: calling each CpG site as 100% or 0% first, then averaging the 100 and 0s in each window

Because pb-cpg tools outputs a probability of methylation at every site (# reads supporting methylation / total number of reads), in this method we convert all probabilities >50% to 100%, and all probabilities <50% to 0% first

```
awk -v OFS='\t' '{if ($4 > 50) {print $1,$2,$3,100} else {print $1,$2,$3,0}}' /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed > /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.0or1bed
```
Next, use bedtools map to get smoothed methylation calls inside of the 500bp and 1000 bp windows
```
bedtools map -a S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed_windows1000.bed -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.0or1bed -c 4 -o mean > /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.0or1.smoothed1000.bedgraph

bedtools map -a S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed_windows500.bed -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.0or1bed -c 4 -o mean > /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.0or1.smoothed500.bedgraph
```
##### Method 2: Just taking the average of the CpG probability produced by pb-cpg tools in each window.

Directly smoothing the methylation calls output by pb-cpg tools
```
bedtools map -a S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed_windows1000.bed -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed -c 4 -o mean > /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.smoothed1000.bedgraph

bedtools map -a S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed_windows500.bed -b /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.bed -c 4 -o mean > /Users/miramastoras/Desktop/Miga_lab/S3CXHL/S3CXH1L.hg002_t2tX.pb_cpg.combined.reference.smoothed500.bedgraph
```

##### Comparing Method 1 and Method 2:

![method1vs2](pics/CDR0or1vsavg.png)

Although the dips are sharper with method 1, we will still have the same regions being called when we pull windows <50% methylated, so it won't make a difference which method we use. I am more inclined toward method 2, because it takes into account the uncertainty in the methylation calls at each site.  Therefore, I will move ahead with method 1.

### 1.2 Choosing a strict and a lenient quantitative definition for a CDR

CDRs are regions of hypomethylation amongst regions of hypermethylation in highly repetitive alpha arrays. They have been associated with the sites of cenpA binding and kinetochore attachment. However, it is not known what the minimum size of a CDR is for cenpa to bind it. In the HG002 chrX array for example, when we look on the scale of the whole array, the 5 previously identified  CDR dips are clear to see.

![whole x array](pics/whole_array_x.png)

But, smaller CDR dips are also clear to see when we look on a closer scale

![tinyCDR](pics/tiny_CDR_in.png)

If we are trying to associate regions which bind cenpA using CUT&run with regions containing CDRs, we need to know what defines a CDR, we need to know the size of the dip we are looking for. Excluding these smaller dips may be excluding real regions which bind cenpa, but it also  may confuse our ability to associate repeats with cenpA, if they are actually below the size required to bind it (which is unknown). In my mind this provides a reason to only search for larger CDRs, because we need a different approach to answer the question about minimum CDR size.

In looking to create a quantitative definition of a CDR, I tested several options and visually inspected them against the smoothed methylation windows created in 1.1 and the S3CXHL bam file.

![quant options](pics/quantCDR_options.png)

They vary between using 500 bp windows and 1000 bp windows, merging adjacent windows or windows separated by 500 or 1000 bp, and setting a minimum window size after merging (no minimum, 1000 and 1500)

For my large CDR definition, I've chosen 1000 bp windows, merging adjacent windows <50% methylated and a minimum CDR size of 1500 bp

for my small CDR definition, I've chosen 1000 bp windows, merging adjacent windows <50% methylated and a minimum CDR size of 1500 bp

#### 1.2A Looking at dips in CHM13

fix chrom names in chm13 primrose bamfile.
Key located here: https://www.ncbi.nlm.nih.gov/assembly/GCA_009914755.3

```
grep -v "^#" GCA_009914755.3_T2T-CHM13v1.1_assembly_report.txt | cut -f1,5 | awk '{print "sed -e " "\s/chr"$1"/"$2"/"}'

cat t2t_cenAnnotation.v2.021921.bed | sed -e 's/chr1/CP068277.2/' | sed -e 's/chr2/CP068276.2/' | sed -e 's/chr3/CP068275.2/' | sed -e 's/chr4/CP068274.2/' | sed -e 's/chr5/CP068273.2/' | sed -e 's/chr6/CP068272.2/' | sed -e 's/chr7/CP068271.2/' | sed -e 's/chr8/CP068270.2/' | sed -e 's/chr9/CP068269.2/' | sed -e 's/chr10/CP068268.2/' | sed -e 's/chr11/CP068267.2/' | sed -e 's/chr12/CP068266.2/' | sed -e 's/chr13/CP068265.2/' | sed -e 's/chr14/CP068264.2/' | sed -e 's/chr15/CP068263.2/' | sed -e 's/chr16/CP068262.2/' | sed -e 's/chr17/CP068261.2/' | sed -e 's/chr18/CP068260.2/' | sed -e 's/chr19/CP068259.2/' | sed -e 's/chr20/CP068258.2/' | sed -e 's/chr21/CP068257.2/' | sed -e 's/chr22/CP068256.2/' | sed -e 's/chrX/CP068255.2/' | sed -e 's/MT/CP068254.1/' > t2t_cenAnnotation.v2.021921.genbank_chrnames.bed
```


Pull out alpha arrays from [chm13-chm13v1.1 primrose](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/a211382e-938f-4aa2-8424-1db8fc75a0cd--CHM13-HIFI-METHYLATION-CHM13v1.1/) (on wakanda)

```
cd /data/mira/CDR_detect/data
grep "hor" t2t_cenAnnotation.v2.021921.bed | grep -v "dhor" > t2t_cenAnnotation.v2.021921.hor.bed
bedtools intersect -abam /scratch/mira/chm13.CHM13_v1.1.bam -b /data/mira/CDR_detect/data/t2t_cenAnnotation.v2.021921.genbank_chrnames.bed -wa > /scratch/mira/chm13.CHM13_v1.1.hor.bam
```

Download and run pb-cpg tools on it to get methylation calls per site
```
conda activate cpg
python /Users/miramastoras/Desktop/Miga_lab/pb-CpG-tools/aligned_bam_to_cpg_scores.py -b /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/chm13.CHM13_v1.1.hor.bam -f /Users/miramastoras/Desktop/IGV_files/HG002_t2t_chrX.fa -o /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/CHM13v1.1.pb_cpg -p model -d /Users/miramastoras/Desktop/Miga_lab/pb-CpG-tools/pileup_calling_model/ -m reference
```
Smooth results into 1000 bp windows
```
# make windows
python3 /Users/miramastoras/Desktop/Miga_lab/CDR_detect/scripts/make_bed_windows.py -b /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/CHM13v1.1.pb_cpg.combined.reference.bed -s 1000
# smooth cpg in windows
bedtools map -a /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/CHM13v1.1.pb_cpg.combined.reference.bed_windows1000.bed -b /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/CHM13v1.1.pb_cpg.combined.reference.bed -c 4 -o mean > /Users/miramastoras/Desktop/Miga_lab/CHM13_hor/CHM13v1.1.pb_cpg.combined.reference.smoothed1000.bedgraph
```

## 2. Benchmark CDR_detect performance across CHM13 and HG002 t2t-X

## 3. Run CDR_detect across HPRC individuals hifi data

## 4. Expand CDR_detect to ONT data
