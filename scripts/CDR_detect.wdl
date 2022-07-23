version 1.0

workflow runCDRdetect{
    input{
        File matHORBed
        File patHORBed
        File secPhaseHifiBam
        String sampleName
        File mappedPrimroseBam
    }

    call getHORReads{
        input:
            matHORBed=matHORBed,
            patHORBed=patHORBed,
            secPhaseHifiBam=secPhaseHifiBam,
            sampleName=sampleName,
            mappedPrimroseBam=mappedPrimroseBam
    }

    call CDRdetect{
        input:
            primroseHORBam=getHORReads.primroseHORBam,
            primroseHORBai=getHORReads.primroseHORBai,
            sampleName=sampleName
    }

    call formatResults{
        input:
            primroseHORBam=getHORReads.primroseHORBam,
            primroseHORBai=getHORReads.primroseHORBai,
            CDRReadnames=CDRdetect.CDRReadnames,
            sampleName=sampleName,
            hifiHORBam=getHORReads.hifiHORBam,
            hifiHORBai=getHORReads.hifiHORBai,
            dipHORBed=getHORReads.dipHORBed
    }

    output{
        File CDRfastq=formatResults.CDRfastq
        File CDRannotations=formatResults.CDRannotations
    }
}

task getHORReads{
    input{
        File matHORBed
        File patHORBed
        File secPhaseHifiBam
        String sampleName
        File mappedPrimroseBam

        String dockerImage = "miramastoras/cdr_detect:latest"
        Int memSizeGB = 24
        Int threads = 2
        Int diskSizeGB = 128
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # get list of readnames mapping to HOR bedfiles in the diploid assemblies
        cat ~{matHORBed} ~{patHORBed} | bedtools sort -i stdin > ~{sampleName}_AS_HOR_dip.srt.bed
        bedtools intersect -abam ~{secPhaseHifiBam} -b ~{sampleName}_AS_HOR_dip.srt.bed -wa > ~{sampleName}_hifi_diploid_HOR.bam
        samtools sort --threads ~{threads} ~{sampleName}_hifi_diploid_HOR.bam > ~{sampleName}_hifi_diploid_HOR.srt.bam
        samtools index -@ ~{threads} ~{sampleName}_hifi_diploid_HOR.srt.bam
        samtools view --threads ~{threads} ~{sampleName}_hifi_diploid_HOR.srt.bam | cut -f1 > ~{sampleName}_HOR.readnames.txt

        # pull HOR readnames from primrose data
        python3 /opt/CDR_detect/scripts/extract_reads.py -b ~{mappedPrimroseBam} -n ~{sampleName}_HOR.readnames.txt -o ~{sampleName}_hifi_primrose_HOR.bam
        samtools sort --threads ~{threads} ~{sampleName}_hifi_primrose_HOR.bam > ~{sampleName}_hifi_primrose_HOR.srt.bam
        samtools index -@ ~{threads} ~{sampleName}_hifi_primrose_HOR.srt.bam
    >>>
    output{
        File dipHORBed = "~{sampleName}_AS_HOR_dip.srt.bed"
        File hifiHORBam = "~{sampleName}_hifi_diploid_HOR.srt.bam"
        File hifiHORBai = "~{sampleName}_hifi_diploid_HOR.srt.bam.bai"
        File primroseHORBam = "~{sampleName}_hifi_primrose_HOR.srt.bam"
        File primroseHORBai = "~{sampleName}_hifi_primrose_HOR.srt.bam.bai"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task CDRdetect{
    input{
        File primroseHORBam
        File primroseHORBai
        String sampleName

        String windowSize=3000
        String windowThresh=0.3
        String readThresh=0.4
        String stepSize=1

        String dockerImage = "miramastoras/cdr_detect:latest"
        Int memSizeGB = 24
        Int threads = 2
        Int diskSizeGB = 128
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        CDR_detect.py -w ~{windowSize} -t ~{readThresh} -x ~{windowThresh} -i ~{stepSize} -b ~{primroseHORBam} -o ~{sampleName}
    >>>
    output{
        File CDRReadnames="~{sampleName}_CDR.txt"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

task formatResults{
    input{
        File hifiHORBam
        File hifiHORBai
        File CDRReadnames
        File dipHORBed
        File primroseHORBam
        File primroseHORBai
        String sampleName

        String dockerImage = "miramastoras/cdr_detect:latest"
        Int memSizeGB = 24
        Int threads = 2
        Int diskSizeGB = 128
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # get locations & annotations for each read
        python3 /opt/CDR_detect/scripts/extract_reads.py -b ~{hifiHORBam} -n ~{CDRReadnames} -o ~{sampleName}_hifi_CDRs.bam
        samtools sort --threads ~{threads} ~{sampleName}_hifi_CDRs.bam > ~{sampleName}_hifi_CDRs.srt.bam
        samtools index -@ ~{threads} ~{sampleName}_hifi_CDRs.srt.bam

        samtools view --threads ~{threads} ~{sampleName}_hifi_CDRs.srt.bam | awk '{print $3"\t"$4"\t"$4+length($10)"\t"$1}' > ~{sampleName}_hifi_diploid_CDRreads.bed
        bedtools map -a ~{sampleName}_hifi_diploid_CDRreads.bed -b ~{dipHORBed} -c 4 -o collapse > ~{sampleName}_hifi_diploid_CDRreads_HUMAS_HMMER_annotations.txt

        # get fastq sequence using primrose data
        python3 /opt/CDR_detect/scripts/extract_reads.py -b ~{primroseHORBam} -n ~{CDRReadnames} -o ~{sampleName}_primrose_CDRs.bam
        bedtools bamtofastq -i ~{sampleName}_primrose_CDRs.bam -fq ~{sampleName}_CDRs.fastq
        gzip ~{sampleName}_CDRs.fastq

    >>>
    output{
        File CDRannotations = "~{sampleName}_hifi_diploid_CDRreads_HUMAS_HMMER_annotations.txt"
        File CDRfastq="~{sampleName}_CDRs.fastq.gz"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}
