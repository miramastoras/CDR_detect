version 1.0

workflow runCDRdetect{
    input{
        File matHORBed
        File patHORBed
        File secPhaseHifiBam
        String sampleName
        # hifi ccs bams
    }
    call getHORReadnames{
        input:
            matHORBed=matHORBed
            patHORBed=patHORBed
            secPhaseHifiBam=secPhaseHifiBam
            sampleName=sampleName
    }
    call mapPrimrose{}
    call getHORPrimrose{
        input{
            mappedPrimroseBam=mapPrimrose.mappedPrimroseBam
            sampleName=sampleName
            HORHifiBamReadnames=getHORReadnames.HORHifiBamReadnames
        }
    }
    call CDRdetect{
        input{
            HORPrimroseBam=getHORPrimrose.HORPrimroseBam
            sampleName=sampleName
        }
    }
    call getCDRfastq{
        input{
            HORPrimroseBam=getHORPrimrose.HORPrimroseBam
            CDRReadnames=CDRdetect.CDRReadnames
            sampleName=sampleName
        }
    }
    call getCDRLocations{
        input{
            hifiBam=secPhaseHifiBam
            CDRReadnames=CDRdetect.CDRReadnames
            sampleName=sampleName
            dipHORBed=getHORReadnames.dipHORBed
        }
    }
    output{
        File CDRfastq=getCDRfastq.CDRfastq
        File CDRannotations=getCDRLocations.CDRannotations
    }
}

task getHORReadnames{
    input{
        File matHORBed
        File patHORBed
        File secPhaseHifiBam
        String sampleName
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        cat ~{matHORBed} ~{patHORBed} | bedtools sort -i stdin > ~{sampleName}_AS_HOR_dip.srt.bed
        bedtools intersect -abam ~{secPhaseHifiBam} -b ~{sampleName}_AS_HOR_dip.srt.bed -wa > ~{sampleName}_hifi_diploid_HOR.bam
        samtools view ~{sampleName}_hifi_diploid_HOR.bam | cut -f1 > ~{sampleName}_hifi_diploid_HOR.readnames.txt
    >>>
    output{
        File dipHORBed = "~{sampleName}_AS_HOR_dip.srt.bed"
        File HORHifiBamReadnames = "~{sampleName}_hifi_diploid_HOR.readnames.txt"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task mapPrimrose{
    input{
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace


    >>>
    output{
       File mappedPrimroseBam = "~{sampleName}_hifi_primrose_hg38.bam"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task getHORPrimrose{
    input{
        File mappedPrimroseBam
        String SampleName
        File HORHifiBamReadnames
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 extract_reads.py -b ~{mappedPrimroseBam} -n ~{HORHifiBamReadnames} -o ~{sampleName}_hifi_primrose_hg38_HOR.bam
        samtools sort ~{sampleName}_hifi_primrose_hg38_HOR.bam > ~{sampleName}_hifi_primrose_hg38_HOR.srt.bam
        samtools index ~{sampleName}_hifi_primrose_hg38_HOR.srt.bam
    >>>
    output{
        File HORPrimroseBam = "~{sampleName}_hifi_primrose_hg38_HOR.srt.bam"
        File HORPrimroseBai = "~{sampleName}_hifi_primrose_hg38_HOR.srt.bam.bai"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task CDRdetect{
    input{
        File HORPrimroseBam
        String sampleName
        String windowSize=1500
        String windowThresh=0.5
        String readThresh=0.4
        String stepSize=1
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 CDR_detect.py -w ~{windowSize} -t ~{readThresh} -x ~{windowThresh} -i ~{stepSize} -b ~{HORPrimroseBam} -o ~{sampleName}_CDR_readnames.txt
    >>>
    output{
        File CDRReadnames="~{sampleName}_CDR_readnames.txt"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task getCDRLocations{
    input{
        File hifiBam
        File CDRReadnames
        String sampleName
        File dipHORBed
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        python3 extract_reads.py -b ~{HORHifiBam} -n ~{CDRReadnames} -o ~{sampleName}_hifi_CDRs.bam
        samtools sort ~{sampleName}_hifi_CDRs.bam > ~{sampleName}_hifi_CDRs.srt.bam
        samtools index ~{sampleName}_hifi_CDRs.srt.bam

        samtools view ~{sampleName}_hifi_CDRs.srt.bam | awk '{print $3"\t"$4"\t"$4+length($10)"\t"$1}' > ~{sampleName}_hifi_diploid_CDRreads.bed
        bedtools map -a ~{sampleName}_hifi_diploid_CDRreads.bed -b ~{dipHORBed} -c 4 -o collapse > ~{sampleName}_hifi_diploid_CDRreads_HUMAS_HMMER_annotations.txt

    >>>
    output{
        File CDRannotations = "~{sampleName}_hifi_diploid_CDRreads_HUMAS_HMMER_annotations.txt"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}

task getCDRfastq{
    input{
        File HORPrimroseBam
        File CDRReadnames
        String sampleName
    }
    command <<<
        python3 extract_reads.py -b ~{HORPrimroseBam} -n ~{CDRReadnames} -o ~{sampleName}_primrose_CDRs.bam
        bedtools bamtofastq -i ~{sampleName}_primrose_CDRs.bam -fq ~{sampleName}_CDRs.fastq
        gzip ~{sampleName}_CDRs.fastq
    >>>
    output{
        File CDRfastq="~{sampleName}_CDRs.fastq.gz"
    }
    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
