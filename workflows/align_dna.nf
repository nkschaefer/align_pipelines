#! /usr/bin/env nextflow
import java.util.zip.GZIPInputStream
import java.nio.file.Files

process align_dna_files{
    time { 36.hour * task.attempt }
    cpus params.threads
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    val(basename),
    val(num),
    file(r1),
    file(r2),
    file(fa),
    file(idx)
    
    output:
    tuple val(lib), file("${basename}_${num}_sorted.bam"), file("${basename}_${num}_sorted.bam.bai")

    script:
    """
    minimap2 -t ${params.threads} -a -x sr -R "@RG\\tID:${basename}\\tSM:${lib}\\tPL:Illumina" ${idx} ${r1} ${r2} \
| samtools sort -n - | samtools fixmate -m - - | samtools sort -o ${basename}_${num}_sorted.bam
    samtools index ${basename}_${num}_sorted.bam    
    """    
   
}

process cat_dna_bams{
    time { 12.hour * task.attempt }
    cpus params.threads
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    file(bams),
    file(bais)
    
    output:
    tuple val(lib), file("${lib}_merged.bam"), file("${lib}_merged.bam.bai")

    script:
    """
    samtools merge -@ ${params.threads} ${lib}_merged.bam ${bams}
    samtools index ${lib}_merged.bam
    """
}

process dna_mkdup{
    time '24h'
    cpus params.threads
    
    input:
    tuple val(lib),
    file(bam),
    file(bai)
    
    publishDir "${params.output_directory}", mode: 'copy'

    output:
    tuple val(lib),
    file("${lib}.bam"),
    file("${lib}.bam.bai")

    script:
    """
    samtools markdup -@ ${params.threads} -m t ${bam} ${lib}.bam
    samtools index ${lib}.bam    
    """
}

process split_fai_regions{
    time '10m'
    cpus 1
    
    input:
    path(fa)
    
    output:
    tuple path(fa), path("*.bed")

    script:
    """
    samtools faidx ${fa}
    ${baseDir}/segment_genome.py -n 30 -d . -f ${fa}.fai 
    """
}

process varcall{
    cpus 1
    time { 36.hour * task.attempt }
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
    tuple path(bams), path(bais), path(fa), path(bed)

    output:
    path("*.vcf")

    script:
    """
    for bam in ${bams}; do
        echo \$bam >> bamlist.txt
    done
    freebayes -f ${fa} -L bamlist.txt -t ${bed} --use-best-n-alleles 4 --min-alternate-count 4 > ${bed}.vcf
    """
}

process cat_vcf{
    cpus 1 
    time 8.hour

    input:
    path(vcfs)
    
    publishDir "${params.output_directory}", mode: 'copy'

    output:
    tuple path("*.vcf.gz"), path("*.vcf.gz.tbi")

    script:
    """
    echo -e \"${vcfs}\" | tr ' ' '\\n' | sort -V | while read fn; do
        if [ ! -e vars.all.vcf ]; then
            cat \$fn > vars.all.vcf
        else
            if [ \$( cat \$fn | grep -v \"#\" | wc -l ) -gt 0 ]; then
                cat \$fn | grep -v \"#\" >> vars.all.vcf
            fi
        fi
    done
    bgzip vars.all.vcf
    tabix -p vcf vars.all.vcf.gz
    """
}

process mm2_idx{
    time 2.hour
    cpus 1

    input:
    path(fa)

    output:
    tuple path(fa), path("*.mm2")

    script:
    """
    minimap2 -d ${fa}.mm2 ${fa}
    """
}

process filt_vcf{
    cpus 1
    time { 36.hour * task.attempt }
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple path(vcf), path(tbi)

    publishDir "${params.output_directory}", mode: 'copy'
    
    output:
    tuple path("*.filt.vcf.gz"), path("*.filt.vcf.gz.tbi"), path("*.filt.hist"), path("*cov.pdf"), path("*hetmiss.pdf") 

    script:
    def vcf_filt = (vcf.toString() =~ /(.+)\.vcf(.gz)?/)[0][1] + ".filt.vcf.gz"
    def vcf_hist = (vcf.toString() =~ /(.+)\.vcf(.gz)?/)[0][1] + ".filt.hist"
    """
    ${baseDir}/vcf_depth_filter -v ${vcf} -o ${vcf_filt}
    tabix -p vcf ${vcf_filt}
    ${baseDir}/plot_filt_stats.R ${vcf_hist}
    """


}

process split_reads{
    input:
    tuple val(libname), val(fnbase), path(r1), path(r2)

    output:
    tuple val(libname),
          val(fnbase),
          path("*_R1_001.*.fastq.gz"),
          path("*_R2_001.*.fastq.gz")

    script:
    """
    ${baseDir}/split_read_files -1 ${r1} -2 ${r2} -o . -n ${params.num_chunks}
    """
}

workflow call_variants{
    take:
    libs
    
    main:
    
    if (!params.dna_ref){
        error("FASTA reference genome is required")
    }
    if (!params.dna_dir){
        error("Genomic DNA reads directory is required")
    }
    
    def fa_genome = Channel.fromPath(params.dna_ref) | mm2_idx
    
    libs.collect().map{ it as Set }.set{ lib_set }
    
    // Get FASTQ file names
    def dna_pairs = Channel.fromPath("${params.dna_dir}/*_R1*.fastq.gz").map{ fn -> 
        def r1 = fn.toString().trim()
        def match = (r1 =~ /(.*)\/(.*)\_R1(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
        def libname = match[2].replaceFirst(/_S(\d+)_L(\d+)/, '')
        def dirn = ""
        if (match[1] != null && match[1] != ""){
            dirn += match[1] + '/'
        }
        def end = ""
        if (match[3] != null){
            end = match[3]
        }
        def gz = ""
        if (match[5] != null){
            gz = match[5]
        }
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def fnbase = match[2]
        return [ libname, fnbase, file(r1), file(r2)]
    }.filter{
        item -> 
        def (key, _) = item
        lib_set.contains(item)
    }

    dna_bams = split_reads(dna_pairs).flatMap{ ln, fsub, files1, files2 -> 
        files1.indices.collect{ i -> [ln, fsub, i, files1[i], files2[i]] }
    }.combine(fa_genome) | align_dna_files

    bams_mkdup = dna_bams.groupTuple() | cat_dna_bams | dna_mkdup
    bamslist = bams_mkdup.map{ tup -> 
        tup[1]
    }
    baislist = bams_mkdup.map{ tup -> 
        tup[2]
    }
    
    // Extract BED files for breaking up genomic seqs (to run var calling in parallel)
    beds = split_fai_regions(Channel.fromPath(params.dna_ref)).flatMap{ fa, bedlist -> 
        bedlist.collect{ bed -> [ fa, bed ] }
    } 
    vcfs = varcall(bamslist.collect().toList().combine(baislist.collect().toList()).combine(beds))
    vcf_concat = cat_vcf(vcfs.collect())
    filtfiles = filt_vcf(vcf_concat)
}
