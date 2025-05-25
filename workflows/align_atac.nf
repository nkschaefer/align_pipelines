#! /usr/bin/env nextflow
import java.util.zip.GZIPInputStream
import java.nio.file.Files

def atac_map = [:]
if (!params.demux_species){
    if (params.atac_map){
        new File(params.atac_map).eachLine { line ->
            def (atac, lib) = line.split('\t')
            atac_map[atac] = lib
        }
    } else{
        if (!params.libs){
            error("libs file required.")
        }
        new File(params.libs).eachLine { line ->
            def linetrim = line.trim()
            atac_map[linetrim] = linetrim
        }
    }
}

def atac_ref_map = [:]
if (params.atac_ref_species){
    new File(params.atac_ref_species).eachLine { line ->
        def (species, ref) = line.split('\t')
        atac_ref_map[species] = ref
    }
}

process preproc_atac_files{
    time '24h'
    
    input:
    tuple val(lib),
    val(basename),
    file(r1),
    file(r2),
    file(r3),
    file(idx),
    file(wl)
    
    output:
    tuple val(lib),
    val(basename),
    file("preproc/*_R1*.fastq.gz"),
    file("preproc/*_R2*.fastq.gz"),
    file(idx)

    script:
    """
    mkdir preproc
    ${baseDir}/atac_fq_preprocess -1 ${r1} -2 ${r2} -3 ${r3} -o preproc -w ${wl}
    """
}
process preproc_atac_files_multiome{
    time '24h'
    
    input:
    tuple val(lib),
    val(basename),
    file(r1),
    file(r2),
    file(r3),
    file(idx),
    file(wl_rna),
    file(wl_atac)
    
    output:
    tuple val(lib),
    val(basename),
    file("preproc/*_R1*.fastq.gz"),
    file("preproc/*_R2*.fastq.gz"),
    file(idx)

    script:
    """
    mkdir preproc
    ${baseDir}/atac_fq_preprocess -1 ${r1} -2 ${r2} -3 ${r3} -o preproc -w ${wl_rna} -W ${wl_atac}
    """
}

process align_atac_files{
    time { 36.hour * task.attempt }
    cpus params.threads
    memory params.memgb + ' GB'
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    val(basename),
    val(num),
    file(idx),
    file(r1),
    file(r2)
    
    output:
    tuple val(lib), file("${basename}_${num}_sorted.bam"), file("${basename}_${num}_sorted.bam.bai")

    script:
    """
    minimap2 -t ${params.threads} -y -a -x sr -R "@RG\\tID:${basename}\\tSM:${lib}\\tPL:Illumina" ${idx} ${r1} ${r2} \
| samtools sort -n - | samtools fixmate -m - - | samtools sort -o ${basename}_${num}_sorted.bam
    samtools index ${basename}_${num}_sorted.bam    
    """    
   
}

process cat_atac_bams{
    time { 12.hour * task.attempt }
    cpus params.threads
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    file(bams),
    file(bais)
    
    output:
    tuple val(lib), file("atac_merged.bam"), file("atac_merged.bam.bai")

    script:
    """
    samtools merge -@ ${params.threads} atac_merged.bam ${bams}
    samtools index atac_merged.bam
    """
}

process atac_mkdup{
    time '24h'
    cpus params.threads
    
    input:
    tuple val(lib),
    file(bam),
    file(bai)
    
    output:
    tuple val(lib),
    file("atac.bam"),
    file("atac.bam.bai")

    script:
    """
    samtools markdup -@ ${params.threads} -m t --barcode-tag CB ${bam} atac.bam
    samtools index atac.bam    
    """
}

process atac_namesort{
    time '12h'
    cpus params.threads

    input:
    tuple val(lib),
    file(bam),
    file(bai)
    
    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: {path -> path }
    
    output:
    tuple val(lib),
    file("atac.bam"),
    file("atac.bam.bai"),
    file("atac_namesort.bam")

    script:
    """
    samtools sort -@ ${params.threads} -n -o atac_namesort.bam ${bam}
    """
}

process atac_fragments{
    time { 12.hour * task.attempt }
    cpus params.threads
    
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3
    
    input:
    tuple val(lib),
    file(bam),
    file(bai),
    file(bam_namesort)
    
    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: {path -> path }
    
    output:
    tuple val(lib),
    file("atac_fragments.tsv.gz"),
    file("atac_fragments.tsv.gz.tbi")

    script:
    """
    sinto fragments -p ${params.threads} --collapse_within -m 30 -t CB \
--use_chrom "." -b ${bam} -f atac_fragments.tsv
    cat atac_fragments.tsv | sort -k1,1V -k2,2n -k3,3n | bgzip > atac_fragments.tsv.gz
    tabix -s 1 -b 2 -e 3 atac_fragments.tsv.gz 
    """
}

/**
 * Checks the first line of a gzipped FASTQ file to see if a cell barcode tag
 * has been appended.
 */
def peek_check_bc(filePath){
    def line = ""
    Files.newInputStream(filePath).withCloseable { stream ->
        new GZIPInputStream(stream).withReader { reader ->
            line = new BufferedReader(reader).readLine()        
        }
    }
    lnsplit = line.split(' ')
    if ( lnsplit[-1] ==~ /^CB:Z:[ACGT]+$/){
        return true
    }
    else{
        return false
    }
}

workflow align_atac_demux_species{
    main:
    
    if (!params.atac_ref_species){
        error("ATAC reference is required")
    }
    if (!params.demux_species){
        error("demux_species output dir is required")
    }
    
    def atac_triples = Channel.fromPath("${params.demux_species}/*/*/ATAC*_R3*.fastq.gz").map{ fn -> 
        def r3 = fn.toString().trim()
        def libn = r3.split('/')[-3]
        def species = r3.split('/')[-2]
        
        def match = (r3 =~ /(.*)\/(.*)\_R3(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
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
        def r1 = dirn + match[2] + "_R1" + end + '.' + match[4] + gz
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def fnbase = match[2]
        if (! atac_ref_map[species]){
            error("Species " + species + " does not have an ATAC reference specified")
        }
        return [ libn + "/" + species, fnbase, file(r1), file(r2), file(r3), file(atac_ref_map[species]) ]
    }
    
    if (params.multiome){
        if (!params.rna_whitelist || !params.atac_whitelist){
            error("Both rna_whitelist and atac_whitelist are required if multiome data")
        }
        def wl_rna = Channel.fromPath(params.rna_whitelist)
        def wl_atac = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files_multiome(atac_triples.combine(wl_rna).combine(wl_atac))
    }
    else{
        if (!params.atac_whitelist){
            error("ATAC whitelist required")
        }
        def wl = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files(atac_triples.combine(wl))
    } 

    def atac_pairs = Channel.fromFilePairs("${params.demux_species}/*/*/ATAC*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        def libn = reads[0].toString().split('/')[-3]
        def species = reads[0].toString().split('/')[-2]
        if (! atac_ref_map[species]){
            error("Species " + species + " does not have an ATAC reference specified")
        }
        [ libn + "/" + species, id, reads[0], reads[1], file(atac_ref_map[species]) ]
    }.filter{ lib, fnbase, r1, r2, idx ->
        return peek_check_bc(r1) && peek_check_bc(r2)
    }
    
    atac_bams = atac_preproc1.concat(atac_pairs) | align_atac_files
    
    atac_bams.groupTuple() | cat_atac_bams | atac_mkdup | atac_namesort | atac_fragments
}

process split_reads_atac{
    input:
    tuple val(libname), val(fnbase), path(r1), path(r2), path(idx)

    output:
    tuple val(libname),
    val(fnbase),
    path(idx),
    path("*_R1_001.*.fastq.gz"),
    path("*_R2_001.*.fastq.gz")

    script:
    """
    ${baseDir}/split_read_files -1 ${r1} -2 ${r2} -o . -n ${params.num_chunks}
    """
}

workflow align_atac{
    take:
    libs
    
    main:
    
    if (!params.atac_ref){
        error("ATAC reference is required")
    }
    if (!params.atac_dir){
        error("ATAC reads directory is required")
    }
    if (params.num_chunks < 2){
        error("num_chunks must be at least 2")
    }

    def idx_atac = Channel.fromPath(params.atac_ref)
    
    // Get non-preprocessed files
    def atac_triples = libs.cross(Channel.fromPath("${params.atac_dir}/*_R3*.fastq.gz").map{ fn ->
        def r3 = fn.toString().trim()
        def match = (r3 =~ /(.*)\/(.*)\_R3(_\d+)?\.(fastq|fq)(\.gz)?/)[0]
        def libname_atac = match[2].replaceFirst(/_S(\d+)_L(\d+)/, '')
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
        def r1 = dirn + match[2] + "_R1" + end + '.' + match[4] + gz
        def r2 = dirn + match[2] + "_R2" + end + '.' + match[4] + gz
        def fnbase = match[2]
        return [ atac_map[libname_atac], fnbase, file(r1), file(r2), file(r3)]
    }).map{ lib, tup -> tup }
    
    if (params.multiome){
        if (!params.atac_whitelist || !params.rna_whitelist){
            error("rna_whitelist and atac_whitelist are required for multiome data.")
        }
        def wl_rna = Channel.fromPath(params.rna_whitelist)
        def wl_atac = Channel.fromPath(params.atac_whitelist)    
        atac_preproc1 = preproc_atac_files_multiome(atac_triples.combine(idx_atac).combine(wl_rna).combine(wl_atac))
    }
    else{
        if (!params.atac_whitelist){
            error("atac_whitelist is required")
        }
        def wl = Channel.fromPath(params.atac_whitelist)
        atac_preproc1 = preproc_atac_files(atac_triples.combine(idx_atac).combine(wl))
    } 
    
    // Get pre-processed files
    def atac_pairs = libs.cross(
        Channel.fromFilePairs("${params.atac_dir}/*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        [ atac_map[id.replaceFirst(/_S(\d+)_L(\d+)/, '')], id, reads[0], reads[1]]
    }).map{ lib, tup -> tup }.filter{ lib, fnbase, r1, r2 ->
        return peek_check_bc(r1) && peek_check_bc(r2)
    }
    atac_preproc2 = atac_pairs.combine(idx_atac)
    
    atac_bams = split_reads_atac(atac_preproc1.concat(atac_preproc2)).flatMap{ 
        ln, fsub, idx, files1, files2 ->
            files2.indices.collect{ i -> [ln, fsub, i, idx, files1[i], files2[i]] }
    } | align_atac_files

    atac_bams.groupTuple() | cat_atac_bams | atac_mkdup | atac_namesort | atac_fragments

}
