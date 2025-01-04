#! /usr/bin/env nextflow
nextflow.enable.dsl=2

if (!params.output_directory){
    error("Output directory is required.")
}
if (!params.fasta){
    error("Genome FASTA is required.")
}
if (!params.genome_base){
    error("genome_base (Base name for reference files) required.")
}

process make_star_index{
    time '24h'
    memory params.memgb + ' GB'
    
    input:
    tuple file(fasta),
    file(gtf)
    
    publishDir "${params.output_directory}", mode: 'copy', saveAs: {path -> path }
    
    output:
    file("${params.genome_base}/*")
    
    script:
    def mem_int = params.memgb.toInteger()
    def mem_int_bytes = (mem_int - 1) * 1024 * 1024 * 1024
    """
    if [ \$( file -L --mime-type -b ${fasta} | grep "gzip" | wc -l ) -gt 0 ]; then 
        zcat ${fasta} > fasta_unzip.fa
    else 
        mv ${fasta} fasta_unzip.fa
    fi
    if [ \$( file -L --mime-type -b ${gtf} | grep "gzip" | wc -l ) -gt 0 ]; then
        zcat ${gtf} > gtf_unzip.gtf
    else
        mv ${gtf} gtf_unzip.gtf
    fi
    
STAR --runMode genomeGenerate \
--sjdbGTFfile gtf_unzip.gtf \
--genomeFastaFiles fasta_unzip.fa \
--genomeDir ${params.genome_base} \
--limitGenomeGenerateRAM ${mem_int_bytes}
    """
}

process make_mm2_index{
    time '8h'
    
    input:
    file(fasta)
    
    publishDir "${params.output_directory}", mode: 'copy'
    
    output:
    file("${params.genome_base}.mm2")
    
    script:
    """
    minimap2 -d ${params.genome_base}.mm2 ${fasta}
    """
}

workflow{
    if (params.gtf){
        // Build STAR reference
        def star_dat = Channel.fromPath(params.fasta).combine(Channel.fromPath(params.gtf))
        make_star_index(star_dat)               
    }
    // Build minimap2 reference
    def mm2_dat = Channel.fromPath(params.fasta)
    make_mm2_index(mm2_dat)
}

