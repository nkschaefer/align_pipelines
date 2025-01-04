#! /usr/bin/env nextflow
nextflow.enable.dsl=2

include { align_rna } from './workflows/align_rna.nf'
include { align_rna_demux_species } from './workflows/align_rna.nf'
include { align_atac } from './workflows/align_atac.nf'
include { align_atac_demux_species} from './workflows/align_atac.nf'

if (!params.output_directory){
    error("Output directory is required.")
}

process copy_species_files{
    input:
    tuple val(lib),
    file(sfile)
    
    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: {path -> path }
    
    output:
    file("${sfile}")    
    
    script:
    """
    """
}

workflow{
    if (params.demux_species){
        if (params.rna_ref_species){
            align_rna_demux_species()
        } else if (params.rna_ref){
            error("If using demux_species output, you must provide rna_ref_species instead of rna_ref")
        }
        if (params.atac_ref_species){
            align_atac_demux_species()
        } else if (params.atac_ref){
            error("If using demux_species output, you must provide atac_ref_species instead of atac_ref")
        }
        
        def species_files = Channel.fromPath("${params.demux_species}/*/species*").map{ fn -> 
            spec = fn.toString().split('/')[-2]
            return [ spec, fn ]
        }
        copy_species_files(species_files)
    
    }
    else{
        if (!params.libs){
            error("Libs file is required.")
        }
        def libs = Channel.fromPath(params.libs).splitText().map{ id -> id.trim() }
        if (params.rna_ref){
            align_rna(libs)
        }
        if (params.atac_ref){
            align_atac(libs)
        }
    }
}
