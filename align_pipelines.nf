#! /usr/bin/env nextflow
nextflow.enable.dsl=2

include { align_rna } from './workflows/align_rna.nf'
include { align_rna_demux_species } from './workflows/align_rna.nf'
include { align_atac } from './workflows/align_atac.nf'
include { align_atac_demux_species} from './workflows/align_atac.nf'

if (!params.output_directory){
    error("Output directory is required.")
}
if (!params.libs){
    error("Libs file is required.")
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
    }
    else{
        def libs = Channel.fromPath(params.libs).splitText().map{ id -> id.trim() }
        if (params.rna_ref){
            align_rna(libs)
        }
        if (params.atac_ref){
            align_atac(libs)
        }
    }
}
