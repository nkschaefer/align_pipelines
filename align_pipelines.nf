#! /usr/bin/env nextflow
nextflow.enable.dsl=2

include { align_rna } from './workflows/align_rna.nf'

if (!params.output_directory){
    error("Output directory is required.")
}
if (!params.libs){
    error("Libs file is required.")
}

workflow{
    
    // index input FASTA 
    fa_indexed = faidx(Channel.fromPath(params.genome))
    
    // Toss out short sequences and ones marked excludable (i.e. 
    // random, chrUn, alt scaffolds)
    if (params.annotation){
        
        fa_ann_filt = fa_indexed.combine(Channel.fromPath(params.annotation)) \
            | rm_short_scaffolds_annotation
        rm_numts_annotation(fa_ann_filt)
        fa_filt = fa_ann_filt.map{ x ->
            x[0]
        }
        fa_ann_filt.map{ x -> 
            x[1]
        } | compress_annotation
    }
    else{
        fa_filt = rm_short_scaffolds_fa(fa_indexed)
    }
    
    fa_filt | compress_fa

    numtmasked_fa = mask_numts_fasta(fa_filt)
    bl_bed = genmap_blacklist_fasta(fa_filt)
    numtmasked_fa.combine(bl_bed) | mask_both_beds
    
}


