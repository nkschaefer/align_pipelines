#! /usr/bin/env nextflow

process map_rna{
    cpus params.threads
    time '36h'
    memory params.memgb + ' GB'
    
    input: 
    tuple val(lib),
        val(rna_idx_base), 
        file(rna_idx_files), 
        file(whitelist), 
        file(reads1), 
        file(reads2)

    publishDir "${params.out}/lib", mode: 'copy', saveAs: { path -> path }

    output:
    tuple val(lib),
        file("*_gex.bam"),
        file("*_gex.bam.bai"),
        file("Barcodes.stats"),
        file("Features.stats"),
        file("Summary.csv"),
        file("UMIperCellSorted.txt"),
        file("raw/*"),
        file("filtered/*")

    script:
    """
    STAR --genomeDir ${rna_idx_base} \
     --readFilesIn ${reads2} ${reads1} \
     --readFilesCommand zcat\
     --clipAdapterType CellRanger4\
     --outBAMsortingThreadN 1 \
     --outFileNamePrefix ${lib} \
     --outSAMattributes NH HI AS nM CR CY UR UY GX GN CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --soloType CB_UMI_Simple \
     --soloCBstart 1 \
     --soloCBlen 16 \
     --soloUMIstart 17 \
     --soloUMIlen 12 \
     --soloCBwhitelist ${whitelist} \
     --outFilterScoreMin 30 \
     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
     --soloUMIfiltering MultiGeneUMI_CR \
     --soloUMIdedup 1MM_CR \
     --soloCellFilter EmptyDrops_CR \
     --soloBarcodeReadLength 0 \
     --limitSjdbInsertNsj 5000000 \
     --soloFeatures GeneFull_Ex50pAS \
     --soloMultiMappers EM 
    
    mv ${lib}Aligned.sortedByCoord.out.bam ${lib}_gex.bam
    samtools index ${lib}_gex.bam
    mv ${lib}Solo.out/Barcodes.stats Barcodes.stats
    cp ${lib}Solo.out/GeneFull_Ex50pAS/Features.stats .
    cp ${lib}Solo.out/GeneFull_Ex50pAS/Summary.csv .
    cp ${lib}Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt .
    gzip ${lib}Solo.out/GeneFull_Ex50pAS/filtered/*
    gzip ${lib}Solo.out/GeneFull_Ex50pAS/raw/*
    mv ${lib}Solo.out/GeneFull_Ex50pAS/filtered .
    mv ${lib}Solo.out/GeneFull_Ex50pAS/raw .
    """
}

workflow align_rna{
    take:
    val(lib)
    
    main:
    
    if (!params.rna_ref){
        error("RNA reference is required")
    }
    if (!params.rna_whitelist){
        error("RNA whitelist is required")
    }
    if (!params.rna_dir){
        error("RNA directory is required")
    }

    rna_idx = Channel.frompath("${params.rna_dir}/**")
    rna_idx.view()

}
