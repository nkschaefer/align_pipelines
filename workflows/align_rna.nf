#! /usr/bin/env nextflow

def rna_ref_map = [:]
if (params.rna_ref_species){
    new File(params.rna_ref_species).eachLine { line ->
        def (species, ref) = line.split('\t')
        rna_ref_map[species] = ref
    }
}

process map_rna{
    cpus params.threads
    time '36h'
    memory params.memgb + ' GB'
    
    input: 
    tuple val(lib),
        file(reads1),
        file(reads2),
        val(refname),
        file(ref),
        file(whitelist)

    publishDir "${params.output_directory}/${lib}", mode: 'copy', saveAs: { path -> path }

    output:
    tuple val(lib),
        file("gex.bam"),
        file("gex.bam.bai"),
        file("Barcodes.stats"),
        file("Features.stats"),
        file("Summary.csv"),
        file("UMIperCellSorted.txt"),
        file("raw/*"),
        file("filtered/*")

    script:
    def reftrunc = refname.split('/')[-1]    
    def r1 = reads1.join(',')
    def r2 = reads2.join(',')
    def lib2 = lib.replace('/', '_')
    """
    if [ ! -d ${reftrunc} ]; then
        mkdir ${reftrunc}
        mv ${ref} ${reftrunc}
    fi

STAR --genomeDir ${reftrunc} \
 --runThreadN ${params.threads} \
 --readFilesIn ${r2} ${r1} \
 --readFilesCommand zcat \
 --clipAdapterType CellRanger4 \
 --outBAMsortingThreadN 1 \
 --outFileNamePrefix ${lib2} \
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
    mv ${lib2}Aligned.sortedByCoord.out.bam gex.bam
    samtools index gex.bam
    mv ${lib2}Solo.out/Barcodes.stats Barcodes.stats
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/Features.stats .
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/Summary.csv .
    cp ${lib2}Solo.out/GeneFull_Ex50pAS/UMIperCellSorted.txt .
    gzip ${lib2}Solo.out/GeneFull_Ex50pAS/filtered/*
    gzip ${lib2}Solo.out/GeneFull_Ex50pAS/raw/*
    mv ${lib2}Solo.out/GeneFull_Ex50pAS/filtered .
    mv ${lib2}Solo.out/GeneFull_Ex50pAS/raw .
    """
}

workflow align_rna_demux_species{
    main:
    
    if (!params.rna_ref_species){
        error("rna_ref_species is required")
    }
    if (!params.demux_species){
        error("demux_species output dir is required")
    }
    if (!params.rna_whitelist){
        error("rna_whitelist is required")
    }
    
    def read_pairs = Channel.fromFilePairs("${params.demux_species}/*/*/GEX_*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        def libn = reads[0].toString().split('/')[-3]
        def species = reads[0].toString().split('/')[-2]
        if (! rna_ref_map[species]){
            error("Species " + species + " does not have an RNA-seq reference specified")
        }
        [libn + "/" + species, reads, rna_ref_map[species], file(rna_ref_map[species] + "/*") ]
    }.groupTuple().map{ lib, reads, refname, ref ->
        def r1s = []
        def r2s = []
        for (elt in reads){
            r1s.add(elt[0])
            r2s.add(elt[1])    
        }
        return [lib, r1s, r2s, refname[0], ref[0] ]
    }.combine(Channel.fromPath(params.rna_whitelist))
    map_rna(read_pairs)
}

workflow align_rna{
    take:
    libs
    
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
    
    def rna_idx = Channel.fromPath("${params.rna_ref}/**").collect().map{ x -> 
        [params.rna_ref, x]}
    
    def read_pairs = Channel.fromFilePairs("${params.rna_dir}/*_S*L*_R{1,2}*.fastq.gz").map{ id, reads ->
        [id.replaceFirst(/_S(\d+)_L(\d+)/, ''), reads]
    }.groupTuple().map{ lib, reads ->
        def r1s = []
        def r2s = []
        for (elt in reads){
            r1s.add(elt[0])
            r2s.add(elt[1])    
        }
        return [lib, r1s, r2s]
    }
    
    map_rna(libs.cross(read_pairs).map{ lib, tup -> tup }.combine(rna_idx).combine(Channel.fromPath(params.rna_whitelist)))
    
}
