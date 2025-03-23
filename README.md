# align_pipelines
Nextflow pipelines for single cell sequencing data alignment

## Installation

### Prerequisites
Make sure you have [`conda`](https://github.com/conda-forge/miniforge)

```
git clone git@github.com:nkschaefer/align_pipelines.git
cd align_pipelines
mamba env create --file=environment.yml
conda activate align_pipelines
make
```

Be sure to activate the conda environment (`conda activate align_pipelines`) before running the pipeline.

## Running

### Cluster settings
If you want to run this on a cluster, create a file in the same directory you plan to run from called `nextflow.config`. This file should contain cluster-specific settings. If you are using UCSF's Wynton cluster, for example, the file should contain these lines:
```
process.executor = "sge"
process.clusterOptions = "-V -S /bin/bash"
process.penv = "smp"
```

### Job settings
Copy the relevant `example[whatever].yml` file (depending on what you want to do) to your working directory and edit the file to reflect your parameters. These files contain comments explaining what each setting in them is. Then run the pipeline as follows:
```
nextflow [/path/to/align_pipelines]/[pipeline].nf -params-file [params.yml]
```
where `[pipeline].nf` is the `.nf` file for the pipeline you want to run and `[params.yml]` is the file you edited.

To run on a cluster, you will probably want to submit as a cluster job. For example, on UCSF's Wynton cluster, you could write a script like this:
```
#! /usr/bin/env bash
#$ -cwd
#$ -l h_rt=36:00:00
#$ -V
conda activate align_pipelines
nextflow [/path/to/align_pipelines]/align_pipelines.nf -params-file [params.yml]
```
and then submit the script using `qsub [script.sh]`

If something goes wrong, you can often pick up where you left off by fixing the issue, running again, and adding the `-resume` option to your `nextflow` command.

When the pipeline finishes, all results will be copied to the `output_directory` that you specify in the `.yml` file. `nextflow` stores all temporary data in a directory called `work` located wherever you launched the pipeline from. This directory can get very big, and everything you need to keep will be copied to the `output_directory` on completion. Therefore, after you finish running the pipeline, you can wipe out temporary stuff by running
```
rm -r work
```
in whatever directory you ran the pipeline.

#### Building references -- RNA-seq/ATAC-seq only
Before you align RNA-seq/ATAC-seq, you will need to build reference data. 

To do this, use the pipeline `make_ref.nf` and copy the `example_make_ref.yml` file. 

You must provide a genome in FASTA format; if you also provide an annotation in GTF format, it will build a [`STAR`](https://github.com/alexdobin/STAR) index for mapping RNA-seq data. Whether or not you provide an annotation, it will also build a [`minimap2`](https://github.com/lh3/minimap2) index for aligning ATAC-seq data. 

`STAR` will usually complain if input files are gzipped, but this pipeline can handle gzipped or uncompressed FASTA and GTF files.

For this pipeline, the `genome_base` parameter you provide will serve as the name of the `STAR` genome directory, and the `minimap2` index will be called `[genome_base].mm2`.

`STAR` indexes are not compatible across versions - this pipeline is set up to install and use `STAR` v 2.7.10b. 

#### demux_species output
This pipeline can also run on the output from [`CellBouncer`](https://github.com/nkschaefer/cellbouncer)'s `demux_species` program. For this, copy the `example_demux_species.yml` file, edit parameters as you need, and run `align_pipelines.yml` using this file. Output will be structured similarly to that from `demux_species`, and all species demultiplexing data will be copied into output directories.

#### Barcode whitelists
To process ATAC-seq or RNA-seq data, this pipeline requires [allowed barcode lists](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist) to determine which reads are associated with cell barcodes. If you are aligning RNA-seq only, there will be one list (`rna_whitelist` parameter). If you are processing scATAC-seq data only, there will also be one list (`atac_whitelist` parameter). If you are processing multiome data (where RNA-seq and ATAC-seq are collected from the same cells and have corresponding barcodes), there are two whitelists (`rna_whitelist` and `atac_whitelist`). These files are set up so that the barcode on line N of one file corresponds to the barcode on the same line of the other file. ATAC-seq reads contain the barcodes in `atac_whitelist`, but these are then converted into the corresponding barcodes from `rna_whitelist` to be reported in output data.

### What if you have more than one directory containing reads of a specific type?
The pipelines expect reads of any given type (e.g. RNA-seq, ATAC-seq, or genomic DNA) to be located in one directory. This isn't always the case. If your reads are in multiple directories, that's okay, as long as the reads don't have identical names. If they do (e.g. if you got more data back from the same exact sequencing run a second time), you'll need to concatenate the files.

Otherwise, you can do this (assuming you have two directories, `/path/to/dir1` and `/path/to/dir2` containing the same type of data:

```
# Create a new directory
mkdir dir3

# Add symlinks to all reads in dir1 to new directory
ls /path/to/dir1/*.fastq.gz | while read filename; do
    ln -s $( realpath $filename ) dir3
done

# Add symlinks to all reads in dir2 to new directory
ls /path/to/dir2/*.fastq.gz | while read filename; do
    ln -s $( realpath $filename ) dir3
done
```

Now you have a new directory, `dir3`, that contains symbolic links to all files in `dir1` and `dir2`. You can just provide `dir3` in the YAML file as the source of reads for your data type.
