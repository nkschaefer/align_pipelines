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

#### Building references
Before you align anything, you will need to build reference data. 

To do this, use the pipeline `make_ref.nf` and copy the `example_make_ref.yml` file. 

You must provide a genome in FASTA format; if you also provide an annotation in GTF format, it will build a [`STAR`](https://github.com/alexdobin/STAR) index for mapping RNA-seq data. Whether or not you provide an annotation, it will also build a [`minimap2`](https://github.com/lh3/minimap2) index for aligning ATAC-seq data. 

For this pipeline, the `genome_base` parameter you provide will serve as the name of the `STAR` genome directory, and the `minimap2` index will be called `[genome_base].mm2`.

`STAR` indexes are not compatible across versions - this pipeline is set up to install and use `STAR` v 2.7.10b. 

#### demux_species output
This pipeline can also run on the output from [`CellBouncer`](https://github.com/nkschaefer/cellbouncer)'s `demux_species` program. For this, copy the `example_demux_species.yml` file, edit parameters as you need, and run `align_pipelines.yml` using this file. Output will be structured similarly to that from `demux_species`, and all species demultiplexing data will be copied into output directories.
