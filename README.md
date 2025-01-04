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
Copy either `example.yml` or `example_demux_species.yml` (depending on whether you want to map data from the sequencer or data output by CellBouncer's `demux_species`) and edit the file to reflect your parameters. Then run the pipeline as follows:
```
nextflow [/path/to/align_pipelines]/align_pipelines.nf -params-file [params.yml]
```
where `[params.yml]` is the file you edited.

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


