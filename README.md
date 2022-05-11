# hocort-eval

Benchmarking and evaluation of the [HoCoRT](https://github.com/ignasrum/hocort) Python package.
Illumina HiSeq and MiSeq, and Nanopore simulated datasets (1% and 50% host contamination) are generated using InSilicoSeq and NanoSim.
The different HoCoRT pipelines are then tested and benchmarked.

First create the conda environment and activate it:
```
conda env create -f environment.yaml
```
```
conda activate hocort-eval
```
Run the evaluation tests by running the following command:
```
snakemake --cores 16
```
The results will be available in the directories named hiseq_reads, miseq_reads and nanopore_reads.
