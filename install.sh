#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda env create -f environment.yml
conda activate hocort-eval && 
pip install git+https://github.com/ignasrum/hocort.git@1.0.0 &&
pip install git+https://github.com/ignasrum/fastq-tools.git@0.1.0
