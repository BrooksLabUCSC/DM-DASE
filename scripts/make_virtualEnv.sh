#!/bin/bash

conda create -n "temp_rEnv" r-essentials r-base &&
conda activate temp_rEnv &&

conda install -c bioconda bioconductor-drimseq &&
conda install -c bioconda bioconductor-stager &&
conda install -c conda-forge r-reticulate &&
conda install -c r r-gridextra &&
conda install -c conda-forge argparse &&
deactivate 
