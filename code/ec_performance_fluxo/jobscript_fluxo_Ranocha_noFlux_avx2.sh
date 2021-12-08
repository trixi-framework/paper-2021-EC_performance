#!/bin/bash
#PBS -N fluxo_Ranocha_noFlux_avx2
#PBS -l select=1:node_type=clx-25
#PBS -l walltime=06:00:00

# Switch to submission directory
cd $PBS_O_WORKDIR

set -xo pipefail

module load tools/cmake
module load tools/hdf5/1.12.0-openmpi-4.1.1-intel-19.1.3
module load bigdata/conda/miniconda-4.10.3

source ~/.bashrc
conda activate my-env

cd serialTests
time python3 runTests.py  --compiler ifort --id Ranocha_noFlux_avx2 --flux 32 --no-precompileFlux
