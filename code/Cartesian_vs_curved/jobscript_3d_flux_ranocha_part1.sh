#!/bin/bash
#PBS -N 3d_flux_ranocha_part1
#PBS -l select=1:node_type=clx-25
#PBS -l walltime=08:00:00

# Switch to submission directory
cd $PBS_O_WORKDIR

time julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 3, flux_ranocha, flux_ranocha)'
