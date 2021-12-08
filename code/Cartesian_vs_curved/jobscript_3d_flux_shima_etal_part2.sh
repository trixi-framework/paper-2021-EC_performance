#!/bin/bash
#PBS -N 3d_flux_shima_etal_part2
#PBS -l select=1:node_type=clx-25
#PBS -l walltime=04:00:00

# Switch to submission directory
cd $PBS_O_WORKDIR

time julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 3, flux_shima_etal, flux_shima_etal)'
