#!/bin/bash
#PBS -N primitive_variables
#PBS -l select=1:node_type=clx-25
#PBS -l walltime=00:30:00

# Switch to submission directory
cd $PBS_O_WORKDIR

time julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl")'
