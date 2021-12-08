#!/bin/bash
#PBS -N elixir_benchmarks
#PBS -l select=1:node_type=clx-25
#PBS -l walltime=02:00:00

# Switch to submission directory
cd $PBS_O_WORKDIR

time julia --check-bounds=no --threads=1 run_elixir_benchmarks.jl
