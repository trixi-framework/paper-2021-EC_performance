# Instructions

This folder contains material to measure memory bandwidth and floating point
operations using [LIKWID](https://github.com/RRZE-HPC/likwid) via it's Julia
interface [LIKWID.jl](https://github.com/JuliaPerf/LIKWID.jl). Thus, you need
Linux and [install LIKWID](https://github.com/RRZE-HPC/likwid#download-build-and-install).

We used LIKWID Version 5.2.0 (commit: 233ab943543480cd46058b34616c174198ba0459)
on Kubuntu 20.04 for the results reported in the article.

The following commands can be run to reproduce all results. The generated data
will be shown in the Julia REPL and saved in `*.dat` files in this directory.

- Vectorization ratio of volume terms:
  ```bash
  julia vectorization_ratio.jl
  ```
- Peak performance of volume terms:
  ```bash
  julia peak_performance.jl
  ```
