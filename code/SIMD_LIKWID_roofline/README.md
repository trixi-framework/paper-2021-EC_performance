# Empirical roofline model

This folder contains material to measure memory bandwidth and floating point
operations using [LIKWID](https://github.com/RRZE-HPC/likwid) via it's Julia
interface [LIKWID.jl](https://github.com/JuliaPerf/LIKWID.jl). Thus, you need
Linux and [install LIKWID](https://github.com/RRZE-HPC/likwid#download-build-and-install).

We used LIKWID Version 5.2.0 (commit: 233ab943543480cd46058b34616c174198ba0459)
on Kubuntu 20.04 for the results reported in the article.

Execute

```bash
julia --check-bounds=no --threads=1 construct_roofline_model.jl
```

The results are saved in `results.txt`.
