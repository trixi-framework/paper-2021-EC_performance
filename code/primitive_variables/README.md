# Instructions

The following command can be run to reproduce all results. The generated data
will be stored in files `primitive_*.dat` in this directory.

```bash
julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl")'
```

Afterwards, you should execute
```bash
julia compute_ratios.jl
```
to postprocess the results.
