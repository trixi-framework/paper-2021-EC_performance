# Instructions

The following command can be run to reproduce all results. The generated data
will be stored in files `overintegration_*.dat` in this directory.

```bash
julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl")'
```
