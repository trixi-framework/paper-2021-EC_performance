# Instructions

The following commands can be run to reproduce all results. The generated data
will be shown in the Julia REPL.

```bash
julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl")'
```

In contrast to other benchmark suites contained in this repository, this one is
significantly less expensive and will finish in a short amount of time.
