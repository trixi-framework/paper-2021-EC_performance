# Instructions

The following two commands can be run independently of each other to reproduce
all results. The generated data will be stored in files `Gauss_*.dat` in this
directory.

- 2D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 2, flux_ranocha, flux_ranocha)'
  ```
- 3D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 3, flux_ranocha, flux_ranocha)'
  ```
- Relative timings
  ```bash
  julia --check-bounds=no --threads=1 run_elixir_benchmarks.jl
  ```


You can also split the commands into multiple independent parts, e.g.,

- 2D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 2, flux_ranocha, flux_ranocha)'
  ```
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 2, flux_ranocha, flux_ranocha)'
  ```
- 3D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 3, flux_ranocha, flux_ranocha)'
  ```
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 3, flux_ranocha, flux_ranocha)'
  ```
- Relative timings
  ```bash
  julia --check-bounds=no --threads=1 run_elixir_benchmarks.jl
  ```

In that case, you can execute
```bash
julia gather_results.jl
```
to gather results from individual runs.
