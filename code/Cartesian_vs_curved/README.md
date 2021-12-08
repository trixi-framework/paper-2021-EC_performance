# Instructions

The following four commands can be run independently of each other to reproduce
all results. The generated data will be stored in files `pids_*.dat` in this
directory.

- 2D, flux of Shima et al.
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 2, flux_shima_etal, flux_shima_etal)'
  ```
- 3D, flux of Shima et al.
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 3, flux_shima_etal, flux_shima_etal)'
  ```
- 2D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 2, flux_ranocha, flux_ranocha)'
  ```
- 3D, flux of Ranocha
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(polydegs, 3, flux_ranocha, flux_ranocha)'
  ```


You can also split the commands into multiple independent parts, e.g.,

- 2D, flux of Shima et al.
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 2, flux_shima_etal, flux_shima_etal)'
  ```
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 2, flux_shima_etal, flux_shima_etal)'
  ```
- 3D, flux of Shima et al.
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 3, flux_shima_etal, flux_shima_etal)'
  ```
  ```bash
  julia --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 3, flux_shima_etal, flux_shima_etal)'
  ```
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

In that case, you can execute
```bash
julia gather_results.jl
```
to gather results from individual runs.
