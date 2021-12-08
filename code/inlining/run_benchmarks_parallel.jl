
using DelimitedFiles

@sync begin
  # 2D, flux of Shima et al.
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 2, flux_shima_etal_notinlined, flux_shima_etal)'`)
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 2, flux_shima_etal_notinlined, flux_shima_etal)'`)


  # 3D, flux of Shima et al.
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 3, flux_shima_etal_notinlined, flux_shima_etal)'`)
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 3, flux_shima_etal_notinlined, flux_shima_etal)'`)

  # 2D, flux of Ranocha
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 2, flux_ranocha_notinlined, flux_ranocha)'`)
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 2, flux_ranocha_notinlined, flux_ranocha)'`)

  # 3D, flux of Ranocha
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(3:13, 3, flux_ranocha_notinlined, flux_ranocha)'`)
  @async run(`$(Base.julia_cmd()) --check-bounds=no --threads=1 -e 'include("run_benchmarks.jl"); run_benchmarks(14:15, 3, flux_ranocha_notinlined, flux_ranocha)'`)
end
