
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## compute relative performance of inlining
using DelimitedFiles
using Measurements

function read_data(filename)
  header = readline(filename)
  data = readdlm(filename, comments=true)
  polydegs = data[:, 1]
  results = Vector{Measurement{eltype(data)}}[]
  col = 2
  while col + 1 <= size(data, 2)
    push!(results, data[:, col] .± data[:, col+1])
    col += 2
  end

  header, polydegs, results
end

function compute_relative_performance(file_new, file_baseline)
  header, polydegs, results_baseline = read_data(file_baseline)
  header_new, polydegs_new, results_new = read_data(file_new)
  @assert header == header_new
  @assert polydegs == polydegs_new
  @assert size(results_baseline) == size(results_new)

  data = zeros(length(polydegs_new),
               1 + 2 * length(results_baseline))
  data[:, 1] = polydegs
  for i in eachindex(results_baseline, results_new)
    ratio = results_new[i] ./ results_baseline[i]
    data[:, 2 * i] = Measurements.value.(ratio)
    data[:, 2 * i + 1] = Measurements.uncertainty.(ratio)
  end

  open(file_baseline[1:end-4] * "_relative.dat", "w") do io
    println(io, header)
    writedlm(io, data)
  end
end


setup = "2D_flux_shima_etal"
compute_relative_performance(
  joinpath(dirname(@__DIR__), "Cartesian_vs_curved", "pids_$(setup).dat"),
  joinpath(@__DIR__, "pids_$(setup)_notinlined.dat"))

setup = "2D_flux_ranocha"
compute_relative_performance(
  joinpath(dirname(@__DIR__), "Cartesian_vs_curved", "pids_$(setup).dat"),
  joinpath(@__DIR__, "pids_$(setup)_notinlined.dat"))

setup = "3D_flux_shima_etal"
compute_relative_performance(
  joinpath(dirname(@__DIR__), "Cartesian_vs_curved", "pids_$(setup).dat"),
  joinpath(@__DIR__, "pids_$(setup)_notinlined.dat"))

setup = "3D_flux_ranocha"
compute_relative_performance(
  joinpath(dirname(@__DIR__), "Cartesian_vs_curved", "pids_$(setup).dat"),
  joinpath(@__DIR__, "pids_$(setup)_notinlined.dat"))
