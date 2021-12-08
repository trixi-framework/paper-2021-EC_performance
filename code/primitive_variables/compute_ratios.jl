
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## compute relative performance
using DelimitedFiles
using Measurements

function compute_relative_performance(filename)
  header = "# Polydeg Primitive/conservative-variables-mean std"
  data = readdlm(filename, comments=true)
  cons = data[:, 2] .± data[:, 3]
  prim = data[:, 4] .± data[:, 5]
  ratio = prim ./ cons
  open(filename[1:end-4] * "_relative.dat", "w") do io
    println(io, header)
    writedlm(io, hcat(data[:, 1], # polydeg
                      Measurements.value.(ratio),
                      Measurements.uncertainty.(ratio)))
  end
end

compute_relative_performance(joinpath(@__DIR__, "primitive_2D_flux_shima_etal.dat"))
compute_relative_performance(joinpath(@__DIR__, "primitive_3D_flux_shima_etal.dat"))
compute_relative_performance(joinpath(@__DIR__, "primitive_2D_flux_ranocha.dat"))
compute_relative_performance(joinpath(@__DIR__, "primitive_3D_flux_ranocha.dat"))
