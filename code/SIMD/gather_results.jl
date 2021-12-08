
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using DelimitedFiles
using NaturalSort

# Gather results scattered across multiple files generated
# from independent benchmark jobs.
function gather_results(filename_start)
  filenames = String[]
  for (root, dirs, files) in walkdir(@__DIR__)
    for filename in files
      if startswith(filename, filename_start) && !contains(filename, "relative")
        push!(filenames, joinpath(root, filename))
      end
    end
  end

  sort!(filenames, lt=natural)
  @assert !isempty(filenames)

  filename = first(filenames)
  header = readline(filename)
  data = readdlm(filename, comments=true)
  for filename in Iterators.drop(filenames, 1)
    new_data = readdlm(filename, comments=true)
    data = vcat(data, new_data)
  end

  @assert last(filename_start) == '_'
  filename = joinpath(@__DIR__, filename_start[1:end-1] * ".dat")
  open(filename, "w") do io
    println(io, header)
    writedlm(io, data)
  end

  filename
end


gather_results("pids_2D_flux_shima_etal_turbo_")
gather_results("pids_3D_flux_shima_etal_turbo_")

gather_results("pids_2D_flux_ranocha_turbo_")
gather_results("pids_3D_flux_ranocha_turbo_")
