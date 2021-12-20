
## prepare IO
filename = joinpath(@__DIR__, "vectorization_ratio.dat")
open(filename, "w") do io
  println(io, "# NumericalFlux\tTreeMesh\tStructuredMesh\tP4estMesh")
end


## gather data for volume terms implemented in Trixi.jl
measured_string = read(`likwid-perfctr -C 0 -g FLOPS_DP -m $(Base.julia_cmd()) --check-bounds=no --threads=1 $(joinpath(@__DIR__, "measure_volume_terms.jl"))`, String)

for numerical_flux in ["flux_shima_etal_turbo", "flux_ranocha_turbo"]
  offset = findfirst("Region TreeMesh-$numerical_flux", measured_string) |> last
  m = match(r"(Vectorization ratio\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  ratio_TreeMesh = parse(Float64, m.captures[2])
  @info "TreeMesh" numerical_flux ratio_TreeMesh

  offset = findfirst("Region StructuredMesh-$numerical_flux", measured_string) |> last
  m = match(r"(Vectorization ratio\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  ratio_StructuredMesh = parse(Float64, m.captures[2])
  @info "StructuredMesh" numerical_flux ratio_StructuredMesh

  offset = findfirst("Region P4estMesh-$numerical_flux", measured_string) |> last
  m = match(r"(Vectorization ratio\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  ratio_P4estMesh = parse(Float64, m.captures[2])
  @info "P4estMesh" numerical_flux ratio_P4estMesh

  # save results to disk
  open(filename, "a") do io
    println(io, numerical_flux, "\t", ratio_TreeMesh, "\t", ratio_StructuredMesh, "\t", ratio_P4estMesh)
  end
end
