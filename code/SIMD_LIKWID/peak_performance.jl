
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using LIKWID


## measure optimistic peakflops (AVX2 FMA or AVX512 FMA if available)
L1_cache_size = LIKWID.get_cpu_topology().cacheLevels[1].size ÷ 1024 # in kB

cpuinfo = LIKWID.get_cpu_info()
if occursin("AVX512", cpuinfo.features)
  likwid_bench_kernel = "peakflops_avx512_fma"
elseif occursin("AVX2", cpuinfo.features)
  likwid_bench_kernel = "peakflops_avx_fma"
else
  likwid_bench_kernel = "peakflops_sse_fma"
end

max_flops_string = read(`likwid-bench -t $likwid_bench_kernel -W N:$(L1_cache_size)kB:1`, String)
max_flops = parse(Float64, match(r"(MFlops/s:\s+)(\d+\.\d+)", max_flops_string).captures[2]) / 1024
@info "PeakPerformance" max_flops


## prepare IO
filename = joinpath(@__DIR__, "peak_performance.dat")
open(filename, "w") do io
  println(io, "# NumericalFlux\tTreeMesh\tStructuredMesh\tP4estMesh\tPeakPerformance")
end


## gather data for volume terms implemented in Trixi.jl
measured_string = read(`likwid-perfctr -C 0 -g MEM_DP -m $(Base.julia_cmd()) --check-bounds=no --threads=1 $(joinpath(@__DIR__, "measure_volume_terms.jl"))`, String)

for numerical_flux in ["flux_shima_etal", "flux_ranocha"]
  offset = findfirst("Region TreeMesh-$numerical_flux", measured_string) |> last
  m = match(r"(DP \[MFLOP/s\]\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  flops_TreeMesh = parse(Float64, m.captures[2]) / 1024
  @info "TreeMesh" numerical_flux flops_TreeMesh

  offset = findfirst("Region StructuredMesh-$numerical_flux", measured_string) |> last
  m = match(r"(DP \[MFLOP/s\]\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  flops_StructuredMesh = parse(Float64, m.captures[2]) / 1024
  @info "StructuredMesh" numerical_flux flops_StructuredMesh

  offset = findfirst("Region P4estMesh-$numerical_flux", measured_string) |> last
  m = match(r"(DP \[MFLOP/s\]\s+\|\s+)(\d+\.\d+)", measured_string, offset)
  flops_P4estMesh = parse(Float64, m.captures[2]) / 1024
  @info "P4estMesh" numerical_flux flops_P4estMesh

  # save results to disk
  open(filename, "a") do io
    println(io, numerical_flux, "\t",
                flops_TreeMesh, "\t",
                flops_StructuredMesh, "\t",
                flops_P4estMesh, "\t",
                max_flops)
    println(io, numerical_flux, "\t",
                flops_TreeMesh / max_flops, "\t",
                flops_StructuredMesh / max_flops, "\t",
                flops_P4estMesh / max_flops, "\t",
                max_flops / max_flops)
  end
end
