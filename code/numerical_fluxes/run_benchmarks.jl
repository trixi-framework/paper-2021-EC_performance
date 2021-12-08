
## activate project environment
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()


## load packages
using LinearAlgebra, Printf, Statistics
using BenchmarkTools
using Trixi


## set up benchmark code
function format_result(result, description)
  text = @sprintf("    & \\SI{%4.1f \\pm %.1f}{\\ns}", time(mean(result)), time(std(result)))
  if !isempty(description)
    text *= "  % " * description
  end
  return text
end

@inline function rotated_precomputed_flux(numerical_flux, u_ll, u_rr, norm_, normal_vector, equations)
  u_ll_rotated = Trixi.rotate_to_x(u_ll, normal_vector, equations)
  u_rr_rotated = Trixi.rotate_to_x(u_rr, normal_vector, equations)

  f = numerical_flux(u_ll_rotated, u_rr_rotated, 1, equations)

  return Trixi.rotate_from_x(f, normal_vector, equations) * norm_
end

@inline function rotated_precomputed_flux(numerical_flux, u_ll, u_rr, norm_, normal_vector, tangent1, tangent2, equations)
  u_ll_rotated = Trixi.rotate_to_x(u_ll, normal_vector, tangent1, tangent2, equations)
  u_rr_rotated = Trixi.rotate_to_x(u_rr, normal_vector, tangent1, tangent2, equations)

  f = numerical_flux(u_ll_rotated, u_rr_rotated, 1, equations)

  return Trixi.rotate_from_x(f, normal_vector, tangent1, tangent2, equations) * norm_
end

function run_benchmarks_2d(numerical_flux, write_or_append="w")
  println("Benchmarking ", numerical_flux, " in 2D")

  equations = CompressibleEulerEquations2D(1.4)

  u_ll = SVector(1.0, 0.5, -0.7, 1.0)
  u_rr = SVector(1.5, -0.2, 0.1, 5.0)

  normal_direction = SVector(1.0, 0.0)

  sleep(0.5)
  f = numerical_flux
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $1, $equations)
  format_result(result, "Cartesian") |> println
  open(joinpath(@__DIR__, "Cartesian_2D.txt"), write_or_append) do io
    println(io, format_result(result, "Cartesian"))
  end

  sleep(0.5)
  f = numerical_flux
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $(Ref(normal_direction))[], $equations)
  format_result(result, "Directional") |> println
  open(joinpath(@__DIR__, "Directional_2D.txt"), write_or_append) do io
    println(io, format_result(result, "Directional"))
  end

  sleep(0.5)
  f = FluxRotated(numerical_flux)
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $(Ref(normal_direction))[], $equations)
  format_result(result, "Rotated") |> println
  open(joinpath(@__DIR__, "Rotated_2D.txt"), write_or_append) do io
    println(io, format_result(result, "Rotated"))
  end

  # uses code from `FluxRotated`
  sleep(0.5)
  norm_ = norm(normal_direction)
  normal_vector = normal_direction / norm_
  result = @benchmark rotated_precomputed_flux($numerical_flux, $(Ref(u_ll))[], $(Ref(u_rr))[],
    $(Ref(norm_))[], $(Ref(normal_vector))[], $equations)
  format_result(result, "Rotated (precomputed)") |> println
  open(joinpath(@__DIR__, "Rotated_precomputed_2D.txt"), "w") do io
    println(io, format_result(result, "Rotated (precomputed)"))
  end

  return nothing
end

function run_benchmarks_3d(numerical_flux, write_or_append="w")
  println("Benchmarking ", numerical_flux, " in 3D")

  equations = CompressibleEulerEquations3D(1.4)

  u_ll = SVector(1.0, 0.5, -0.7, 0.1, 1.0)
  u_rr = SVector(1.5, -0.2, 0.1, -0.3, 5.0)

  normal_direction = SVector(1.0, 0.0, 0.0)

  sleep(0.5)
  f = numerical_flux
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $1, $equations)
  format_result(result, "Cartesian") |> println
  open(joinpath(@__DIR__, "Cartesian_3D.txt"), write_or_append) do io
    println(io, format_result(result, "Cartesian"))
  end

  sleep(0.5)
  f = numerical_flux
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $(Ref(normal_direction))[], $equations)
  format_result(result, "Directional") |> println
  open(joinpath(@__DIR__, "Directional_3D.txt"), write_or_append) do io
    println(io, format_result(result, "Directional"))
  end

  sleep(0.5)
  f = FluxRotated(numerical_flux)
  result = @benchmark $f($(Ref(u_ll))[], $(Ref(u_rr))[], $(Ref(normal_direction))[], $equations)
  format_result(result, "Rotated") |> println
  open(joinpath(@__DIR__, "Rotated_3D.txt"), write_or_append) do io
    println(io, format_result(result, "Rotated"))
  end

  # uses code from `FluxRotated`
  sleep(0.5)
  norm_ = norm(normal_direction)
  normal_vector = normal_direction / norm_
  tangent1 = SVector(normal_direction[2], normal_direction[3], -normal_direction[1])
  tangent1 -= dot(normal_vector, tangent1) * normal_vector
  tangent1 = normalize(tangent1)
  tangent2 = normalize(cross(normal_direction, tangent1))
  result = @benchmark rotated_precomputed_flux($numerical_flux, $(Ref(u_ll))[], $(Ref(u_rr))[],
    $(Ref(norm_))[], $(Ref(normal_vector))[], $(Ref(tangent1))[], $(Ref(tangent2))[], $equations)
  format_result(result, "Rotated (precomputed)") |> println
  open(joinpath(@__DIR__, "Rotated_precomputed_3D.txt"), write_or_append) do io
    println(io, format_result(result, "Rotated (precomputed)"))
  end

  return nothing
end


## run benchmarks
run_benchmarks_2d(flux_shima_etal, "w")
run_benchmarks_2d(flux_ranocha, "a")
run_benchmarks_2d(flux_lax_friedrichs, "a")
run_benchmarks_2d(flux_hll, "a")

run_benchmarks_3d(flux_shima_etal, "w")
run_benchmarks_3d(flux_ranocha, "a")
run_benchmarks_3d(flux_lax_friedrichs, "a")
run_benchmarks_3d(flux_hll, "a")

