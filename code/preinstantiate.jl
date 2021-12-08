function preinstantiate(path=".")
  benchmark_dirs = ["Cartesian_vs_curved", "Gauss", "inlining", "numerical_fluxes",
                    "overintegration", "precompute_logs", "primitive_variables"]
  current_dir = pwd()

  for dir in benchmark_dirs
    cd(joinpath(path, dir))
    println("#"^100)
    println("### Instantiating $dir... ###")
    run(`julia -e 'import Pkg; Pkg.activate("."); Pkg.instantiate()'`)
    cd(current_dir)
  end
  println("#"^100)

  cd(current_dir)
end
