# Running benchmarks

Documentation for how benchmarks were run on HLRS' Vulcan cluster.

## Job submission
All jobscripts are named `jobscript_...sh` and should be submitted using `qsub`
from the folder in which they are located. For example, for `jobscript.sh`,
execute
```bash
qsub jobscript.sh
```
The Trixi.jl tests can all be submitted simultaneously, while the FLUXO results
need to be submitted serially (or a separate `serialTests` directory needs to be
created for each run).

## Trixi.jl

To pre-instantiate all benchmark directories (since we have no internet
connection on the compute nodes), run the following command first:
```bash
 julia -e 'include("preinstantiate.jl"); preinstantiate()'
```
All resulting files for each run were gathered in the corresponding `run_XXX`
folder. Previous result files have **not been modified**.

For all tests, the pre-built binaries for Julia 1.7.0 were used.

To postprocess results from Vulcan and obtain the same directory structure as
from a local run following the instructions in the individual README.md files,
execute
```bash
julia postprocess_results_from_vulcan.jl
```

## FLUXO
All resulting *auxiliary* files for each run were gathered in the corresponding
`run_XXX` folder in the `serialTests` directory. Previous  result files **have
been overwritten**.

The following configuration was used for FLUXO:
* Source: https://github.com/project-fluxo/fluxo
* Branch: `performance`
* Commit hash: `01b8bfc62a566dde12769989c87630754deab7e1`
* Intel Fortran Compiler 19.1.3
