# Running EC performance tests with FLUXO

## Requirements
* Fortran compiler
* HDF5 with fortran support
* CMake (minimum required version is 3.5.1)
* python (>=3) with numpy

## Instructions

To run the performance tests, follow the instructions:

* Move to this directory and Clone the fluxo repository:
  ```
  cd 2021_EC_performance/code/ec_performance_fluxo
  git clone git@github.com:project-fluxo/fluxo.git
  ```
* Move to the branch of fluxo were the performance improvements are implemented:
  ```
  cd fluxo
  git checkout performance
  cd ..
  ```
* *Optional:* For the tests, you can use the Cartesian mesh with 8² elements under `MeshGeneration/CartBoxPeriodic_008_008_008_mesh.h5`. If you prefer to generate the mesh, you can do it using [hopr](https://gitlab.com/project-fluxo/hopr):
  ```
  cd MeshGeneration
  hopr parameter_hopr.ini
  cd ..
  ```
* Go to the folder `serialTests` and run the tests using the python script:
  ```
  cd serialTests
  python3 runTests.py [options]
  ```
  Following optional arguments are available:
  ```
  -h, --help           show this help message and exit
  --compiler COMPILER  Fortran compiler for FLUXO. If empty, the Fortran compiler found by cmake is used.
  --id ID              Test Identifier. Used for the build directories ("build/build_Euler_<id>...") and the results files ("serialTests/results/results_<id>_.dat")
  --flux FLUX          Volume and surface flux used for the tests.
                           32: [DEFAULT] Ranocha flux
                           40: Shima et al. flux
  --precompileN        [DEFAULT] Compile fluxo with fixed polynomial degree.
  --no-precompileN     [optional] deactivate --precompileN.
  --precompileFlux     [DEFAULT] Compile fluxo with volume flux.
  --no-precompileFlux  [optional] deactivate --precompileN.
  --extra_options EXTRA_OPTIONS
                        Additional options for cmake. For example, to compile with Gauss nodes use --extra_options=" -DFLUXO_DISC_NODETYPE=GAUSS"

  ```
  
  In particular, the tests for the EC performance can be launched with the following commands:
  ```
  python3 runTests.py  --compiler ifort --id Ranocha --flux 32
  python3 runTests.py  --compiler ifort --id Ranocha_noN_noFlux --flux 32 --no-precompileN --no-precompileFlux
  python3 runTests.py  --compiler ifort --id Ranocha_noN --flux 32 --no-precompileN
  python3 runTests.py  --compiler ifort --id Ranocha_noFlux --flux 32 --no-precompileFlux
  python3 runTests.py  --compiler ifort --id Shima --flux 40
  python3 runTests.py  --compiler ifort --id Ranocha_Gauss --flux 32 --extra_options=" -DFLUXO_DISC_NODETYPE=GAUSS"
  ```
    
All results were obtained with [this version of fluxo](https://github.com/project-fluxo/fluxo/tree/01b8bfc62a566dde12769989c87630754deab7e1).
