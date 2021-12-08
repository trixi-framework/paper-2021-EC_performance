#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to run performance tests with fluxo
Created on Thu Oct 28 14:29:59 2021

@author: andresrueda
"""
import os
import numpy as np
import argparse
import sys

#==========================================================================
# open the file (default "std.out") and parse for `findstr`, 
# in the last `n_tail` lines, 
# returns if found and line of the file
#==========================================================================
def find_last_occurence(findstr,n_tail=40,stdout_filepath="std.out",**kwargs):
   with open(stdout_filepath,'r') as fp: 
      lines=fp.readlines()
   for line in reversed(lines[-n_tail:]) :
      if findstr in line : 
         return True,line
   return False, ' ' 

def read_PID(stdout_filepath="std.out"):
    [found,line] = find_last_occurence("CALCULATION TIME PER TSTEP/DOF",stdout_filepath=stdout_filepath)
    assert found, ('read_PID: did not find PID in '+stdout_filepath)
    
    fluxoPID = float((line.split("[")[1]).split("sec ]")[0])
    nStages = int(line.split("nRKstages:")[1])
    
    return fluxoPID/nStages

"""
Compiles fluxo. Must be called from serialTests folder...
"""
def compileFluxo(precompileFlux, volFlux, build_dir="build_Euler",extra_options=" ", compiler=""):
    
    # Create build directory
    exec_path = "../builds/"+build_dir
    os.system("mkdir -p "+exec_path)
    
    #Move to that directory
    os.chdir(exec_path)
    
    # Compile the code with the desired options
    print("Configuring and compiling fluxo")
    options = " -DCMAKE_BUILD_TYPE=Release -DFLUXO_ENABLE_MPI=NO "+extra_options
    if compiler != "":
        options += " -DCMAKE_Fortran_COMPILER=" +compiler +" "
    
    if precompileFlux:
        options += " -DFLUXO_EQN_VOLFLUX="+str(volFlux)
    os.system("cmake ../../fluxo/ "+options+"  > std_config.outerr 2>&1")
    os.system("make -j > std_compile.outerr 2>&1")
    print("Done!... ")
    
    # Go back to working directory and update exec_path
    os.chdir("../../serialTests/")
    exec_path += "/bin/fluxo"
    
    return exec_path

def RunPerformanceTest(test_id, polyDeg_ini = 3, polyDeg_end = 15, volFlux=32, precompileN = True, precompileFlux = True, extra_options="", compiler=""):
    
    # Create results folder
    os.system("mkdir -p results")
    output_file_name = 'results/results_'+test_id+'.dat'
    
    # run tests
    with open(output_file_name, "w") as output_file:
        output_file.write("# polyDeg ,   meanPID ,    stdDevPID\n")
    
    # Create builds directory
    os.system("mkdir -p ../builds")
    
    # Build code if precompileN = False
    if not precompileN:
        build_dir="build_Euler_"+test_id
        exec_path = compileFluxo(precompileFlux, volFlux, extra_options = extra_options, compiler=compiler, build_dir=build_dir)
    
    # Define parameter file
    prm_path = "parameter_"+test_id+".ini"
    
    polyDegs = range(polyDeg_ini,polyDeg_end+1)
    PIDs = np.zeros(polyDeg_end-polyDeg_ini+1)
    PIDs_std = np.zeros(polyDeg_end-polyDeg_ini+1)
    locPIDs = np.zeros(5)
    
    for polyDeg in polyDegs:
        
        if precompileN: 
            build_dir="build_Euler_"+test_id+"_N"+str(polyDeg)
            exec_path = compileFluxo(precompileFlux, volFlux, build_dir=build_dir, extra_options = extra_options+" -DFLUXO_POLYNOMIAL_DEGREE="+str(polyDeg), compiler=compiler)
            
        # Generate copy of template parameter file
        os.system("cp parameter_template.ini "+prm_path)
        os.system("sed -i 's/@polyDeg@/"+str(polyDeg)+"/g' "+prm_path)
        os.system("sed -i 's/@volFlux@/"+str(volFlux)+"/g' "+prm_path)
        
        # Run fluxo
        cmd = exec_path.strip()+' '+prm_path.strip()
        
        for i in range(5):
            os.system(cmd+' 2>std'+test_id+'.err 1>std'+test_id+'.out')
            
            locPIDs[i] = read_PID('std'+test_id+'.out')
             
        PIDs[polyDeg-polyDeg_ini] = np.mean(locPIDs)
        PIDs_std[polyDeg-polyDeg_ini] = np.std(locPIDs)
        print("  ->  N=%d : %11.5E (+-%11.5E) Time/RHS/DOF[s]" % (polyDeg,PIDs[polyDeg-polyDeg_ini],PIDs_std[polyDeg-polyDeg_ini]))
        
        with open(output_file_name, "a") as output_file:
            output_file.write("%3d , %11.5E , %11.5E \n" % (polyDeg,PIDs[polyDeg-polyDeg_ini],PIDs_std[polyDeg-polyDeg_ini]))
        
        # Delete build directory
        if precompileN: 
            os.system("rm -r ../builds/"+build_dir)
        
        # Delete parameter file
        os.system("rm "+prm_path)
        
    output_file.close()
    
    # Delete build directory
    if not precompileN:
        os.system("rm -r ../builds/"+build_dir)
        
##############
# MAIN PROGRAM
##############


# Check python version
######################
assert (sys.version_info.major == 3),'python>=3.0 is needed!'

# Define arguments
##################
parser = argparse.ArgumentParser(description='Tool to run serial performance tests with FLUXO (Euler equations)', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--compiler', type=str, default='', help='Fortran compiler for FLUXO. If empty, the Fortran compiler found by cmake is used.' )

parser.add_argument('--id', type=str, default='_', help='Test Identifier. Used for the build directories ("build/build_Euler_<id>...") and the results files ("serialTests/results/results_<id>_.dat")' )

parser.add_argument('--flux', type=int, default=32, help='Volume and surface flux used for the tests.\n'
                                                        '    32: [DEFAULT] Ranocha flux\n'
                                                        '    40: Shima et al. flux\n')

parser.add_argument('--precompileN', dest='precompileN', action='store_true', help='[DEFAULT] Compile fluxo with fixed polynomial degree.')
parser.add_argument('--no-precompileN', dest='precompileN', action='store_false', help='[optional] deactivate --precompileN.')
parser.set_defaults(precompileN=True)

parser.add_argument('--precompileFlux', dest='precompileFlux', action='store_true', help='[DEFAULT] Compile fluxo with volume flux.')
parser.add_argument('--no-precompileFlux', dest='precompileFlux', action='store_false', help='[optional] deactivate --precompileN.')
parser.set_defaults(precompileFlux=True)

parser.add_argument('--extra_options', type=str, default="", help='Additional options for cmake. For example, to compile with Gauss nodes use --extra_options=" -DFLUXO_DISC_NODETYPE=GAUSS"')

args = parser.parse_args()


# Run test
##########

RunPerformanceTest(args.id,volFlux=args.flux,compiler=args.compiler, precompileN = args.precompileN, precompileFlux=args.precompileFlux, extra_options=args.extra_options)
