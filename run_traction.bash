#!/bin/bash
#SBATCH -J run_g6_CEAR_dyn_vtu
#SBATCH -o g6_CEAR_dyn.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 28:00:00
ibrun $WORK/aspect/build/aspect g6_CEAR.prm
