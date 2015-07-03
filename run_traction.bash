#!/bin/bash
#SBATCH -J run_g5_CEAR_vtu
#SBATCH -o g5_CEAR.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 08:00:00
ibrun $WORK/aspect/build/aspect g5_CEAR.prm
