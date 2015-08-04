#!/bin/bash
#SBATCH -J run_g5_uniVisc
#SBATCH -o g5_uniVisc.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 28:00:00
ibrun $WORK/aspect/build/aspect g5_traction_uniVisc.prm
