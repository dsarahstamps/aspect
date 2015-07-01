#!/bin/bash
#SBATCH -J run_tractions_g5_wShearStress_vtu
#SBATCH -o g5_tractions_wShearStress_vtu.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 08:00:00
ibrun $WORK/aspect/build/aspect g5_traction_vector.prm
