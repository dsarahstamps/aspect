#!/bin/bash
#SBATCH -J run_g6_traction_uniTopo
#SBATCH -o g6_traction_uniTopo.o%j
#SBATCH -n 64 
#SBATCH -p normal 
#SBATCH -t 28:00:00
ibrun $WORK/aspect/build/aspect g6_traction_vector_uniTopo.prm
