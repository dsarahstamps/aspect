#!/bin/bash
#SBATCH -J run_africa_litho1.0
#SBATCH -o africa_litho1.0_crust1.0_g4.o%j
#SBATCH -n 20 
#SBATCH -p normal 
#SBATCH -t 08:00:00
ibrun $HOME/packages/aspect/build/aspect africa_litho1.0_harmAvg.prm
