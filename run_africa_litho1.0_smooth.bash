#!/bin/bash
#SBATCH -J run_africa_litho1.0_smooth_harmAvg.prm
#SBATCH -o africa_litho1.0_smooth_harmAvg_g3.o%j
#SBATCH -n 16 
#SBATCH -p normal 
#SBATCH -t 00:15:00
ibrun $HOME/packages/aspect/build/aspect africa_litho1.0_harmAvg.prm
