#!/bin/bash
#SBATCH -J g4_run_africa_litho1.0_crust1.0_smooth_harmAvg_g4
#SBATCH -o africa_africa_litho1.0_crust1.0_smooth_harmAvg_g4.o%j
#SBATCH -n 24 
#SBATCH -p normal 
#SBATCH -t 20:00:00
ibrun $HOME/packages/aspect/build/aspect africa_litho1.0_crust1.0_smooth_harmAvg.prm
