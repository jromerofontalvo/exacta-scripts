#!/bin/bash

#SBATCH -n 1                # Number of cores
#SBATCH -N 1                 # Ensure that all cores are on one machine
#SBATCH -t 7-00:00:00           # Runtime in D-HH:MM
#SBATCH -p aspuru-guzik		     # Partition to submit to
#SBATCH --mem=8000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o /n/home08/jromerofontalvo/libor/H2O/test7/outputs/H2O-2.15-libor-LBFGS-exact-0-read-CAS-0-4.out      # File to which STDOUT will be written
#SBATCH -e /n/home08/jromerofontalvo/libor/H2O/test7/outputs/H2O-2.15-libor-LBFGS-exact-0-read-CAS-0-4.err      # File to which STDERR will be written

# User specific aliases and functions
source new-modules.sh
module load Anaconda/2.1.0-fasrc01

cd /n/home08/jromerofontalvo/paper_calculations/

integralFolder=/n/home08/jromerofontalvo/libor/H2O/integrals
guessFolder=/n/home08/jromerofontalvo/libor/H2O/guesses2/guessDMRGGen
initialState=/n/home08/jromerofontalvo/libor/H2O/data/r_2.15/mps.out

python variationalCASlibor.py H2O 2.15-libor 1 "ucc" 0 4 $integralFolder/H2O-2.15-libor.int $guessFolder/H2O-2.15-libor-gen.guess LBFGS exact 0 read 0 11001100 $initialState