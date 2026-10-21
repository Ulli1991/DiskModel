#!/bin/bash
# Submit from the repository root:  sbatch scripts/script_rusty.sh
#SBATCH --job-name="make_ics"
#SBATCH --output="out.log"
#####SBATCH --error="err.log"
#SBATCH -p preempt --qos=preempt
#SBATCH --nodes=1
##########SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 168:00:00    
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=uli@mpa-garching.mpg.de
#SBATCH --constraint=rome

module purge

module load slurm/20.02.5
module load gcc/7.4.0
module load openmpi4/4.0.5
module load lib/fftw3/3.3.8-openmpi4 
module load python2/2.7.16  

#./MakeDiskGalaxy params/input_myMW.txt 
#./MakeDiskGalaxy params/input_LMC_highres_MMH.txt 
./MakeDiskGalaxy params/input_dwarf_01solar_mass.txt 
