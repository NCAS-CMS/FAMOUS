#This is a general example for a slurm config with both the account and queue
#Account info
#SBATCH --account=n02
#Queue info
#SBATCH --partition=standard
#SBATCH --qos=standard

#Extra options required for srun
#export SRUN_OPTS=""

#Override standard mpif90 command, if required.
#export MPIF90_UM=""

#If shared exec required set MPIF90_SHARED
#export MPIF90_SHARED=1

#Set LD_LIBRARY_PATH if required.
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH: