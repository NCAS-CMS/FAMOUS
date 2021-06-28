#SBATCH --account=n02-cms
#SBATCH --partition=standard
#SBATCH --qos=standard

#Extra options required for srun
export SRUN_OPTS="--cpu-bind=cores"

#Extra stuff to account for $HOME not being available in slurm so
#.profile etc cannot be used in batch.
export DATA_DIR=$WORKDIR
export UMDIR=/work/n02/n02/simon/FAMOUS
. $UMDIR/setvars_4.5

#Override standard mpif90 command for singularity container.
export MPIF90_UM="singularity exec -B /work/n02/n02,/opt/cray,/usr/lib64:/usr/lib/host /work/n02/n02/simon/famous_sing/lfric_env.sif /container/mpich/bin/mpif90"

#Shared exec required so set MPIF90_SHARED
export MPIF90_SHARED=1

#Set LD path to MPICH for singularity built exec.
#Needed due to new module system
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/mpich/8.0.16/ofi/gnu/9.1/lib-abi-mpich
