#SBATCH --account=<your account>
#SBATCH --partition=standard
#SBATCH --qos=standard


#Extra stuff to account for $HOME not being available in slurm so
#.profile etc cannot be used in batch.
export DATA_DIR=/work/n02/n02/$USER
export UMDIR=/work/n02/n02/simon/FAMOUS
. $UMDIR/setvars_4.5

#Extra options required for srun
export SRUN_OPTS="--cpu-bind=cores"

#Container
export SINGULARITY_CONTAINER=/work/n02/n02/simon/famous_sing/famous.sif
#Or download latest version with
#singularity pull famous.sif library://simonwncas/default/famous.sif:latest 
#and set SINGULARITY_CONTAINER

#If running exec inside container use the following two lines
export SRUN_OPTS="--cpu-bind=cores singularity exec -B /work/n02/n02,/work/n02/shared,/opt/cray,/usr/lib64:/usr/lib/host,/var/spool/slurmd $SINGULARITY_CONTAINER"
export SINGULARITYENV_LOCAL_LD_LIBRARY_PATH=/opt/cray/pe/mpich/8.0.16/ofi/gnu/9.1/lib-abi-mpich:/opt/cray/libfabric/1.11.0.0.233/lib64:/opt/cray/pe/pmi/6.0.7/lib

#Override standard mpif90 command for singularity container.
export MPIF90_UM="singularity exec -B /work/n02/n02,/work/n02/shared,/opt/cray,/usr/lib64:/usr/lib/host $SINGULARITY_CONTAINER /container/mpich/bin/mpif90"

#Shared exec required so set MPIF90_SHARED
export MPIF90_SHARED=1

#Set LD_LIBRARY_PATH to MPICH for singularity built exec.
#Needed due to new module system
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/cray/pe/mpich/8.0.16/ofi/gnu/9.1/lib-abi-mpich
