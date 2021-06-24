#Location of unpacked tarball/git directory.
export UMDIR=/where/the/installation/is/um_famous

#Location of where the model will write its output files to.

export DATA_DIR=${DATA_DIR:-/your/data/directory}

#If umui_runs needs locating away from $HOME
#export UMUI_RUNS=$DATA_DIR

#Only use this line if using own MPI installation
#Note: If using slurm use the config file the slurm dir instead
#export MPI_DIR=/path/to/mpi
#Unless you are using modules. Ensure that the modules are automatically
#loaded when you log on, and uncomment the next line
#export MPI_MODULES=Y

. $UMDIR/setvars_4.5
