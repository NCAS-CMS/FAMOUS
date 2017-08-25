#Location of unpacked tarball/git directory.
export UMDIR=/where/the/installation/is/um_famous

#Location of where the model will write its output files to.

export DATA_DIR=${DATA_DIR:-/your/data/directory}

#Only use this line if using own MPI installation
#export MPI_DIR=/path/to/mpi
#Unless you are using modules. Ensure that the modules are automatically
#loaded when you log on, and uncomment the next line
#export MPI_MODULES=Y

. $UMDIR/setvars_4.5
