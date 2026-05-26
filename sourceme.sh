module purge
module load cmake
module load openblas/0.3.10-single
module load gsl
module load boost
module load miniconda

export LD_LIBRARY_PATH="$PWD/src:$LD_LIBRARY_PATH"
