module purge

# wm2 default modules verified for this checkout.
module load cmake/3.31.9
module load OpenBLAS/0.3.17
module load gsl/2.7.0
module load boost/1.83.0

# Original settings from this checkout before the wm2 adjustment.
# They are kept for reference, but are not used by default on wm2.
# module load cmake
# module load openblas/0.3.10-single
# module load gsl
# module load boost
# module load miniconda
#
# Do not enable miniconda by default here:
# pyIMSRG is built against the system Python unless CMake is told otherwise.
# Loading miniconda here changes python3 and makes the existing module fail to import.

export OPENBLAS_NUM_THREADS=1

export GSL_ROOT_DIR=/lustre/software/gsl/2.7.0/gcc_8.5.0
export CMAKE_PREFIX_PATH="$GSL_ROOT_DIR:$CMAKE_PREFIX_PATH"

export PYTHONPATH="$PWD/build:$PYTHONPATH"
export LD_LIBRARY_PATH="$PWD/build:$LD_LIBRARY_PATH"

# Original path from this checkout before the wm2 adjustment.
# export LD_LIBRARY_PATH="$PWD/src:$LD_LIBRARY_PATH"
