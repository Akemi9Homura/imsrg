module purge

# Module names differ between machines, so branch on the machine. wm2 has the
# /lustre/home/2401110128 tree; point7 does not.
if [ -d /lustre/home/2401110128 ]; then
    # wm2 default modules verified for this checkout.
    module load cmake/3.31.9
    module load OpenBLAS/0.3.17
    module load gsl/2.7.0
    module load boost/1.83.0

    export GSL_ROOT_DIR=/lustre/software/gsl/2.7.0/gcc_8.5.0
    export CMAKE_PREFIX_PATH="$GSL_ROOT_DIR:$CMAKE_PREFIX_PATH"
else
    # point7 modules. The checkout's imsrg++ links libopenblas.so.0, which the
    # single-threaded OpenBLAS module provides; gsl/2.7.1 supplies libgsl.so.27.
    module load cmake/3.25.2
    module load openblas/0.3.10-single
    module load gsl/2.7.1
    module load boost/1.81.0
fi

# Do not load miniconda by default: pyIMSRG is built against the system Python,
# and loading miniconda changes python3 and breaks importing the existing module.

export OPENBLAS_NUM_THREADS=1

export PYTHONPATH="$PWD/build:$PYTHONPATH"
export LD_LIBRARY_PATH="$PWD/build:$LD_LIBRARY_PATH"
