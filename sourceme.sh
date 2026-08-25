if ! type module >/dev/null 2>&1 && [ -r /opt/modules/init/bash ]; then
    source /opt/modules/init/bash
fi

# Module names differ between machines.  The development workstation uses the
# system toolchain and has no Environment Modules installation; the two Slurm
# clusters keep their existing explicit module selections.
if type module >/dev/null 2>&1 && [ -d /lustre/home/2401110128 ]; then
    module purge
    # wm2 default modules verified for this checkout.
    module load cmake/3.31.9
    module load OpenBLAS/0.3.17
    module load gsl/2.7.0
    module load boost/1.83.0

    export GSL_ROOT_DIR=/lustre/software/gsl/2.7.0/gcc_8.5.0
    export CMAKE_PREFIX_PATH="$GSL_ROOT_DIR:$CMAKE_PREFIX_PATH"
elif type module >/dev/null 2>&1 && [ -d /opt/library/modulefiles ]; then
    module use /opt/library/modulefiles
    module purge
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

# A root-level CMake configure places targets in build/src, while the cluster
# checkouts historically configure src/ directly and place them in build/.
# Keep both layouts available without changing which executable job scripts use.
export PYTHONPATH="$PWD/build/src:$PWD/build:$PYTHONPATH"
export LD_LIBRARY_PATH="$PWD/build/src:$PWD/build:$LD_LIBRARY_PATH"
