# Rapid MR-IMSRG prototype

This directory contains the deliberately small m-scheme prototype described in
`docs/MR-IMSRG-快速结果计划.md`.  The C++ bridge does not implement another
interaction reader or NCSM solver: it builds against the validated
`/home/mengziyan/shell-model-obs` `Hamiltonian` and `simpleFCI` libraries.

Build the input bridge:

```bash
cmake -S prototype/mrimsrg -B prototype/mrimsrg/build \
  -DSHELL_MODEL_OBS_ROOT=/home/mengziyan/shell-model-obs \
  -DCMAKE_BUILD_TYPE=Release
cmake --build prototype/mrimsrg/build --target mrimsrg_prepare -j2
```

Prepare and check the first He4 reference:

```bash
prototype/mrimsrg/build/mrimsrg_prepare \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --output prototype/mrimsrg/data/He4_Nrefmax0 --Z 2 --N 2 --nrefmax 0
python3 prototype/mrimsrg/densities.py prototype/mrimsrg/data/He4_Nrefmax0
```

Generated Hamiltonians, wavefunctions and densities belong under `data/` and
are ignored by Git.  The bridge refuses to overwrite an existing bundle.

Run the dependency-free unit tests with:

```bash
python3 -m unittest discover -s prototype/mrimsrg/tests -v
```

The bridge requires NumPy and SciPy on the Python side.  The first checked
integration run used `shell-model-obs` revision `1687f16` and reproduced the
full He4 `emax=2` ground-state benchmark as `-20.3388325043 MeV`; contraction
of the independently constructed `gamma1/gamma2` gave the same value.
