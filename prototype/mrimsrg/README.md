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
cmake --build prototype/mrimsrg/build --target mrimsrg_prepare mrimsrg_validate -j2
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

`commutator.py` implements the natural-orbital MR-IMSRG(2) commutator with
`lambda3=0`.  Its equations were regenerated independently with QCombo 0.2.0:
for a random four-orbit correlated reference, direct summation of the QCombo
0B/1B/2B expressions agreed with the numerical implementation to
`1.8e-15`.  The committed tests additionally compare the correlated 1B--2B
zero-body term and the single-reference single/double-excitation blocks with
explicit Fock-space matrix commutators. Direct calls reject a non-diagonal
`gamma1` instead of silently treating its diagonal as natural occupations.
For the required `He4/O16, Nrefmax=2` flows, `flow.py` diagonalizes only the
connected off-diagonal blocks of `gamma1`, evaluates the published equations
in that temporary basis, and transforms every generator back before applying
the `Delta e` mask in the original HO basis. Final and checkpoint Hamiltonians
are transformed back to the unchanged input orbit ordering; the exact
orthogonal matrix is saved as `formula_basis_vectors.npy`. This is an internal
basis-covariant evaluation, not a natural-orbital NCSM optimization or a
replacement of the HO quantum-number labels.

`generator.py` implements one generator only: the modified White generator
used in the IM-NCSM literature, restricted by the relaxed `Delta e != 0`
one- and two-body masks. Directional `D1/D2=<Psi|H:A:|Psi>` matrix elements
retain the published terms linear in `lambda2` and set `lambda3=0`; their
anti-Hermitian combination is monitored as the decoupling residual. The
leading Epstein--Nesbet denominators are Mongelli Eqs. (5.213)--(5.214), with
the displayed `O(lambda2)` denominator corrections omitted and a
sign-preserving `1e-6 MeV` cutoff. In the Slater limit they reduce term by
term to the denominator in `src/Generator.cc`. QCombo positive-product
expansions fix the directional index order, while random correlated tests
compare the anti-Hermitian residual directly with
`<Psi|[H,:A:]|Psi>` from the independently verified commutator.

`flow.py` directly integrates `dH/ds=[eta,H]` with an adaptive DOP853
stepper.  Like the existing C++ `IMSRGSolver`, it updates `eta` from every
trial Hamiltonian and stops on the masked `D-D^dagger` residual norm. The ODE state stores only independent
Hermitian one-body and antisymmetric/Hermitian pair-space two-body elements;
this mirrors the established operator storage and prevents redundant tensor
components from developing symmetry-violating numerical modes.  On the first
real He4 diagnostic, the masked residual fell to `4.79e-3` of its initial value
by `s=0.35`, while all reconstructed tensor symmetry errors remained exactly
zero.

Run and materialize one calculation with:

```bash
PYTHONPATH=prototype/mrimsrg python3 prototype/mrimsrg/run_mrimsrg.py \
  prototype/mrimsrg/data/He4_Nrefmax0_final \
  prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --checkpoint-s 1.0
```

The output directory contains `gamma1/gamma2/lambda2`, initial and final
MR-normal-ordered tensors, and the final ordinary vacuum-normal-ordered
`E0+t+V` tensors.  `metadata.json` records the fixed interaction identity,
reference metadata, `lambda3=0`, generator/mask, ODE settings and every
accepted-step residual.  Existing output directories are never overwritten.
The output also materializes the ordinary vacuum Hamiltonian at `s=0` and at
the requested intermediate `--checkpoint-s` under `checkpoints/`; the root
payload is the final point.  Thus a production run retains the three points
required by the acceptance plan without keeping every large ODE state.
The command returns status 2 after saving if `smax` is reached before the
residual target, allowing a Slurm job to distinguish a diagnostic from an
accepted result.

Read a materialized ordinary Hamiltonian back into the existing NCSM solver:

```bash
prototype/mrimsrg/build/mrimsrg_validate \
  --interaction /home/mengziyan/Forces/N2LO_opt/TwBME_N2LO_opt_hw20_emax2_e2max4.minipack \
  --flow-output prototype/mrimsrg/data/He4_flow_rtol1e-6 \
  --Z 2 --N 2 --nmax 8
```

The validator uses the original interaction only to initialize the identical
orbit/channel ordering, replaces its m-scheme matrix elements with the saved
vacuum `E0+t+V`, and then calls the existing `simpleFCI`.  At `s=0` this path
reproduces the complete 3060-dimensional He4 result as
`-20.3388325043 MeV, J=0`.

The bridge requires NumPy and SciPy on the Python side.  The first checked
integration run used `shell-model-obs` revision `1687f16` and reproduced the
full He4 `emax=2` ground-state benchmark as `-20.3388325043 MeV`; contraction
of the independently constructed `gamma1/gamma2` gave the same value.
The dense m-scheme readback interface is revision `4b10afa` of that dependency.

On point7, generate and inspect exactly one Slurm job at a time with:

```bash
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --submit
python3 prototype/mrimsrg/gen_flow_job.py --nucleus He4 --nrefmax 2
```

The dedicated generator is for this Python prototype; it deliberately does
not alter the existing production `gen_job.py` used by `imsrg++`.  It fixes the
documented point7 paths, one node/task, 64 cores, QOS/log conventions, sources
`sourceme.sh`, and prints the complete generated script before an optional
submission.  Reinvoke it for one different nucleus or tolerance at a time;
there is no scan or batch-generation mode.
Use an explicit single-run label such as `--label smax5` when preserving a
short diagnostic and generating a longer flow with otherwise identical
tolerances; an existing materialized output is still never overwritten.

Reference metadata records the original development-machine interaction path.
When the exact fixed file is relocated, pass `--interaction <path>`; the runner
does not trust the path alone and still requires the frozen SHA-256 before it
starts a flow.  The point7 generator uses the task-specific copy under
`/tns/mengziyan/mr-imsrg-inputs/` because the shared force tree has larger
N2LOopt spaces but not this exact emax2 file.
The production path is compatible with point7's system Python 3.9 and the
generated batch script initializes Environment Modules explicitly before it
sources the machine-aware repository environment.
Compute nodes do not carry the login node's `/usr/local` NumPy/SciPy install,
so this standalone prototype calls the shared
`/opt/library/miniconda-3.12.9/bin/python3` interpreter explicitly and logs its
versions.  It does not activate or load miniconda and does not import the
system-Python pyIMSRG module.
